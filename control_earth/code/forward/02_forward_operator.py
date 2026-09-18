"""Phase 2 -- forward operator (METHODOLOGY.md section 4.2).

Truth (spatial, from 01_truth_builder) -> clean truth Stokes (lmax 120,
PREM Love numbers shared with the pipeline) -> degraded synthetic L2 for
three centres x two truncations, plus the synthetic planet's own auxiliary
observation products (TN-14, GravIS geocentre), emitted in the exact
formats the unchanged pipeline parses.

Degradation chain per (centre, truncation, month):
  coeffs = static(real centre mean, D016)
         + truth_stokes(t)          [hydrology + Bagge GIA x dt, D011]
         + noise(t)                 [formal-sigma template x sigma_scale,
                                     white + AR(1)-across-degree stripe
                                     chains per (order, parity), D015]
  C20/C30 <- truth + GRACE-like corruption (AR(1) in time + S2-alias tone)
  degree 0, 1 <- 0                  [L2 contract]

Synthetic auxiliaries:
  TN-14   rows at the REAL TN-14 epochs; C20/C30 = truth + SLR-like noise;
          C30 NaN exactly where real TN-14 has NaN (preserves the pre-2012
          C30 axis degeneracy).
  GravIS  truth degree-1 + small noise, written into a copy of the real
          file's structure.

All amplitudes live in FORWARD_PARAMS (provisional until the Phase-3
reproduce-reality gate; every value registered in decision_register.md
D015/D016).

Run
---
    python \
        code/forward/02_forward_operator.py [--spec ...] [--params ...]
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import pathlib
import sys
import time as _time

import numpy as np
import xarray as xr

PIPELINE_CODE = pathlib.Path(__file__).resolve().parents[3] / "code"
sys.path.insert(0, str(PIPELINE_CODE))

import gravity_toolkit as gt  # noqa: E402
from gravity_toolkit.gen_stokes import gen_stokes  # noqa: E402
from gravity_toolkit.read_GRACE_harmonics import read_GRACE_harmonics  # noqa: E402

from grace_pipeline import config as cfg  # noqa: E402
from grace_pipeline.io_shm import enumerate_months  # noqa: E402
from grace_pipeline.io_aux import load_tn14  # noqa: E402
from grace_pipeline.pipeline_global import GLOBAL_AOI  # noqa: E402

CONTROL_ROOT = pathlib.Path(__file__).resolve().parents[2]
REPO_ROOT = CONTROL_ROOT.parent
TRUTH_DIR = CONTROL_ROOT / "datasets" / "truth"
GIA_DIR = CONTROL_ROOT / "datasets" / "gia_truth"
L2_DIR = CONTROL_ROOT / "datasets" / "synthetic_l2" / "pilot"
AUX_DIR = L2_DIR / "auxiliary"

LMAX_TRUTH = 120
_GT_DATA = pathlib.Path(gt.__file__).parent / "data" / "love_numbers"

# ---------------------------------------------------------------------------
# Provisional degradation amplitudes (D015/D016). The Phase-3 calibration
# gate tunes sigma_scale (and, if needed, stripe_weight/stripe_rho) until the
# synthetic recipe-spread signatures match the real ones; values here are
# starting points, not results.
# ---------------------------------------------------------------------------
FORWARD_PARAMS = {
    "sigma_scale": 12.0,        # formal sigmas are known-optimistic; the
                                # stripe-dominated calibrated error sits an
                                # order of magnitude above formal at high l.
    "stripe_weight": 0.7,       # fraction of noise variance in correlated
                                # (order, parity) chains -> N-S stripes
    "stripe_rho": 0.9,          # AR(1) correlation across degree in a chain
    "c20_grace_sigma": 2.0e-10,  # GRACE-like C20 excess error
    "c20_grace_rho_t": 0.5,     # month-to-month AR(1) of the C20 error
    "c20_alias_amp": 1.0e-10,   # S2-alias tone (161-day) amplitude
    "c30_grace_sigma": 1.5e-10,
    "slr_c20_sigma": 3.0e-11,   # synthetic TN-14 noise (SLR-like, small)
    "slr_c30_sigma": 5.0e-11,
    "geocenter_sigma": 2.0e-11,  # synthetic GravIS degree-1 noise
    "alias_period_days": 161.0,
    "centre_indep_frac": 1.0,   # fraction of stripe-noise VARIANCE that is
                                # centre-specific (D026). 1.0 = fully
                                # independent centres (iterations 1-2). Real
                                # CSR/JPL/GFZ process the same satellite
                                # data, so their errors share a large common
                                # component; the centre axis carries only
                                # ~4% of real cross-recipe variance.
}


def _love(lmax: int):
    return gt.read_love_numbers(str(_GT_DATA), LMAX=lmax, FORMAT="tuple")


def real_sampled_months() -> list:
    ba = {(m.year, m.month): m
          for m in enumerate_months("CSR", 60, months=cfg.APPENDIX_MONTHS)}
    bb = {(m.year, m.month): m
          for m in enumerate_months("CSR", 96, months=cfg.APPENDIX_MONTHS)}
    return [(ba[ym], bb[ym]) for ym in sorted(set(ba) & set(bb))]


def read_gfc_rates(path: pathlib.Path, lmax: int) -> tuple[np.ndarray,
                                                           np.ndarray]:
    """Bagge et al. (2023) ICGEM .gfc -> (clm_rate, slm_rate) per year."""
    c = np.zeros((lmax + 1, lmax + 1))
    s = np.zeros((lmax + 1, lmax + 1))
    for line in path.read_text().splitlines():
        if not line.startswith("gfc"):
            continue
        _, l, m, cv, sv = line.split()[:5]
        l, m = int(l), int(m)
        if l <= lmax:
            c[l, m] = float(cv)
            s[l, m] = float(sv)
    return c, s


# ---------------------------------------------------------------------------
# Noise generator (D015)
# ---------------------------------------------------------------------------

def stripe_mix(lmax: int, rng: np.random.Generator,
               p: dict) -> tuple[np.ndarray, np.ndarray]:
    """Standardized (unit-variance) correlated noise mix.

    White component: independent per coefficient. Stripe component: AR(1)
    chains across degree within each (order, parity) family -- the exact
    correlation structure Swenson destriping and DDK regularisation exist
    to remove, so the filter axis faces a genuine signal-vs-stripe
    trade-off (METHODOLOGY section 4.2). Scaling by the per-centre formal-
    sigma template and sigma_scale happens at the call site, where the
    common and centre-specific mixes are composed (D026).
    """
    w = p["stripe_weight"]
    rho = p["stripe_rho"]
    out = []
    for _ in range(2):                                  # clm, slm
        white = rng.standard_normal((lmax + 1, lmax + 1))
        corr = np.zeros((lmax + 1, lmax + 1))
        for m in range(lmax + 1):
            for parity in (0, 1):
                ls = np.arange(max(m, 2) + ((max(m, 2) + parity) % 2),
                               lmax + 1, 2)
                if ls.size == 0:
                    continue
                z = rng.standard_normal(ls.size)
                chain = np.empty(ls.size)
                chain[0] = z[0]
                for k in range(1, ls.size):
                    chain[k] = rho * chain[k - 1] + \
                        np.sqrt(1 - rho ** 2) * z[k]
                corr[ls, m] = chain
        out.append(np.sqrt(1 - w) * white + np.sqrt(w) * corr)
    return out[0], out[1]


def centre_noise(sigma_c: np.ndarray, sigma_s: np.ndarray, lmax: int,
                 seed: int, centre_idx: int, tr_idx: int, month_idx: int,
                 p: dict) -> tuple[np.ndarray, np.ndarray]:
    """Compose the (centre, truncation, month) noise realization (D026).

    noise = sigma_scale * sig_template * [sqrt(1-f) * common + sqrt(f) * own]

    with f = centre_indep_frac. The common mix is deterministic per
    (truncation, month) and IDENTICAL across centres -- the synthetic
    analogue of CSR/JPL/GFZ sharing the same satellite observations -- so
    inter-centre spread scales as sqrt(f) while each centre's total noise
    variance is unchanged ((1-f) + f = 1). Slot-seeded RNGs keep every
    realization reproducible independent of emission order.
    """
    f = p["centre_indep_frac"]
    scale = p["sigma_scale"]
    rng_common = np.random.default_rng([seed, 1000 + tr_idx, month_idx])
    rng_own = np.random.default_rng(
        [seed, 2000 + 100 * centre_idx + tr_idx, month_idx])
    com_c, com_s = stripe_mix(lmax, rng_common, p)
    own_c, own_s = stripe_mix(lmax, rng_own, p)
    mix_c = np.sqrt(1 - f) * com_c + np.sqrt(f) * own_c
    mix_s = np.sqrt(1 - f) * com_s + np.sqrt(f) * own_s
    return scale * sigma_c * mix_c, scale * sigma_s * mix_s


# ---------------------------------------------------------------------------
# GSM writer (same format as the contract audit)
# ---------------------------------------------------------------------------

_HEADER = """header:
  dimensions:
    degree: {lmax}
    order: {lmax}
  non-standard_attributes:
    product_id: GSM-2
    normalization: fully normalized
    permanent_tide_flag: inclusive
  global_attributes:
    title: Control Earth synthetic GSM ({tag})
    summary: Synthetic Level-2 coefficients for the Control Earth pilot
      planet. Static background from real {centre} monthly means; truth
      hydrology + Bagge et al. (2023) GIA + calibrated noise on top.
      NOT REAL DATA. See control_earth/metadata/decision_register.md.
    institution: control_earth project
# End of YAML header
"""


def write_gsm(path: pathlib.Path, clm, slm, eclm, eslm, lmax: int,
              tag: str, centre: str) -> None:
    parts = [_HEADER.format(lmax=lmax, tag=tag, centre=centre)]
    for l in range(2, lmax + 1):
        for m in range(0, l + 1):
            parts.append(
                f"GRCOF2 {l:4d} {m:4d} {clm[l, m]: .12E} {slm[l, m]: .12E} "
                f"{eclm[l, m]:.4E} {eslm[l, m]:.4E} "
                f"00000000.0000 00000000.0000 nnnn\n")
    path.write_text("".join(parts))


# ---------------------------------------------------------------------------
# Synthetic auxiliaries
# ---------------------------------------------------------------------------

def emit_tn14(truth_c20: dict, truth_c30: dict, months: list,
              rng: np.random.Generator, p: dict) -> pathlib.Path:
    """Synthetic TN-14 at the real TN-14 epochs (D016).

    truth_c20/c30 map (year, month) -> ABSOLUTE truth coefficient. Rows
    for epochs outside the synthetic record keep the real TN-14 values
    (harmless: the pipeline only looks rows up by GSM start epoch).
    """
    real = load_tn14()
    real_lines = cfg.TN14_FILE.read_text().splitlines(keepends=True)
    ym_by_mjd = {}
    for ba, _bb in months:
        ym_by_mjd[round(ba.mjd_start)] = (ba.year, ba.month)

    out, in_data = [], False
    for line in real_lines:
        if not in_data:
            out.append(line)
            if line.startswith("Product:"):
                in_data = True
            continue
        parts = line.split()
        if len(parts) < 10:
            out.append(line)
            continue
        mjd = float(parts[0])
        # nearest synthetic month within the TN-14 matching tolerance
        cand = [ym for k, ym in ym_by_mjd.items() if abs(k - mjd) <= 15]
        if not cand:
            out.append(line)
            continue
        ym = cand[0]
        c20 = truth_c20[ym] + rng.normal(0.0, p["slr_c20_sigma"])
        if parts[5] == "NaN" or ym not in truth_c30:
            c30_txt = "NaN"
        else:
            c30_txt = f"{truth_c30[ym] + rng.normal(0.0, p['slr_c30_sigma']):.10E}"
        out.append(f"{parts[0]} {parts[1]} {c20:.10E} {parts[3]} {parts[4]} "
                   f"{c30_txt} {parts[6]} {parts[7]} {parts[8]} {parts[9]}\n")
    path = AUX_DIR / cfg.TN14_FILE.name
    path.write_text("".join(out))
    return path


def emit_gravis(truth_d1: dict, rng: np.random.Generator,
                p: dict) -> pathlib.Path:
    """Synthetic GravIS: truth degree-1 + small noise, real file structure.

    truth_d1 maps (year, month) -> (c10, c11, s11). Real rows without a
    synthetic month keep their original values.
    """
    real_lines = cfg.GRAVIS_GEOCENTER_FILE.read_text(
        encoding="utf8").splitlines(keepends=True)
    out = []
    for line in real_lines:
        parts = line.split()
        is_data = (len(parts) == 11 and parts[0].replace(".", "").isdigit())
        if not is_data:
            out.append(line)
            continue
        mjd = float(parts[0])
        t = (np.datetime64("1858-11-17") +
             np.timedelta64(int(round(mjd)), "D"))
        ym = (int(str(t.astype("datetime64[Y]"))),
              int(str(t.astype("datetime64[M]")).split("-")[1]))
        if ym not in truth_d1:
            out.append(line)
            continue
        c10, c11, s11 = truth_d1[ym]
        n = rng.normal(0.0, p["geocenter_sigma"], size=3)
        out.append(
            f"{parts[0]} {parts[1]} "
            f"{c10 + n[0]: .10E} {parts[3]} {parts[4]} "
            f"{c11 + n[1]: .10E} {parts[6]} {parts[7]} "
            f"{s11 + n[2]: .10E} {parts[9]} {parts[10]}\n")
    path = AUX_DIR / cfg.GRAVIS_GEOCENTER_FILE.name
    path.write_text("".join(out))
    return path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", type=pathlib.Path,
                    default=TRUTH_DIR / "truth_spec_pilot.json")
    ap.add_argument("--params", type=pathlib.Path, default=None,
                    help="JSON overriding FORWARD_PARAMS (calibration loop)")
    args = ap.parse_args()
    spec = json.loads(args.spec.read_text())
    planet = spec["planet_id"]
    # L2 tree per planet (D030): synthetic_l2/<planet_id>. The module-level
    # constants are rebound so emit_tn14/emit_gravis write into the same tree.
    global L2_DIR, AUX_DIR
    L2_DIR = CONTROL_ROOT / "datasets" / "synthetic_l2" / planet
    AUX_DIR = L2_DIR / "auxiliary"
    p = dict(FORWARD_PARAMS)
    if args.params is not None:
        p.update(json.loads(args.params.read_text()))
    rng = np.random.default_rng(int(spec["seed"]))

    months = real_sampled_months()
    n_t = len(months)
    aoi = GLOBAL_AOI

    # ---- truth spatial -> truth Stokes (lmax 120) --------------------------
    ds = xr.open_dataset(TRUTH_DIR / f"truth_tws_{planet}.nc")
    truth = ds["tws_anom"].to_numpy()                     # (t, lat, lon)
    assert truth.shape[0] == n_t, "truth stack / month list mismatch"

    gia_c, gia_s = read_gfc_rates(GIA_DIR / spec["gia_truth_gfc"], LMAX_TRUTH)
    t_ref = np.datetime64(spec["trend_ref_epoch"])

    love = _love(LMAX_TRUTH)
    print(f"[forward {planet}] SH analysis of {n_t} months at "
          f"lmax {LMAX_TRUTH}")
    t0 = _time.time()
    truth_clm = np.zeros((n_t, LMAX_TRUTH + 1, LMAX_TRUTH + 1))
    truth_slm = np.zeros_like(truth_clm)
    plm = None
    for i in range(n_t):
        h = gen_stokes(truth[i].T, aoi.lon, aoi.lat, LMIN=0,
                       LMAX=LMAX_TRUTH, UNITS=1, PLM=plm, LOVE=love)
        if plm is None and hasattr(h, "PLM"):
            plm = h.PLM  # reuse Legendre polynomials if exposed
        dt_yr = float((months[i][0].time_mid - t_ref)
                      / np.timedelta64(1, "D")) / 365.25
        truth_clm[i] = h.clm + gia_c * dt_yr
        truth_slm[i] = h.slm + gia_s * dt_yr
        if (i + 1) % 50 == 0:
            print(f"  SH {i+1}/{n_t} ({(_time.time()-t0)/60:.1f} min)")
    print(f"[forward {planet}] SH analysis done "
          f"({(_time.time()-t0)/60:.1f} min)")

    # ---- static backgrounds + sigma templates (real record) ---------------
    # Pass over the real files once per (centre, truncation): accumulate the
    # static mean field and keep per-month formal sigmas for the noise draw.
    if L2_DIR.exists():
        import shutil
        shutil.rmtree(L2_DIR)
    AUX_DIR.mkdir(parents=True)
    for c in cfg.CENTRES:
        (L2_DIR / c.lower()).mkdir()

    # C20/C30 GRACE-like corruption: one AR(1)-in-time series per centre
    # (the S2-alias tone is common across centres, phase-locked to time).
    t_mid_days = np.array([(m[0].time_mid - np.datetime64("2002-01-01"))
                           / np.timedelta64(1, "D") for m in months])
    alias = np.sin(2 * np.pi * t_mid_days / p["alias_period_days"])
    c20_err = {}
    c30_err = {}
    for c in cfg.CENTRES:
        e20 = np.empty(n_t)
        e30 = np.empty(n_t)
        e20[0] = rng.normal(0, p["c20_grace_sigma"])
        e30[0] = rng.normal(0, p["c30_grace_sigma"])
        for i in range(1, n_t):
            r = p["c20_grace_rho_t"]
            e20[i] = r * e20[i - 1] + np.sqrt(1 - r ** 2) * \
                rng.normal(0, p["c20_grace_sigma"])
            e30[i] = r * e30[i - 1] + np.sqrt(1 - r ** 2) * \
                rng.normal(0, p["c30_grace_sigma"])
        c20_err[c] = e20 + p["c20_alias_amp"] * alias
        c30_err[c] = e30 + 0.5 * p["c20_alias_amp"] * alias

    truth_c20_abs, truth_c30_abs, truth_d1 = {}, {}, {}
    static_c20_anchor = None

    for centre in cfg.CENTRES:
        token = cfg.CENTRE_TOKEN[centre]
        for tr_idx, (lmax, mnem) in enumerate(((60, "BA01"), (96, "BB01"))):
            real_months = {(m.year, m.month): m for m in enumerate_months(
                centre, lmax, months=cfg.APPENDIX_MONTHS)}
            # static mean + sigma templates over the synthetic month set;
            # months the real centre record lacks reuse the nearest real
            # file's sigmas (registered as part of D014).
            static_c = np.zeros((lmax + 1, lmax + 1))
            static_s = np.zeros((lmax + 1, lmax + 1))
            sig_c, sig_s, have = {}, {}, []
            for ym, m in real_months.items():
                if ym not in {(a.year, a.month) for a, _b in months}:
                    continue
                d = read_GRACE_harmonics(m.file_path, lmax)
                static_c += d["clm"]
                static_s += d["slm"]
                sig_c[ym], sig_s[ym] = d["eclm"], d["eslm"]
                have.append(ym)
            static_c /= len(have)
            static_s /= len(have)
            print(f"[forward {planet}] {centre} {mnem}: static from "
                  f"{len(have)} real months "
                  f"({(_time.time()-t0)/60:.1f} min)")
            if centre == "CSR" and lmax == 60:
                static_c20_anchor = float(static_c[2, 0])
                static_c30_anchor = float(static_c[3, 0])

            for i, (ba, bb) in enumerate(months):
                m = (ba, bb)[tr_idx]
                ym = (m.year, m.month)
                # nearest sigma template if this centre lacks the month
                if ym in sig_c:
                    ec, es = sig_c[ym], sig_s[ym]
                else:
                    k = min(have, key=lambda a: abs(
                        (a[0] - ym[0]) * 12 + a[1] - ym[1]))
                    ec, es = sig_c[k], sig_s[k]
                nc, ns = centre_noise(ec, es, lmax, int(spec["seed"]),
                                      cfg.CENTRES.index(centre), tr_idx,
                                      i, p)
                clm = static_c + truth_clm[i, :lmax + 1, :lmax + 1] + nc
                slm = static_s + truth_slm[i, :lmax + 1, :lmax + 1] + ns
                # low-degree truth bookkeeping (CSR BA01 anchors TN-14/GravIS)
                if centre == "CSR" and lmax == 60:
                    truth_c20_abs[ym] = static_c20_anchor + truth_clm[i, 2, 0]
                    truth_c30_abs[ym] = static_c30_anchor + truth_clm[i, 3, 0]
                    truth_d1[ym] = (truth_clm[i, 1, 0], truth_clm[i, 1, 1],
                                    truth_slm[i, 1, 1])
                # GRACE-like corrupted low degrees (truth + centre error)
                clm[2, 0] = static_c[2, 0] + truth_clm[i, 2, 0] + \
                    c20_err[centre][i]
                clm[3, 0] = static_c[3, 0] + truth_clm[i, 3, 0] + \
                    c30_err[centre][i]
                # L2 contract: degrees 0-1 zero
                clm[0:2, :] = 0.0
                slm[0:2, :] = 0.0
                # The month objects come from the CSR record (D014), so their
                # filenames carry the CSR token UTCSR. Substitute this centre's
                # token so the pipeline's centre filter (io_shm token match)
                # and centre-specific destriping aux resolve correctly
                # (D023). CSR -> UTCSR is a no-op.
                name = m.file_path.name.replace("UTCSR", token)
                write_gsm(L2_DIR / centre.lower() / name, clm, slm, ec, es,
                          lmax, tag=f"{planet} {ym[0]}-{ym[1]:02d} {mnem}",
                          centre=centre)
            print(f"[forward {planet}] {centre} {mnem}: {n_t} GSM emitted "
                  f"({(_time.time()-t0)/60:.1f} min)")

    # ---- synthetic auxiliaries ---------------------------------------------
    tn14_path = emit_tn14(truth_c20_abs, truth_c30_abs, months, rng, p)
    gravis_path = emit_gravis(truth_d1, rng, p)
    print(f"[forward {planet}] auxiliaries: {tn14_path.name}, "
          f"{gravis_path.name}")

    # ---- manifest -----------------------------------------------------------
    files = sorted(L2_DIR.rglob("GSM-2_*"))
    manifest = {
        "planet_id": planet,
        "params": p,
        "n_gsm_files": len(files),
        "n_months": n_t,
        "gia_truth": spec["gia_truth_gfc"],
        "tn14_sha256": hashlib.sha256(tn14_path.read_bytes()).hexdigest(),
        "gravis_sha256": hashlib.sha256(gravis_path.read_bytes()).hexdigest(),
        "created": dt.datetime.utcnow().isoformat() + "Z"}
    (L2_DIR / "forward_manifest.json").write_text(
        json.dumps(manifest, indent=2))
    print(f"[forward {planet}] DONE: {len(files)} GSM files, "
          f"{(_time.time()-t0)/60:.1f} min total")
    return 0


if __name__ == "__main__":
    sys.exit(main())
