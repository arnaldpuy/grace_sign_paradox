"""Phase 0 -- contract audit (METHODOLOGY.md section 4.3).

Proves that a trivial synthetic L2 record round-trips through the UNCHANGED
GRACE pipeline: emitted by a minimal forward operator, ingested by
`io_shm.py`, processed by `pipeline_global.py`, and recovered with the trend
we planted, for the null recipe.

Interpretation registered as D009: because every recipe applies a smoothing
filter (there is no "no filter" level), exact recovery of the raw planted
trend is impossible by construction. The audit criterion is therefore
    recovered(null recipe) == trend of the FILTERED truth
where "filtered truth" is computed with the pipeline's own `apply_filter` +
`synthesise_to_aoi` on the planted Stokes rate. This tests the file format,
the month enumeration, the anomaly logic, the axis identities of the null
levels and the trend arithmetic -- everything except the physics, which the
audit deliberately keeps trivial (pure linear trend, no noise).

Stages
------
A. Internal check: gen_stokes(blob) -> harmonic_summation round-trip sanity.
B. Emit synthetic GSM files (CSR, BA01 + BB01) for the first 24 appendix
   months: C(t) = C_static(real 2002-04 CSR GSM) + R_lm * dt_years.
C. Invoke the unchanged pipeline via runtime path redirection (D005),
   filters=("G300",), window="appendix", max_months=24.
D. Compare recovered null-recipe per-cell trends against the predicted
   filtered-truth trend field. PASS iff max|diff| < TOL_CM_YR.

Run
---
    cd <repository>/control_earth
    python \
        code/forward/00_contract_audit.py
"""
from __future__ import annotations

import json
import pathlib
import shutil
import sys
import time as _time

import numpy as np

# --- pipeline package (READ-ONLY source; imported, never modified) ----------
PIPELINE_CODE = pathlib.Path(__file__).resolve().parents[3] / "code"
sys.path.insert(0, str(PIPELINE_CODE))

import gravity_toolkit as gt  # noqa: E402
from gravity_toolkit.gen_stokes import gen_stokes  # noqa: E402

from grace_pipeline import config as cfg  # noqa: E402
from grace_pipeline.filters_sh import apply_filter  # noqa: E402
from grace_pipeline.io_shm import load_shm  # noqa: E402
from grace_pipeline.synthesis import synthesise_to_aoi, _cmwe_factor  # noqa: E402
from grace_pipeline.pipeline_global import (  # noqa: E402
    GLOBAL_AOI, compute_centre_global)

CONTROL_ROOT = pathlib.Path(__file__).resolve().parents[2]
REPO_ROOT = CONTROL_ROOT.parent
AUDIT_L2 = CONTROL_ROOT / "datasets" / "synthetic_l2" / "audit"
AUDIT_OUT = CONTROL_ROOT / "datasets" / "output" / "audit_cube"
RESULT_JSON = CONTROL_ROOT / "datasets" / "output" / "audit_result.json"

N_MONTHS = 24                      # first 24 appendix months (2002-04 ...)
RATE_PEAK_CM_YR = 5.0              # planted peak trend, cm EWH / yr
TOL_CM_YR = 1e-3                   # pass tolerance on |recovered - predicted|
BLOB_LAT, BLOB_LON, BLOB_SIGMA_DEG = 27.0, 79.0, 8.0   # IGP-ish Gaussian blob

_GT_DATA = pathlib.Path(gt.__file__).parent / "data" / "love_numbers"


def _love(lmax: int):
    return gt.read_love_numbers(str(_GT_DATA), LMAX=lmax, FORMAT="tuple")


# ---------------------------------------------------------------------------
# Stage A -- truth rate pattern and analysis/synthesis sanity
# ---------------------------------------------------------------------------

def build_rate_harmonics(lmax: int):
    """Gaussian blob (peak RATE_PEAK_CM_YR cm/yr) -> Stokes rate R_lm.

    The blob is only a generator of a plausible band-limited pattern; the
    SH-domain object returned here IS the truth rate (D009), so no
    band-limitation error enters the audit.
    """
    # gen_stokes expects data (nlon, nlat) on cell-centre lon/lat vectors.
    aoi = GLOBAL_AOI
    lon_g, lat_g = np.meshgrid(aoi.lon, aoi.lat)           # (nlat, nlon)
    d2 = (lat_g - BLOB_LAT) ** 2 + (np.minimum(
        np.abs(lon_g - BLOB_LON), 360 - np.abs(lon_g - BLOB_LON))) ** 2
    blob = RATE_PEAK_CM_YR * np.exp(-d2 / (2 * BLOB_SIGMA_DEG ** 2))
    harm = gen_stokes(blob.T, aoi.lon, aoi.lat, LMIN=0, LMAX=lmax,
                      UNITS=1, LOVE=_love(lmax))
    # Zero degrees 0-1: the L2 contract excludes them (real GSM files start
    # at degree 2; degree-1 lives in the geocentre product). Registered in
    # D007/D009.
    harm.clm[0:2, :] = 0.0
    harm.slm[0:2, :] = 0.0
    return harm, blob


def stage_a_sanity(lmax: int = 60) -> dict:
    harm, blob = build_rate_harmonics(lmax)
    back = synthesise_to_aoi(harm, GLOBAL_AOI, lmax)       # (nlat, nlon), cm
    aoi = GLOBAL_AOI
    # This stage guards against CONVENTION errors (transposes, lon wrapping),
    # not truncation physics: global-domain correlation is diluted by Gibbs
    # ringing far from the blob (measured 0.969 at lmax 60), so score (i) the
    # recovered peak lands within one cell of the planted centre and (ii) the
    # near-field (25 deg) correlation is ~1. A transpose/wrap error fails both
    # catastrophically.
    i, j = np.unravel_index(np.nanargmax(back), back.shape)
    peak_dist_deg = float(max(abs(aoi.lat[i] - BLOB_LAT),
                              abs(aoi.lon[j] - BLOB_LON)))
    lon_g, lat_g = np.meshgrid(aoi.lon, aoi.lat)
    near = (np.abs(lat_g - BLOB_LAT) < 25) & (np.abs(lon_g - BLOB_LON) < 25)
    corr = float(np.corrcoef(back[near], blob[near])[0, 1])
    peak_ratio = float(np.nanmax(back) / np.nanmax(blob))
    ok = bool((corr > 0.999) and (0.75 < peak_ratio < 1.25)
              and (peak_dist_deg <= aoi.res))
    return {"stage": "A", "local_corr": corr, "peak_ratio": peak_ratio,
            "peak_dist_deg": peak_dist_deg, "pass": ok}


# ---------------------------------------------------------------------------
# Stage B -- emit synthetic GSM files
# ---------------------------------------------------------------------------

_HEADER_TEMPLATE = """header:
  dimensions:
    degree: {lmax}
    order: {lmax}
  non-standard_attributes:
    product_id: GSM-2
    normalization: fully normalized
    permanent_tide_flag: inclusive
  global_attributes:
    title: Control Earth synthetic GSM ({tag})
    summary: Synthetic Level-2 coefficients emitted by the Control Earth
      forward operator (contract audit stage). Static background from the
      real CSR {static_name}; planted Stokes rate on top. NOT REAL DATA.
    institution: control_earth project
# End of YAML header
"""


def real_sampled_months(n_target: int) -> list:
    """The synthetic record must mirror the REAL mission sampling.

    First audit run failed on synthetic 2002-05: real GRACE has no May-2002
    solution, so TN-14 carries no SLR row near that epoch and the pipeline
    (correctly) refuses the month. The synthetic L2 therefore reuses the
    real CSR record's months AND filename date-spans (short arcs included),
    which is also what METHODOLOGY section 4.2 prescribes ("impose the
    actual monthly sampling"). Registered as D010.

    Returns real `Month` objects common to BA01 and BB01, first n_target.
    """
    from grace_pipeline.io_shm import enumerate_months
    ba = {(m.year, m.month): m
          for m in enumerate_months("CSR", 60, months=cfg.APPENDIX_MONTHS)}
    bb = {(m.year, m.month): m
          for m in enumerate_months("CSR", 96, months=cfg.APPENDIX_MONTHS)}
    common = sorted(set(ba) & set(bb))[:n_target]
    return [(ba[ym], bb[ym]) for ym in common]


def write_gsm(path: pathlib.Path, clm: np.ndarray, slm: np.ndarray,
              lmax: int, tag: str, static_name: str) -> None:
    lines = [_HEADER_TEMPLATE.format(lmax=lmax, tag=tag,
                                     static_name=static_name)]
    for l in range(2, lmax + 1):
        for m in range(0, l + 1):
            lines.append(
                f"GRCOF2 {l:4d} {m:4d} {clm[l, m]: .12E} {slm[l, m]: .12E} "
                f"0.0000E+00 0.0000E+00 00000000.0000 00000000.0000 nnnn\n")
    path.write_text("".join(lines))


def stage_b_emit() -> list:
    month_pairs = real_sampled_months(N_MONTHS)
    # Reference epoch for the planted linear trend: mid of the emitted window
    # (trend recovery is invariant to this choice; it only sets the anomaly
    # zero-point).
    t_mid = np.datetime64("2003-04-15")

    real_csr = cfg.CENTRE_DIRS["CSR"]      # capture BEFORE redirection
    csr_dir = AUDIT_L2 / "csr"
    if AUDIT_L2.exists():
        shutil.rmtree(AUDIT_L2)
    csr_dir.mkdir(parents=True)
    # Empty jpl/gfz dirs so accidental centre runs fail loudly, not wrongly.
    (AUDIT_L2 / "jpl").mkdir()
    (AUDIT_L2 / "gfz").mkdir()

    for tr_idx, (lmax, mnem) in enumerate(((60, "BA01"), (96, "BB01"))):
        static_name = f"GSM-2_2002095-2002120_GRAC_UTCSR_{mnem}_0600"
        static = load_shm(real_csr / static_name, lmax=lmax)
        rate, _ = build_rate_harmonics(lmax)
        # gen_stokes returned Stokes for a cm-valued field; the blob is cm/yr,
        # so R_lm is Stokes per year by construction.
        for pair in month_pairs:
            m = pair[tr_idx]
            dt_yr = float((m.time_mid - t_mid)
                          / np.timedelta64(1, "D")) / 365.25
            clm = static.clm + rate.clm * dt_yr
            slm = static.slm + rate.slm * dt_yr
            # Reuse the REAL file's name (same date span, short arcs incl.),
            # so month enumeration, midpoints and TN-14 lookups behave
            # exactly as on the real record (D010).
            write_gsm(csr_dir / m.file_path.name, clm, slm, lmax,
                      tag=f"audit {m.year}-{m.month:02d} {mnem}",
                      static_name=static_name)
    n = len(list(csr_dir.iterdir()))
    print(f"[stage B] emitted {n} synthetic GSM files "
          f"({len(month_pairs)} months x 2 truncations) -> {csr_dir}")
    return month_pairs


# ---------------------------------------------------------------------------
# Stage C -- run the unchanged pipeline on the synthetic record
# ---------------------------------------------------------------------------

def stage_c_recover(month_pairs: list) -> pathlib.Path:
    # D005: runtime path redirection only; zero pipeline-code changes.
    cfg.CENTRE_DIRS["CSR"] = AUDIT_L2 / "csr"
    AUDIT_OUT.mkdir(parents=True, exist_ok=True)
    # max_months slices the appendix window by POSITION; real months are
    # sparse in it, so cover up to the slot of the last emitted month.
    last_ym = (month_pairs[-1][0].year, month_pairs[-1][0].month)
    n_slots = cfg.APPENDIX_MONTHS.index(last_ym) + 1
    t0 = _time.time()
    out = compute_centre_global("CSR", filters=("G300",), window="appendix",
                                out_dir=AUDIT_OUT, max_months=n_slots)
    print(f"[stage C] pipeline done in {(_time.time()-t0)/60:.1f} min -> {out}")
    return out


# ---------------------------------------------------------------------------
# Stage D -- score the null recipe against the filtered-truth prediction
# ---------------------------------------------------------------------------

def stage_d_score(cube_path: pathlib.Path) -> dict:
    import netCDF4 as nc4
    nc = nc4.Dataset(cube_path, "r")
    filt_names = [str(s) for s in nc.variables["filter"][:]]
    gia_names = [str(s) for s in nc.variables["gia_model"][:]]
    c20_names = [str(s) for s in nc.variables["c20_treatment"][:]]
    c30_names = [str(s) for s in nc.variables["c30_treatment"][:]]
    gc_names = [str(s) for s in nc.variables["geocenter"][:]]
    i_f = filt_names.index("G300")
    i_g = gia_names.index("none")
    i_c20 = c20_names.index("original")
    i_c30 = c30_names.index("original")
    i_gc = gc_names.index("none")

    times = nc.variables["time"][:]
    t_yr = (times - times[0]) / 365.25
    results = {}
    for tr_idx, lmax in enumerate(cfg.TRUNCATIONS):
        y = nc.variables["ewh_anom"][:, tr_idx, i_f, i_g, i_c20, i_c30, i_gc,
                                     :, :]
        y = np.asarray(np.ma.filled(y, np.nan))             # (t, lat, lon)
        # Drop the all-NaN slots (appendix window positions with no synthetic
        # month — real-mission sampling gaps, D010).
        have = np.isfinite(y).any(axis=(1, 2))
        y, tv = y[have], t_yr[have]
        tm = tv - tv.mean()
        ym = y - y.mean(axis=0, keepdims=True)
        slope = (tm[:, None, None] * ym).sum(0) / (tm ** 2).sum()

        rate, _ = build_rate_harmonics(lmax)
        pred = synthesise_to_aoi(apply_filter(rate, "G300", lmax),
                                 GLOBAL_AOI, lmax)          # cm/yr field
        diff = slope - pred
        results[f"BA01_{lmax}" if lmax == 60 else f"BB01_{lmax}"] = {
            "max_abs_diff_cm_yr": float(np.nanmax(np.abs(diff))),
            "rms_diff_cm_yr": float(np.sqrt(np.nanmean(diff ** 2))),
            "peak_recovered": float(np.nanmax(slope)),
            "peak_predicted": float(np.nanmax(pred)),
        }
    nc.close()

    ok = all(v["max_abs_diff_cm_yr"] < TOL_CM_YR for v in results.values())
    return {"stage": "D", "tolerance_cm_yr": TOL_CM_YR,
            "per_truncation": results, "pass": bool(ok)}


def main() -> int:
    report = {"n_months": N_MONTHS, "rate_peak_cm_yr": RATE_PEAK_CM_YR}
    a = stage_a_sanity()
    print(f"[stage A] local_corr={a['local_corr']:.6f} "
          f"peak_ratio={a['peak_ratio']:.4f} "
          f"peak_dist={a['peak_dist_deg']:.2f}deg pass={a['pass']}")
    report["stage_a"] = a
    if not a["pass"]:
        RESULT_JSON.write_text(json.dumps(report, indent=2))
        print("STAGE A FAILED -- convention error; stopping.")
        return 1

    month_pairs = stage_b_emit()
    cube = stage_c_recover(month_pairs)
    d = stage_d_score(cube)
    report["stage_d"] = d
    RESULT_JSON.parent.mkdir(parents=True, exist_ok=True)
    RESULT_JSON.write_text(json.dumps(report, indent=2))
    print(json.dumps(d, indent=2))
    print("CONTRACT AUDIT:", "PASS" if d["pass"] else "FAIL")
    return 0 if d["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
