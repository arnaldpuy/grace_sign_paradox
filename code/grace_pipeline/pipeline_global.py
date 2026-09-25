"""Global Tier-A pipeline: streaming variant-matrix synthesis at 0.5 deg.

A fork of pipeline.py that produces the per-centre EWH cube on a global
0.5 deg grid (340 x 720 cells), with three architectural changes required
at this scale:

1.  **Streaming writes.** The regional pipeline allocates the full cube
    in RAM up front (~1 GB for IGP at 0.25 deg). Global at 0.5 deg would
    require ~111 GB in RAM -- too big. We open the output NetCDF once
    and write each month's slice as soon as it's computed.

2.  **float32 + zlib L4 storage of ewh_anom only.** The regional cube
    stores both ewh_abs (full synthesised EWH ~1e9 cm) and ewh_anom
    (anomaly ~±25 cm). At float32 the ewh_abs magnitude floors precision
    at ~100 cm, larger than the inter-recipe hydrology variance we want
    to study. So Tier-A stores ewh_anom only, at float32 (precision
    ~1e-6 cm relative to ±25 cm anomaly range). 4x compression of the
    smooth field gives ~43 GB total for 3 centres x POC.

3.  **Two-pass synthesis.** Computing ewh_anom requires the time mean
    per (truncation, filter, gia, c20, c30, geocenter, lat, lon),
    which is only known after all months are processed. So Pass 1
    synthesises + accumulates per-recipe-cell sum/count; Pass 2
    re-synthesises and writes (ewh - mean) as float32. Doubles wall
    clock vs. a single-pass alternative; in exchange we get
    precision-safe anomalies with no temp disk needed.

The regional `pipeline.py` is left untouched (zero regression risk for
the IGP/CV/HP case-study cubes).

Usage
-----
    cd <repository>/code
    python -m grace_pipeline.pipeline_global --centre CSR
    # Smoke test on 2 months:
    python -m grace_pipeline.pipeline_global --centre CSR --max_months 2

Output: data/processed/grace/global/grace_ewh_global_<centre>.nc
"""
from __future__ import annotations

import argparse
import datetime as dt
import subprocess
import sys
import time as _time
from pathlib import Path
from typing import Iterable

import netCDF4 as nc4
import numpy as np
import gravity_toolkit as gt

from . import config as cfg
from .corrections_grid import subtract_gia_at_epoch
from .corrections_sh import replace_low_degrees
from .filters_sh import apply_filter
from .geocenter import apply_geocenter, load_geocenter
from .io_aux import load_gia_rate, load_tn14
from .io_shm import Month, enumerate_months, load_shm
from .synthesis import synthesise_to_aoi


GLOBAL_AOI = cfg.AOI(lat_min=-85.0, lat_max=85.0,
                     lon_min=-180.0, lon_max=180.0, res=0.5)
GLOBAL_OUT_DIR = cfg.OUT_DIR / "global"


def _git_commit() -> str:
    try:
        out = subprocess.check_output(
            ["git", "-C", str(cfg.REPO_ROOT),
             "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
        )
        return out.decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return "unknown"


def _create_skeleton(path: Path, centre: str, n_t: int, aoi: cfg.AOI,
                      times: np.ndarray, ran_filters: Iterable[str],
                      window: str = "poc") -> nc4.Dataset:
    """Create the output NetCDF with all dims + the (empty) ewh_anom variable.

    `window` is recorded as a global attribute and (when not "poc") also
    triggers the appendix-window caveats attribute (mission gap, late-mission
    GRACE noise, C30 axis collapse pre-MJD 55987)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    nc = nc4.Dataset(path, mode="w", format="NETCDF4")

    nc.createDimension("time", n_t)
    nc.createDimension("truncation", len(cfg.TRUNCATIONS))
    nc.createDimension("filter", len(cfg.FILTERS))
    nc.createDimension("gia_model", len(cfg.GIA_MODELS))
    nc.createDimension("c20_treatment", len(cfg.C20_TREATMENTS))
    nc.createDimension("c30_treatment", len(cfg.C30_TREATMENTS))
    nc.createDimension("geocenter", len(cfg.GEOCENTER_TREATMENTS))
    nc.createDimension("lat", aoi.lat.size)
    nc.createDimension("lon", aoi.lon.size)

    v_time = nc.createVariable("time", "f8", ("time",))
    v_time.units = "days since 1970-01-01"
    v_time[:] = (times.astype("datetime64[D]").astype("int64") -
                  np.datetime64("1970-01-01").astype("datetime64[D]").astype("int64"))

    v_lat = nc.createVariable("lat", "f4", ("lat",))
    v_lat[:] = aoi.lat.astype(np.float32)
    v_lat.units = "degrees_north"
    v_lon = nc.createVariable("lon", "f4", ("lon",))
    v_lon[:] = aoi.lon.astype(np.float32)
    v_lon.units = "degrees_east"

    for dim, values in [
        ("truncation",   list(cfg.TRUNCATIONS)),
        ("filter",       list(cfg.FILTERS)),
        ("gia_model",    list(cfg.GIA_MODELS)),
        ("c20_treatment", list(cfg.C20_TREATMENTS)),
        ("c30_treatment", list(cfg.C30_TREATMENTS)),
        ("geocenter",    list(cfg.GEOCENTER_TREATMENTS)),
    ]:
        v = nc.createVariable(dim, str, (dim,))
        for i, val in enumerate(values):
            v[i] = str(val)

    chunks = (1, 1,
              len(cfg.FILTERS), len(cfg.GIA_MODELS),
              len(cfg.C20_TREATMENTS), len(cfg.C30_TREATMENTS),
              len(cfg.GEOCENTER_TREATMENTS),
              aoi.lat.size, aoi.lon.size)
    ewh = nc.createVariable(
        "ewh_anom", "f4",
        ("time", "truncation", "filter", "gia_model",
         "c20_treatment", "c30_treatment", "geocenter", "lat", "lon"),
        zlib=True, complevel=4, fill_value=np.float32("nan"),
        chunksizes=chunks)
    ewh.long_name = "equivalent water height anomaly"
    ewh.units = "cm"
    ewh.comment = (
        "Tier-A global 0.5 deg cube. ewh_abs is NOT stored on disk: at "
        "float32 the ~1e9 cm absolute-EWH magnitude floors precision at "
        "~100 cm, which is larger than the inter-recipe hydrology "
        "variance we study. ewh_anom is the absolute EWH minus its "
        "time-mean over the reference window, taken per "
        "(truncation, filter, gia_model, c20_treatment, c30_treatment, "
        "geocenter, lat, lon). Anomaly values are O(±25 cm); float32 "
        "precision is ~1e-6 cm, fully sufficient for the 1e-2 cm signal "
        "level. Reference window = full POC."
    )

    mission_label = ("GRACE-FO" if window == "poc"
                     else "GRACE + GRACE-FO")
    nc.title = f"Tier-A {mission_label} EWH anomaly cube (global) -- {centre}"
    nc.centre = centre
    nc.aoi = f"global 0.5 deg, lat [{aoi.lat_min}, {aoi.lat_max}], lon [-180, 180]"
    nc.window = window
    nc.gravity_toolkit_version = gt.__version__
    nc.code_git_commit = _git_commit()
    nc.processing_date = dt.datetime.utcnow().isoformat(timespec="seconds") + "Z"
    nc.tn14_file = cfg.TN14_FILE.name
    nc.gia_files = ",".join(p.name for p in cfg.GIA_FILES.values())
    nc.filters_attempted = ",".join(cfg.FILTERS)
    nc.filters_realised = ",".join(ran_filters)
    nc.Conventions = "CF-1.10"
    nc.reference_period_note = (
        f"ewh_anom = ewh_abs - ewh_abs.mean('time'); reference period is "
        f"the full {window} window. To rebase against a different reference "
        "period, the synthesis must be re-run (ewh_abs is not stored on "
        "disk for precision reasons; see ewh_anom comment).")
    if window == "appendix":
        nc.appendix_caveats = (
            "Full-mission (2002-04 -> 2025-12) window. Three known wrinkles, "
            "all handled cleanly downstream by the NaN-aware OLS / metric "
            "reductions: "
            "(1) C30 axis is degenerate pre-MJD 55987 (~April 2012) because "
            "TN-14 has no C30 over this period; corrections_sh.py silently "
            "retains the GSM C30, so 'replaced' and 'original' levels of the "
            "c30_treatment dimension produce IDENTICAL fields for months "
            "before April 2012. Within W-GRACE (2002-04 -> 2017-06) the "
            "1,920-recipe factorial collapses to ~960 distinct outputs and "
            "the C30 variance share -> 0 in the variance decomposition. "
            "(2) Mission-bridging gap (2017-07 -> 2018-05) is encoded by "
            "the absence of entries in APPENDIX_MONTHS; no NaN sentinels "
            "are written, time dim has 272 entries. "
            "(3) Late-mission GRACE (2016-2017) has degraded solutions and "
            "two months ship two short-arc files; io_shm.py selects the arc "
            "whose midpoint is closest to day-15 of the month.")
    return nc


def _iter_variants(centre: str, filters: tuple,
                    window_months: list,
                    tn14, gia_rates: dict, geocenter_series: dict,
                    aoi: cfg.AOI):
    """Yield (t_idx, tr_idx, f_idx, g_idx, c20_idx, c30_idx, gc_idx, ewh_slice)
    for every (month, recipe) combination. Used by both Pass 1 and Pass 2."""
    month_to_t = {ym: i for i, ym in enumerate(window_months)}
    for tr_idx, trunc in enumerate(cfg.TRUNCATIONS):
        months = enumerate_months(centre, trunc, months=window_months)
        for m in months:
            t_idx = month_to_t[(m.year, m.month)]
            harm_raw = load_shm(m.file_path, lmax=trunc)
            for c20_idx, c20 in enumerate(cfg.C20_TREATMENTS):
                for c30_idx, c30 in enumerate(cfg.C30_TREATMENTS):
                    harm_corr = replace_low_degrees(
                        harm_raw, tn14, m, c20, c30)
                    for gc_idx, gc_name in enumerate(cfg.GEOCENTER_TREATMENTS):
                        harm_gc = apply_geocenter(
                            harm_corr, geocenter_series[gc_name], m)
                        for f_idx, filt in enumerate(cfg.FILTERS):
                            if filt not in filters:
                                continue
                            harm_filt = apply_filter(harm_gc, filt, trunc)
                            ewh = synthesise_to_aoi(harm_filt, aoi, trunc)
                            for g_idx, gia in enumerate(cfg.GIA_MODELS):
                                ewh_corr = subtract_gia_at_epoch(
                                    ewh, gia_rates[gia], m,
                                    cfg.GIA_REF_EPOCH, aoi)
                                yield (t_idx, tr_idx, f_idx, g_idx,
                                        c20_idx, c30_idx, gc_idx, ewh_corr,
                                        m)


def compute_centre_global(centre: str,
                           filters: Iterable[str] = cfg.FILTERS,
                           window: str = "poc",
                           out_dir: Path = GLOBAL_OUT_DIR,
                           max_months: int | None = None) -> Path:
    """Two-pass synthesis: Pass 1 accumulates time means; Pass 2 writes anomalies."""
    aoi = GLOBAL_AOI
    filters = tuple(filters)
    if not set(filters).issubset(cfg.FILTERS):
        raise ValueError(f"Unknown filter(s) in {filters!r}")

    window_months = cfg.months_for_window(window)
    if max_months is not None:
        window_months = window_months[:max_months]
    n_t = len(window_months)
    times = np.array(
        [np.datetime64(f"{y:04d}-{m:02d}-15") for (y, m) in window_months],
        dtype="datetime64[D]")
    n_tr = len(cfg.TRUNCATIONS)
    n_f = len(cfg.FILTERS)
    n_g = len(cfg.GIA_MODELS)
    n_c20 = len(cfg.C20_TREATMENTS)
    n_c30 = len(cfg.C30_TREATMENTS)
    n_gc = len(cfg.GEOCENTER_TREATMENTS)
    nlat, nlon = aoi.lat.size, aoi.lon.size

    tn14 = load_tn14()
    gia_rates: dict[str, "object | None"] = {"none": None}
    for name in cfg.GIA_MODELS:
        if name != "none":
            gia_rates[name] = load_gia_rate(name)
    geocenter_series = {name: load_geocenter(name)
                         for name in cfg.GEOCENTER_TREATMENTS}

    # ---------- Pass 1: synthesise + accumulate per-recipe-cell mean ----------
    # Shape excludes time: (truncation, filter, gia_model, c20, c30, geocenter,
    # lat, lon). 1.77 GB at float64+int32; comfortably in RAM.
    print(f"[{centre} global] Pass 1: accumulating time means "
          f"({n_t} months x {sum(1 for _ in filters)} active filters)",
          flush=True)
    running_sum = np.zeros((n_tr, n_f, n_g, n_c20, n_c30, n_gc, nlat, nlon),
                            dtype=np.float64)
    running_cnt = np.zeros_like(running_sum, dtype=np.int32)

    t0 = _time.time()
    months_done = set()
    for (t_idx, tr_idx, f_idx, g_idx, c20_idx, c30_idx, gc_idx,
         ewh_slice, m) in _iter_variants(
            centre, filters, window_months, tn14, gia_rates,
            geocenter_series, aoi):
        finite = np.isfinite(ewh_slice)
        running_sum[tr_idx, f_idx, g_idx, c20_idx, c30_idx, gc_idx] += \
            np.where(finite, ewh_slice, 0.0)
        running_cnt[tr_idx, f_idx, g_idx, c20_idx, c30_idx, gc_idx] += \
            finite.astype(np.int32)
        key = (m.year, m.month, tr_idx)
        if key not in months_done:
            months_done.add(key)
            elapsed = (_time.time() - t0) / 60
            print(f"[{centre} global] [P1] {cfg.TRUNC_MNEM[cfg.TRUNCATIONS[tr_idx]]} "
                   f"{m.year}-{m.month:02d}   elapsed {elapsed:.1f} min",
                   flush=True)

    print(f"[{centre} global] Pass 1 done in "
           f"{(_time.time() - t0)/60:.1f} min", flush=True)

    with np.errstate(invalid="ignore", divide="ignore"):
        running_mean = np.where(running_cnt > 0,
                                  running_sum / running_cnt,
                                  np.nan)
    del running_sum, running_cnt

    # ---------- Pass 2: re-synthesise + write anomalies -----------------------
    # KEY: write one full (t_idx, tr_idx) chunk at a time. NetCDF chunking is
    # (1, 1, all variants, all lat, all lon); a per-recipe sub-region write
    # would force decompress-recompress on every variant (~640x overhead).
    # Non-POC windows get a filename suffix so the appendix cube cannot clobber
    # the headline POC cube.
    suffix = "" if window == "poc" else f"_{window}"
    out_path = out_dir / f"grace_ewh_global_{centre.lower()}{suffix}.nc"
    nc = _create_skeleton(out_path, centre, n_t, aoi, times, filters,
                          window=window)
    ewh_var = nc.variables["ewh_anom"]
    print(f"[{centre} global] Pass 2: writing anomalies -> {out_path.name}",
          flush=True)

    chunk_shape = (n_f, n_g, n_c20, n_c30, n_gc, nlat, nlon)
    chunk_slab = np.full(chunk_shape, np.nan, dtype=np.float32)
    current_key: tuple | None = None

    def _flush_chunk(key, slab):
        if key is None:
            return
        t_idx_k, tr_idx_k = key
        ewh_var[t_idx_k, tr_idx_k, :, :, :, :, :, :, :] = slab

    t0 = _time.time()
    months_done = set()
    for (t_idx, tr_idx, f_idx, g_idx, c20_idx, c30_idx, gc_idx,
         ewh_slice, m) in _iter_variants(
            centre, filters, window_months, tn14, gia_rates,
            geocenter_series, aoi):
        chunk_key = (t_idx, tr_idx)
        if chunk_key != current_key:
            _flush_chunk(current_key, chunk_slab)
            chunk_slab = np.full(chunk_shape, np.nan, dtype=np.float32)
            current_key = chunk_key
            key = (m.year, m.month, tr_idx)
            if key not in months_done:
                months_done.add(key)
                elapsed = (_time.time() - t0) / 60
                print(f"[{centre} global] [P2] {cfg.TRUNC_MNEM[cfg.TRUNCATIONS[tr_idx]]} "
                       f"{m.year}-{m.month:02d}   elapsed {elapsed:.1f} min",
                       flush=True)

        mean_slice = running_mean[tr_idx, f_idx, g_idx, c20_idx,
                                    c30_idx, gc_idx]
        chunk_slab[f_idx, g_idx, c20_idx, c30_idx, gc_idx] = \
            (ewh_slice - mean_slice).astype(np.float32)

    # Flush final chunk
    _flush_chunk(current_key, chunk_slab)

    nc.close()
    print(f"[{centre} global] Pass 2 done in "
           f"{(_time.time() - t0)/60:.1f} min", flush=True)
    return out_path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--centre", required=True, choices=cfg.CENTRES)
    ap.add_argument("--filters", nargs="+", default=list(cfg.FILTERS),
                     choices=cfg.FILTERS)
    ap.add_argument("--window", default="appendix", choices=("poc", "appendix"),
                     help="appendix (default): full 2002-2025 GRACE+GRACE-FO "
                          "record, the paper's headline window. poc: legacy "
                          "GRACE-FO-only window, kept for the robustness "
                          "comparison.")
    ap.add_argument("--out_dir", type=Path, default=GLOBAL_OUT_DIR)
    ap.add_argument("--max_months", type=int, default=None,
                     help="Validation only: cap number of months per "
                          "truncation. Use ~2 for a quick smoke run.")
    args = ap.parse_args(argv)

    tag = f"{args.centre}/global/{args.window}"
    months_set = cfg.months_for_window(args.window)
    n_variants = cfg.variant_matrix_size_per_month()
    n_actual_months = (args.max_months
                        if args.max_months is not None else len(months_set))
    t0 = _time.time()
    print(f"[{tag}] starting global pipeline over "
          f"{n_actual_months} months x {n_variants} variants = "
          f"{n_actual_months * n_variants:,} fields "
          f"(global 0.5 deg, float32 + zlib L4, two-pass)", flush=True)
    out_path = compute_centre_global(args.centre, filters=args.filters,
                                       window=args.window,
                                       out_dir=args.out_dir,
                                       max_months=args.max_months)
    elapsed = _time.time() - t0
    size_gb = out_path.stat().st_size / 1e9
    print(f"[{tag}] finished in {elapsed/60:.1f} min, wrote "
           f"{out_path.name} ({size_gb:.2f} GB)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
