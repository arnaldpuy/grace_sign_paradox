"""Tier-A per-aquifer reduction.

Reduces the global EWH anomaly cube to per-aquifer metrics for the
73-aquifer Jasechko cohort. Outputs two CSVs that downstream R-side
visualisation scripts read directly:

    datasets/output/tier_a/jasechko_aquifer_metrics.csv
        One row per (centre, aquifer). Columns include:
          - aquifer_median_trend_cmyr
          - aquifer_iqr_trend_cmyr
          - dominance_R (= iqr / |median|)
          - sign_agree_frac
          - n_cells, n_recipes
          - frac_below_-0.5, -1.0, -1.5, -2.0  (threshold-crossings)

    datasets/output/tier_a/jasechko_per_recipe_trends.csv
        One row per (centre, aquifer, recipe). Aquifer-mean trend
        per recipe. Used downstream by the well-skill ranking and the
        gradient regression.

Performance notes:
  - For each aquifer we read its bbox sub-cube (~340 x 720 -> ~50 x 50
    cells max), apply the mask, then vectorise OLS slope computation
    across all (recipe, in-polygon-cell) pairs.
  - The bbox read is one NetCDF slab op per aquifer per centre
    (73 x 3 = 219 reads total), each at most ~50 MB compressed.
  - Vectorised OLS via the closed-form
      slope = sum((t - tmean) * (y - ymean)) / sum((t - tmean)**2)
    runs in NumPy, no Python loops over (recipe, cell).

Usage
-----
    cd <repository>/code
    python -m grace_pipeline.tier_a_per_aquifer
"""
from __future__ import annotations

import argparse
import sys
import time as _time
from pathlib import Path

import numpy as np
import pandas as pd
import netCDF4 as nc4

from . import config as cfg


_CODE_REPO = cfg.REPO_ROOT

DEFAULT_CUBE_DIR = cfg.OUT_DIR / "global"
DEFAULT_MASKS_CSV = _CODE_REPO / "datasets" / "jasechko_2024" / "cell_masks_0p5deg.csv"
DEFAULT_COHORT_CSV = _CODE_REPO / "datasets" / "jasechko_2024" / "cohort_50K_n73.csv"
DEFAULT_OUT_DIR = _CODE_REPO / "datasets" / "output" / "tier_a"

DEFAULT_THRESHOLDS_CMYR = (-0.5, -1.0, -1.5, -2.0)

# CLI window key -> (cube suffix, time bounds, default out subdir). Default "appendix"
# is the headline window; `poc` = legacy GRACE-FO cube;
# `w_grace`/`w_fo` slice it for the bibliometric per-window trends.
WINDOW_REGISTRY: dict[str, dict] = {
    "poc":      {"suffix": "",          "tmin": None, "tmax": None,
                  "out_subdir": "tier_a"},
    "appendix": {"suffix": "_appendix", "tmin": None, "tmax": None,
                  "out_subdir": "tier_a_full"},
    "w_full":   {"suffix": "_appendix", "tmin": None, "tmax": None,
                  "out_subdir": "tier_a_full"},
    "w_grace":  {"suffix": "_appendix",
                  "tmin": None, "tmax": np.datetime64("2017-06-30"),
                  "out_subdir": "tier_a_grace_only"},
    "w_fo":     {"suffix": "_appendix",
                  "tmin": np.datetime64("2018-06-01"), "tmax": None,
                  "out_subdir": "tier_a_fo"},
}


def _ols_slopes(y: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Closed-form OLS slope across the first (time) axis of y.

    y shape: (n_t, n_recipes, n_cells)
    t shape: (n_t,)
    Returns (n_recipes, n_cells) slopes (cm/yr).
    """
    # Mask NaNs: an entry is "valid" iff finite. Build per-(recipe, cell)
    # means of t and y over only-the-valid entries.
    finite = np.isfinite(y)             # (n_t, R, C)
    n_ok = finite.sum(axis=0)           # (R, C)
    # NaN-aware sums
    y_filled = np.where(finite, y, 0.0)
    t_b = t[:, None, None]              # broadcast to (n_t, 1, 1)
    t_filled = np.where(finite, t_b, 0.0)
    sum_y = y_filled.sum(axis=0)        # (R, C)
    sum_t = t_filled.sum(axis=0)        # (R, C)
    mean_y = np.where(n_ok > 0, sum_y / np.maximum(n_ok, 1), np.nan)
    mean_t = np.where(n_ok > 0, sum_t / np.maximum(n_ok, 1), np.nan)
    # Centred sums
    y_c = np.where(finite, y - mean_y[None, :, :], 0.0)
    t_c = np.where(finite, t_b - mean_t[None, :, :], 0.0)
    num = (y_c * t_c).sum(axis=0)
    den = (t_c * t_c).sum(axis=0)
    slope = np.where((n_ok >= 5) & (den > 0), num / np.maximum(den, 1e-30),
                     np.nan)
    return slope


def reduce_centre(cube_path: Path, masks: pd.DataFrame,
                   thresholds: tuple = DEFAULT_THRESHOLDS_CMYR,
                   tmin: "np.datetime64 | None" = None,
                   tmax: "np.datetime64 | None" = None,
                   ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Reduce one global cube to per-aquifer metrics + per-recipe trends.

    The cube is chunked by (1 time, 1 truncation, full-recipe-axes,
    full-spatial). Reading per-aquifer bboxes forces full-chunk
    decompression for every aquifer (=14h total). Instead we read
    each chunk ONCE per (time, truncation), then dispatch the spatial
    extract to every aquifer in a single inner loop. Memory budget:
    one chunk = 1.88 GB. Per-aquifer accumulator across time: ~76 MB
    for the largest aquifer, ~3-5 GB total over all 73.

    `tmin` / `tmax` (inclusive bounds, np.datetime64[D]) restrict the
    reduction to a sub-range of the cube's time axis. Used by Phase 5
    of Path C to compute W-GRACE and W-FO window-matched trends from
    the same full-mission appendix cube.
    """
    nc = nc4.Dataset(cube_path, "r")
    times_jd = nc.variables["time"][:]
    times_dt = (np.datetime64("1970-01-01")
                + times_jd.astype("int64").view("timedelta64[D]"))
    keep_t = np.ones(len(times_jd), dtype=bool)
    if tmin is not None:
        keep_t &= (times_dt >= tmin)
    if tmax is not None:
        keep_t &= (times_dt <= tmax)
    t_indices = np.where(keep_t)[0]
    if len(t_indices) == 0:
        raise ValueError(
            f"time slice [{tmin}, {tmax}] selects 0 months out of "
            f"{len(times_jd)} in cube {cube_path.name}")
    times_jd = times_jd[t_indices]
    t_years = (times_jd - times_jd[0]) / 365.25
    n_t = len(t_years)
    n_tr = len(nc.dimensions["truncation"])
    n_f = len(nc.dimensions["filter"])
    n_g = len(nc.dimensions["gia_model"])
    n_c20 = len(nc.dimensions["c20_treatment"])
    n_c30 = len(nc.dimensions["c30_treatment"])
    n_gc = len(nc.dimensions["geocenter"])
    n_recipes = n_tr * n_f * n_g * n_c20 * n_c30 * n_gc

    trunc_vals = [str(s) for s in nc.variables["truncation"][:]]
    filt_vals = [str(s) for s in nc.variables["filter"][:]]
    gia_vals = [str(s) for s in nc.variables["gia_model"][:]]
    c20_vals = [str(s) for s in nc.variables["c20_treatment"][:]]
    c30_vals = [str(s) for s in nc.variables["c30_treatment"][:]]
    gc_vals = [str(s) for s in nc.variables["geocenter"][:]]
    recipe_keys = pd.DataFrame(
        [(tr, f, g, c20, c30, gc)
         for tr in trunc_vals
         for f in filt_vals
         for g in gia_vals
         for c20 in c20_vals
         for c30 in c30_vals
         for gc in gc_vals],
        columns=["truncation", "filter", "gia_model", "c20_treatment",
                 "c30_treatment", "geocenter"])
    recipe_keys["recipe_idx_per_centre"] = np.arange(len(recipe_keys))

    # Extract centre by matching against the known centre set; robust to
    # suffixed filenames like grace_ewh_global_csr_appendix.nc (Path C).
    _known_centres = {"csr", "jpl", "gfz"}
    _hits = [p for p in cube_path.stem.lower().split("_")
              if p in _known_centres]
    if len(_hits) != 1:
        raise ValueError(
            f"Cannot determine centre from {cube_path.name}: "
            f"got candidates {_hits!r}")
    centre = _hits[0].upper()
    t0_total = _time.time()

    # Group aquifer cell indices by aquifer_idx and convert to flat indices
    nlat = len(nc.dimensions["lat"])
    nlon = len(nc.dimensions["lon"])
    aq_info: list[dict] = []
    for aq_idx, group in masks.groupby("aquifer_idx", sort=True):
        lat_idx = group["lat_idx"].to_numpy()
        lon_idx = group["lon_idx"].to_numpy()
        # Per-recipe accumulator over in-polygon cells: one (n_t, n_recipes, n_cells)
        # tensor with the full flattened recipe index.
        aq_info.append({
            "aquifer_idx": int(aq_idx),
            "Study_area": group["Study_area"].iloc[0],
            "lat_idx": lat_idx,
            "lon_idx": lon_idx,
            "n_cells": len(lat_idx),
            "buffer": np.full((n_t, n_recipes, len(lat_idx)),
                                np.nan, dtype=np.float32),
        })

    ewh = nc.variables["ewh_anom"]
    n_recipes_per_tr = n_f * n_g * n_c20 * n_c30 * n_gc
    print(f"[{centre}] reading {n_t}x{n_tr} = {n_t * n_tr} cube chunks "
          f"({len(aq_info)} aquifers x {sum(a['n_cells'] for a in aq_info)} "
          "total in-polygon cells)", flush=True)

    t_read_start = _time.time()
    for new_t_idx, t_idx in enumerate(t_indices):
        for tr_idx in range(n_tr):
            # Read one full chunk: shape (n_f, n_g, n_c20, n_c30, n_gc, nlat, nlon)
            chunk = ewh[t_idx, tr_idx, :, :, :, :, :, :, :]
            chunk = np.asarray(chunk).reshape(n_recipes_per_tr, nlat, nlon)
            # Recipe-index offset within the full recipe enumeration:
            # recipe_keys order = (truncation outer, then f, g, c20, c30, gc)
            r_lo = tr_idx * n_recipes_per_tr
            r_hi = r_lo + n_recipes_per_tr
            for a in aq_info:
                # Spatial extract: chunk[recipes, lat_idx_arr, lon_idx_arr]
                a["buffer"][new_t_idx, r_lo:r_hi, :] = chunk[
                    :, a["lat_idx"], a["lon_idx"]]
        if (new_t_idx + 1) % 10 == 0 or new_t_idx == n_t - 1:
            el = (_time.time() - t_read_start) / 60
            print(f"  [{centre}] {new_t_idx + 1}/{n_t} time steps read "
                   f"(elapsed {el:.1f} min)", flush=True)

    nc.close()
    print(f"[{centre}] all chunks read in "
           f"{(_time.time() - t_read_start)/60:.1f} min, computing slopes",
           flush=True)

    aquifer_rows: list[dict] = []
    per_recipe_rows: list[dict] = []
    for a in aq_info:
        slopes = _ols_slopes(a["buffer"], t_years)  # (n_recipes, n_cells)
        # Aquifer-level metrics
        sf = slopes[np.isfinite(slopes)]
        if sf.size == 0:
            print(f"  WARNING: aq {a['aquifer_idx']} ({a['Study_area']}): "
                  f"all-NaN, skipping")
            continue
        median_cmyr = float(np.median(sf))
        q25, q75 = np.percentile(sf, [25, 75])
        iqr_cmyr = float(q75 - q25)
        sign_pos = (slopes > 0).sum(axis=0) / np.maximum(
            np.isfinite(slopes).sum(axis=0), 1)
        sign_neg = (slopes < 0).sum(axis=0) / np.maximum(
            np.isfinite(slopes).sum(axis=0), 1)
        sign_agree = np.maximum(sign_pos, sign_neg).mean()
        per_recipe_mean_slope = np.nanmean(slopes, axis=1)

        row = {
            "centre": centre,
            "aquifer_idx": a["aquifer_idx"],
            "Study_area": a["Study_area"],
            "n_cells": a["n_cells"],
            "n_recipes": int(n_recipes),
            "aquifer_median_trend_cmyr": median_cmyr,
            "aquifer_iqr_trend_cmyr": iqr_cmyr,
            "dominance_R": (iqr_cmyr / abs(median_cmyr)
                              if abs(median_cmyr) > 1e-6 else np.inf),
            "sign_agree_frac": float(sign_agree),
            "recipe_mean_trend_cmyr": float(np.nanmean(per_recipe_mean_slope)),
        }
        for thr in thresholds:
            row[f"frac_below_{thr}"] = float(
                np.mean(per_recipe_mean_slope < thr))
        aquifer_rows.append(row)
        for ridx, slope in enumerate(per_recipe_mean_slope):
            per_recipe_rows.append({
                "centre": centre,
                "aquifer_idx": a["aquifer_idx"],
                "Study_area": a["Study_area"],
                "recipe_idx_per_centre": int(ridx),
                "trend_cm_yr": float(slope) if np.isfinite(slope) else np.nan,
            })
        # Free the buffer
        a["buffer"] = None

    df_aq = pd.DataFrame(aquifer_rows)
    df_pr = pd.DataFrame(per_recipe_rows)
    df_pr = df_pr.merge(recipe_keys, on="recipe_idx_per_centre", how="left")
    print(f"[{centre}] full reduction done in "
           f"{(_time.time() - t0_total)/60:.1f} min", flush=True)
    return df_aq, df_pr


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cube_dir", type=Path, default=DEFAULT_CUBE_DIR)
    ap.add_argument("--masks_csv", type=Path, default=DEFAULT_MASKS_CSV)
    ap.add_argument("--cohort_csv", type=Path, default=DEFAULT_COHORT_CSV)
    ap.add_argument("--out_dir", type=Path, default=None,
                     help="Output directory. Defaults to datasets/output/"
                          "<window-default> -- e.g. tier_a (poc), "
                          "tier_a_full (appendix), tier_a_grace_only (w_grace).")
    ap.add_argument("--window", default="appendix",
                     choices=list(WINDOW_REGISTRY.keys()),
                     help="appendix / w_full (default): full GRACE+GRACE-FO "
                          "cube, the paper's headline window. poc: legacy "
                          "GRACE-FO-only cube. w_grace: appendix cube sliced "
                          "to <=2017-06-30 (pre-GRACE-FO papers). w_fo: "
                          "appendix cube sliced to >=2018-06-01 (GRACE-FO "
                          "era only).")
    ap.add_argument("--centres", nargs="+", default=["csr", "jpl", "gfz"])
    args = ap.parse_args(argv)

    win = WINDOW_REGISTRY[args.window]
    cube_suffix = win["suffix"]
    tmin = win["tmin"]
    tmax = win["tmax"]
    if args.out_dir is None:
        args.out_dir = DEFAULT_OUT_DIR.parent / win["out_subdir"]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[tier_a_per_aquifer] window={args.window!r}: "
          f"cube_suffix={cube_suffix!r} t in [{tmin}, {tmax}] "
          f"-> {args.out_dir}")
    masks = pd.read_csv(args.masks_csv)

    all_aq = []
    all_pr = []
    for centre in args.centres:
        cube_path = args.cube_dir / f"grace_ewh_global_{centre}{cube_suffix}.nc"
        if not cube_path.exists():
            print(f"SKIP {centre}: {cube_path} not found")
            continue
        df_aq, df_pr = reduce_centre(cube_path, masks, tmin=tmin, tmax=tmax)
        all_aq.append(df_aq)
        all_pr.append(df_pr)

    if not all_aq:
        print("No cubes found; nothing to write.")
        return 1

    df_aq_full = pd.concat(all_aq, ignore_index=True)
    df_pr_full = pd.concat(all_pr, ignore_index=True)
    aq_path = args.out_dir / "jasechko_aquifer_metrics.csv"
    pr_path = args.out_dir / "jasechko_per_recipe_trends.csv"
    df_aq_full.to_csv(aq_path, index=False)
    df_pr_full.to_csv(pr_path, index=False)
    print(f"Wrote: {aq_path} ({aq_path.stat().st_size/1e6:.2f} MB, "
          f"{len(df_aq_full)} rows)")
    print(f"Wrote: {pr_path} ({pr_path.stat().st_size/1e6:.2f} MB, "
          f"{len(df_pr_full)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
