"""Tier-A spatial-leakage robustness (Supplementary Materials, leakage test).

GRACE smoothing spreads signal across basin boundaries. If the sign paradox
were a leakage artefact, stripping the leakage-prone boundary cells (erosion)
should sharpen the verdict, and adding a contaminating ring (dilation) should
blur it. We test this directly: recompute every aquifer's 1,920-recipe trend
ensemble on three masks --

    eroded   : interior cells only (binary erosion, 8-connectivity) -- the
               least boundary-contaminated footprint;
    original : the cohort polygon mask (as used everywhere else);
    buffered : original + a one-cell ring (binary dilation, 8-connectivity) --
               the most boundary-contaminated footprint.

Because the OLS trend is a linear functional of the time series and the
per-aquifer trend is the spatial mean of per-cell trends (exactly as in
tier_a_per_aquifer.py: per-cell OLS slope, then nanmean over the aquifer's
cells), all three masks are evaluated from a single per-cell slope field --
no need to re-trend per mask. We stream each centre cube once, accumulate
per-(recipe, cell) OLS sufficient statistics over the union of the dilated
masks, form per-cell slopes, then aggregate per aquifer per mask version.

The `original` version reproduces jasechko_per_recipe_trends.csv to machine
precision (a built-in validation), confirming the per-cell pipeline matches
the canonical reducer.

Output (one row per aquifer x mask_version):
    datasets/output/robustness/leakage_mask_sensitivity.csv

Usage
-----
    cd <repository>/code
    python -m grace_pipeline.tier_a_leakage            # all 3 centres, 276 cohort
    python -m grace_pipeline.tier_a_leakage --centres csr   # smoke (1 centre)

Needs the three per-centre recipe cubes (grace_ewh_global_<centre>_appendix.nc
in $GRACE_DATA_ROOT/processed/grace/global/, track 3 of the README).
"""
from __future__ import annotations

import argparse
import sys
import time as _time
from pathlib import Path

import numpy as np
import pandas as pd
import netCDF4 as nc4
from scipy import ndimage

from . import config as cfg

REPO = cfg.REPO_ROOT
DEFAULT_CUBE_DIR = cfg.OUT_DIR / "global"
DEFAULT_MASKS_CSV = REPO / "datasets" / "jasechko_2024" / "cell_masks_0p5deg_276.csv"
DEFAULT_CANON_PR = cfg.existing_table(
    REPO / "datasets" / "output" / "tier_a_full" / "jasechko_per_recipe_trends.csv")
DEFAULT_OUT = REPO / "datasets" / "output" / "robustness" / "leakage_mask_sensitivity.csv"

MASK_VERSIONS = ("eroded", "original", "buffered")
# 4-connectivity (rook): a cell's interior/boundary status is set by its N/S/E/W
# neighbours. Retains an eroded interior for 159/276 aquifers (vs 89 for the
# 8-connectivity queen rule), keeping good coverage of the small Tier-2 basins
# the leakage concern is really about.
STRUCT = ndimage.generate_binary_structure(2, 1)


def build_masks(masks_df: pd.DataFrame, nlat: int, nlon: int):
    """Return per-aquifer cell lists for each mask version + the union of cells.

    For each aquifer we rasterise its cells onto an (nlat, nlon) bool grid, then
    erode / dilate by one cell (8-connectivity). Returns:
      aq_masks : list of dicts {aquifer_idx, Study_area, version -> (rows, cols)}
      union_rows, union_cols : int arrays indexing every cell any aquifer needs
      col_of   : dict (r, c) -> column index into the per-cell slope array
    """
    aq_masks = []
    needed = set()
    for aq_idx, g in masks_df.groupby("aquifer_idx", sort=True):
        grid = np.zeros((nlat, nlon), dtype=bool)
        rr = g["lat_idx"].to_numpy()
        cc = g["lon_idx"].to_numpy()
        grid[rr, cc] = True
        eroded = ndimage.binary_erosion(grid, structure=STRUCT, border_value=0)
        buffered = ndimage.binary_dilation(grid, structure=STRUCT, border_value=0)
        versions = {"eroded": eroded, "original": grid, "buffered": buffered}
        cells = {}
        for v, m in versions.items():
            r, c = np.where(m)
            cells[v] = (r, c)
            needed.update(zip(r.tolist(), c.tolist()))
        aq_masks.append({
            "aquifer_idx": int(aq_idx),
            "Study_area": g["Study_area"].iloc[0],
            "cells": cells,
        })
    union = np.array(sorted(needed))
    union_rows, union_cols = union[:, 0], union[:, 1]
    col_of = {(int(r), int(c)): i for i, (r, c) in enumerate(union)}
    return aq_masks, union_rows, union_cols, col_of


def slopes_for_centre(cube_path: Path, union_rows: np.ndarray,
                      union_cols: np.ndarray) -> tuple[np.ndarray, list]:
    """Stream one cube, return (n_recipes, n_union_cells) per-cell OLS slopes.

    NaN-aware closed-form OLS via sufficient statistics, matching
    tier_a_per_aquifer._ols_slopes (n_ok >= 5 & den > 0)."""
    nc = nc4.Dataset(cube_path, "r")
    times_jd = nc.variables["time"][:].astype("float64")
    t = (times_jd - times_jd[0]) / 365.25
    n_t = len(t)
    n_tr = len(nc.dimensions["truncation"])
    n_per_tr = (len(nc.dimensions["filter"]) * len(nc.dimensions["gia_model"])
                * len(nc.dimensions["c20_treatment"])
                * len(nc.dimensions["c30_treatment"])
                * len(nc.dimensions["geocenter"]))
    n_recipes = n_tr * n_per_tr
    nC = len(union_rows)

    # Sufficient statistics per (recipe, cell), NaN-aware.
    S_n = np.zeros((n_recipes, nC), dtype=np.float64)
    S_t = np.zeros((n_recipes, nC), dtype=np.float64)
    S_y = np.zeros((n_recipes, nC), dtype=np.float64)
    S_ty = np.zeros((n_recipes, nC), dtype=np.float64)
    S_tt = np.zeros((n_recipes, nC), dtype=np.float64)

    ewh = nc.variables["ewh_anom"]
    t0 = _time.time()
    for ti in range(n_t):
        for tr in range(n_tr):
            chunk = np.asarray(ewh[ti, tr, :, :, :, :, :, :, :]).reshape(
                n_per_tr, ewh.shape[-2], ewh.shape[-1])
            y = chunk[:, union_rows, union_cols]          # (n_per_tr, nC)
            fin = np.isfinite(y)
            yf = np.where(fin, y, 0.0)
            r0 = tr * n_per_tr
            r1 = r0 + n_per_tr
            tv = t[ti]
            S_n[r0:r1] += fin
            S_t[r0:r1] += fin * tv
            S_y[r0:r1] += yf
            S_ty[r0:r1] += yf * tv
            S_tt[r0:r1] += fin * (tv * tv)
        if (ti + 1) % 40 == 0 or ti == n_t - 1:
            print(f"    {ti + 1}/{n_t} time steps ({(_time.time()-t0)/60:.1f} min)",
                  flush=True)
    nc.close()

    den = S_n * S_tt - S_t * S_t
    num = S_n * S_ty - S_t * S_y
    with np.errstate(invalid="ignore", divide="ignore"):
        slope = np.where((S_n >= 5) & (den > 0), num / den, np.nan)
    return slope, None


def metrics_from_trends(trends: np.ndarray) -> dict:
    """Per-aquifer recipe-ensemble metrics over the 1,920 per-recipe trends.

    Matches recipe_subspace_metrics_fun (R): median, IQR, |median|, dominance R,
    sign-agreement (max of pos/neg fraction), spans-zero (min<0 & max>0)."""
    v = trends[np.isfinite(trends)]
    if v.size == 0:
        return dict(n_recipes=0, median_trend=np.nan, iqr=np.nan,
                    abs_median=np.nan, dominance_R=np.nan, sign_agree=np.nan,
                    spans_zero=np.nan, paradox_gt10=np.nan)
    med = float(np.median(v))
    q25, q75 = np.percentile(v, [25, 75])
    iqr = float(q75 - q25)
    absm = abs(med)
    fpos = float(np.mean(v > 0))
    fneg = float(np.mean(v < 0))
    minority = min(fpos, fneg)
    return dict(
        n_recipes=int(v.size),
        median_trend=med, iqr=iqr, abs_median=absm,
        dominance_R=(iqr / absm if absm > 1e-6 else np.inf),
        sign_agree=max(fpos, fneg),
        spans_zero=bool((v.min() < 0) and (v.max() > 0)),
        paradox_gt10=bool(minority > 0.10),
    )


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cube_dir", type=Path, default=DEFAULT_CUBE_DIR)
    ap.add_argument("--masks_csv", type=Path, default=DEFAULT_MASKS_CSV)
    ap.add_argument("--canon_pr", type=Path, default=DEFAULT_CANON_PR)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--centres", nargs="+", default=["csr", "jpl", "gfz"])
    args = ap.parse_args(argv)

    masks_df = pd.read_csv(args.masks_csv)
    # Grid size from the first cube.
    first_cube = args.cube_dir / f"grace_ewh_global_{args.centres[0]}_appendix.nc"
    with nc4.Dataset(first_cube, "r") as nc:
        nlat = len(nc.dimensions["lat"])
        nlon = len(nc.dimensions["lon"])
    aq_masks, urows, ucols, col_of = build_masks(masks_df, nlat, nlon)
    print(f"[leakage] {len(aq_masks)} aquifers; {len(urows)} union cells "
          f"(grid {nlat}x{nlon})", flush=True)

    # Per-(aquifer, version) recipe-trend accumulator across centres.
    # trends_acc[aq_idx][version] -> list of per-centre np arrays (n_recipes,)
    trends_acc = {a["aquifer_idx"]: {v: [] for v in MASK_VERSIONS} for a in aq_masks}
    # Keep original-version per-(centre, recipe) trends for validation.
    valid_rows = []

    for centre in args.centres:
        cube = args.cube_dir / f"grace_ewh_global_{centre}_appendix.nc"
        if not cube.exists():
            print(f"SKIP {centre}: {cube} missing", flush=True)
            continue
        print(f"[{centre}] streaming {cube.name} ...", flush=True)
        ct0 = _time.time()
        slope, _ = slopes_for_centre(cube, urows, ucols)   # (n_recipes, nC)
        print(f"[{centre}] slopes computed ({(_time.time()-ct0)/60:.1f} min); "
              "aggregating per aquifer", flush=True)
        for a in aq_masks:
            for v in MASK_VERSIONS:
                rr, cc = a["cells"][v]
                if len(rr) == 0:
                    tr_v = np.full(slope.shape[0], np.nan)
                else:
                    cols = np.array([col_of[(int(r), int(c))] for r, c in zip(rr, cc)])
                    sub = slope[:, cols]
                    with np.errstate(invalid="ignore"):
                        tr_v = np.nanmean(sub, axis=1)      # (n_recipes,) per centre
                trends_acc[a["aquifer_idx"]][v].append(tr_v)
                if v == "original":
                    for ridx, val in enumerate(tr_v):
                        valid_rows.append((centre, a["aquifer_idx"], ridx, val))

    # Build per-aquifer per-version metrics over the pooled (3-centre) ensemble.
    n_cells = {a["aquifer_idx"]: {v: len(a["cells"][v][0]) for v in MASK_VERSIONS}
               for a in aq_masks}
    study = {a["aquifer_idx"]: a["Study_area"] for a in aq_masks}
    rows = []
    for aq_idx, byv in trends_acc.items():
        for v in MASK_VERSIONS:
            pooled = np.concatenate(byv[v]) if byv[v] else np.array([])
            m = metrics_from_trends(pooled)
            m.update(aquifer_idx=aq_idx, Study_area=study[aq_idx],
                     mask_version=v, n_cells=n_cells[aq_idx][v])
            rows.append(m)
    out = pd.DataFrame(rows)
    cols_order = ["mask_version", "aquifer_idx", "Study_area", "n_cells",
                  "n_recipes", "median_trend", "iqr", "abs_median",
                  "dominance_R", "sign_agree", "spans_zero", "paradox_gt10"]
    out = out[cols_order].sort_values(["mask_version", "aquifer_idx"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"\nWrote {args.out} ({len(out)} rows)", flush=True)

    # ---- Validation: original version vs canonical per-recipe trends ----
    if args.canon_pr.exists() and valid_rows:
        vdf = pd.DataFrame(valid_rows, columns=["centre", "aquifer_idx",
                                                "recipe_idx_per_centre", "trend_new"])
        canon = pd.read_csv(args.canon_pr)
        canon["centre"] = canon["centre"].str.lower()
        m = canon.merge(vdf, on=["centre", "aquifer_idx", "recipe_idx_per_centre"],
                        how="inner")
        d = (m["trend_cm_yr"] - m["trend_new"]).abs()
        print(f"\n[validation] original vs canonical per-recipe trends: "
              f"n={len(m)}, max|diff|={np.nanmax(d):.3e}, "
              f"mean|diff|={np.nanmean(d):.3e} cm/yr", flush=True)

    # ---- Quick headline summary ----
    print("\n=== sign paradox (% spans zero) by mask version ===", flush=True)
    for v in MASK_VERSIONS:
        sub = out[out.mask_version == v]
        sz = sub["spans_zero"].dropna()
        print(f"  {v:9s}: aquifers={len(sub)}  with-interior={int(sub['n_cells'].gt(0).sum())}"
              f"  spans-zero={100*sz.mean():.1f}%  median R={sub['dominance_R'].replace(np.inf,np.nan).median():.2f}",
              flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
