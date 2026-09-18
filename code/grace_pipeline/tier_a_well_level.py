"""Tier-A well-level test (W1 / W2 / W3) -- mirrors POC Result 4 / Result 10.

Per-well, per-recipe time-series alignment between GRACE EWH and Jasechko
2024 annual depth-to-water at the global cohort. Three tests:

  W1: per-recipe global mean Pearson r between annual GRACE EWH and
      annual well DTW across cohort wells with >=5 overlapping years
      in the chosen analysis window (poc 2018-2025 | appendix 2002-2025 |
      w_grace 2002-2017 | w_fo 2018-2025).
  W2: per-recipe Spearman rho across wells between (GRACE trend at
      well cell) and (well's Theil-Sen trend in m/yr).
  W3: cell dominance R vs per-well |GRACE - well| trend-residual IQR
      Spearman rho -- global analog of POC Result 10 Test 1 (which
      was +0.52 in the IGP).

Sign convention: well DepthToWater_m is POSITIVE for deepening (water
level drops). GRACE EWH is POSITIVE for mass gain. So the well-DTW
trend and the GRACE-EWH trend have OPPOSITE signs when the modalities
agree. We flip the well sign in the analysis (multiply by -1) so a
positive Pearson r / Spearman rho means GRACE and wells agree.

Output: datasets/output/tier_a/
    well_level_per_recipe_skill.csv         per (centre, recipe) global stats
    well_level_global_best_worst.csv        top-5 / bottom-5 ranking
    well_level_per_well_residual.csv        per (well, recipe) for W3
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
DEFAULT_WELLS_CSV = _CODE_REPO / "datasets" / "output" / "tier_a" / "wells_in_cohort.csv"
DEFAULT_DTW_CSV = _CODE_REPO / "datasets" / "jasechko_2024" / "groundwater_levels" / "AnnualDepthToGroundwater.csv"
DEFAULT_AQUIFER_METRICS_CSV = _CODE_REPO / "datasets" / "output" / "tier_a" / "jasechko_aquifer_metrics.csv"
DEFAULT_OUT_DIR = _CODE_REPO / "datasets" / "output" / "tier_a"

# Same window registry as tier_a_per_aquifer.py. Kept as a duplicate (rather
# than imported) so each script remains independently runnable from the CLI.
WINDOW_REGISTRY: dict[str, dict] = {
    "poc":      {"suffix": "",          "tmin": None, "tmax": None,
                  "out_subdir": "tier_a"},
    "appendix": {"suffix": "_appendix", "tmin": None, "tmax": None,
                  "out_subdir": "tier_a_full"},
    "w_full":   {"suffix": "_appendix", "tmin": None, "tmax": None,
                  "out_subdir": "tier_a_full"},
    "w_grace":  {"suffix": "_appendix",
                  "tmin": np.datetime64("1900-01-01"),
                  "tmax": np.datetime64("2017-06-30"),
                  "out_subdir": "tier_a_grace_only"},
    "w_fo":     {"suffix": "_appendix",
                  "tmin": np.datetime64("2018-06-01"),
                  "tmax": np.datetime64("2100-01-01"),
                  "out_subdir": "tier_a_fo"},
}


def _ols_slope(y: np.ndarray, t: np.ndarray) -> np.ndarray:
    """OLS slope along first axis. y shape (n_t, ...); returns (...) ."""
    finite = np.isfinite(y)
    n_ok = finite.sum(axis=0)
    y_f = np.where(finite, y, 0.0)
    t_b = t[(slice(None),) + (None,) * (y.ndim - 1)]
    t_f = np.where(finite, t_b, 0.0)
    s_y = y_f.sum(axis=0); s_t = t_f.sum(axis=0)
    mean_y = np.where(n_ok > 0, s_y / np.maximum(n_ok, 1), np.nan)
    mean_t = np.where(n_ok > 0, s_t / np.maximum(n_ok, 1), np.nan)
    yc = np.where(finite, y - mean_y, 0.0)
    tc = np.where(finite, t_b - mean_t, 0.0)
    num = (yc * tc).sum(axis=0)
    den = (tc * tc).sum(axis=0)
    return np.where((n_ok >= 3) & (den > 0), num / np.maximum(den, 1e-30),
                     np.nan)


def reduce_centre_well_level(cube_path: Path, wells_df: pd.DataFrame,
                              dtw_df: pd.DataFrame, min_overlap_years: int = 5,
                              tmin: "np.datetime64 | None" = None,
                              tmax: "np.datetime64 | None" = None,
                              ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per (well, recipe) Pearson r + trend; aggregated per recipe.

    `tmin`/`tmax` (inclusive bounds, np.datetime64[D]) restrict the
    well-level skill computation to a sub-range of the cube's time axis.
    Used for Path C W-GRACE / W-FO window-matched well-skill panels.
    """
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
    print(f"\n[well-level {centre}] reading cube {cube_path.name}", flush=True)
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
    n_t = len(times_jd)
    n_tr = len(nc.dimensions["truncation"])
    n_f = len(nc.dimensions["filter"])
    n_g = len(nc.dimensions["gia_model"])
    n_c20 = len(nc.dimensions["c20_treatment"])
    n_c30 = len(nc.dimensions["c30_treatment"])
    n_gc = len(nc.dimensions["geocenter"])
    n_recipes = n_tr * n_f * n_g * n_c20 * n_c30 * n_gc
    n_rec_per_tr = n_recipes // n_tr

    # Convert NetCDF time (days since 1970-01-01) to year-month
    times_dt = pd.to_datetime(times_jd, unit="D", origin="1970-01-01")
    times_year = times_dt.year.values

    nlat = len(nc.dimensions["lat"])
    nlon = len(nc.dimensions["lon"])

    # Unique cells these wells occupy
    cells = wells_df[["lat_idx", "lon_idx"]].drop_duplicates().reset_index(drop=True)
    cells["cell_idx"] = range(len(cells))
    n_cells = len(cells)
    lat_arr = cells["lat_idx"].to_numpy()
    lon_arr = cells["lon_idx"].to_numpy()
    print(f"[well-level {centre}] {n_cells:,} unique cells to read", flush=True)

    # Buffer for monthly EWH at each unique cell, all recipes
    # Shape (n_t, n_recipes, n_cells) float32, ~300 MB at 1218 cells
    buf = np.full((n_t, n_recipes, n_cells), np.nan, dtype=np.float32)

    ewh = nc.variables["ewh_anom"]
    t_start = _time.time()
    for new_t_idx, t_idx in enumerate(t_indices):
        for tr_idx in range(n_tr):
            chunk = np.asarray(
                ewh[t_idx, tr_idx, :, :, :, :, :, :, :]
            ).reshape(n_rec_per_tr, nlat, nlon)
            r_lo = tr_idx * n_rec_per_tr
            r_hi = r_lo + n_rec_per_tr
            buf[new_t_idx, r_lo:r_hi, :] = chunk[:, lat_arr, lon_arr]
        if (new_t_idx + 1) % 20 == 0 or new_t_idx == n_t - 1:
            el = (_time.time() - t_start) / 60
            print(f"  [{centre}] cube read {new_t_idx+1}/{n_t}, {el:.1f} min",
                  flush=True)
    nc.close()
    print(f"[well-level {centre}] cube buffer ready in "
          f"{(_time.time() - t_start)/60:.1f} min", flush=True)

    # Aggregate monthly -> annual per (recipe, cell)
    years_unique = np.unique(times_year)
    n_y = len(years_unique)
    annual_ewh = np.full((n_y, n_recipes, n_cells), np.nan, dtype=np.float32)
    for i, yr in enumerate(years_unique):
        sel = (times_year == yr)
        annual_ewh[i] = np.nanmean(buf[sel], axis=0)
    del buf
    print(f"[well-level {centre}] annual aggregation done "
          f"({n_y} years: {years_unique[0]}-{years_unique[-1]})", flush=True)

    # GRACE per-recipe annual trend at each cell
    t_years_grace = years_unique.astype(np.float64) - years_unique[0]
    grace_slopes_per_cell = _ols_slope(annual_ewh, t_years_grace)  # (recipes, cells)

    # For each well, pull its annual DTW series + the cell's annual EWH per recipe
    # Join wells with their cell_idx
    wjoin = wells_df.merge(cells, on=["lat_idx", "lon_idx"], how="left")
    # Filter to wells with >=min_overlap_years in our window
    dtw_in = dtw_df[(dtw_df.IntegerYear >= years_unique[0]) &
                     (dtw_df.IntegerYear <= years_unique[-1])]
    yrs_per_well = dtw_in.groupby("StnID").IntegerYear.nunique()
    keep_ids = set(yrs_per_well[yrs_per_well >= min_overlap_years].index)
    wjoin = wjoin[wjoin.StnID.isin(keep_ids)]
    print(f"[well-level {centre}] {len(wjoin):,} wells with "
           f">={min_overlap_years}-yr overlap (of {len(wells_df):,})",
           flush=True)

    # Per-well loop: collect Pearson r + trend per recipe.
    # Vectorize across recipes within a well.
    per_well_records = []
    per_well_trend_records = []
    well_meta_records = []

    dtw_groups = dict(iter(dtw_in.groupby("StnID")))
    t_loop_start = _time.time()

    for k, row in enumerate(wjoin.itertuples(index=False)):
        stn = row.StnID
        cell_i = row.cell_idx
        well_yrs = dtw_groups[stn].sort_values("IntegerYear")
        # FLIP SIGN: well DTW positive = deepening; we want positive=mass-loss
        # to align with GRACE EWH negative = mass-loss. Multiply DTW by -1 so
        # positive r means "they agree".
        well_y = well_yrs.IntegerYear.to_numpy()
        well_v = -well_yrs.DepthToWater_m.to_numpy()
        # Align to our years
        idx_in_grid = np.searchsorted(years_unique, well_y)
        valid = (idx_in_grid >= 0) & (idx_in_grid < n_y) & \
                 (years_unique[np.clip(idx_in_grid, 0, n_y - 1)] == well_y)
        idx_in_grid = idx_in_grid[valid]
        well_v = well_v[valid]
        if len(well_v) < 3:
            continue

        # Annual GRACE at this cell, all recipes: shape (overlap_yrs, n_recipes)
        grace_at_well = annual_ewh[idx_in_grid, :, cell_i]   # (n_yr, n_recipes)

        # Pearson r per recipe (vectorized)
        gx = grace_at_well - np.nanmean(grace_at_well, axis=0)
        wy = well_v - np.mean(well_v)
        num = np.nansum(gx * wy[:, None], axis=0)
        denG = np.sqrt(np.nansum(gx * gx, axis=0))
        denW = np.sqrt(np.sum(wy * wy))
        with np.errstate(invalid="ignore", divide="ignore"):
            r = num / np.maximum(denG * denW, 1e-30)
        r[~np.isfinite(r)] = np.nan

        # Well Theil-Sen trend; vectorized OLS works too (annual)
        if len(well_v) >= 3:
            t_w = well_y[valid] - well_y[valid][0]
            well_slope = (np.sum((t_w - t_w.mean()) * (well_v - well_v.mean())) /
                          max(np.sum((t_w - t_w.mean())**2), 1e-30))
        else:
            well_slope = np.nan

        # GRACE per-recipe slope already computed for the full series; use that.
        grace_slope_at_well = grace_slopes_per_cell[:, cell_i]   # (n_recipes,)

        # Store per (well, recipe) Pearson r (compactly) -- only keep recipe r
        per_well_records.append({
            "centre": centre,
            "StnID": stn,
            "Study_area": row.Study_area,
            "lat": row.Lat,
            "lon": row.Lon,
            "n_yr_overlap": int(np.sum(valid)),
            "well_trend_neg_dtw_per_yr": float(well_slope),
            # mean / max / min of r across recipes is reported in summary;
            # full r vector is huge but useful -- emit as compact pickle? keep means here
            "recipe_r_mean": float(np.nanmean(r)),
            "recipe_r_best": float(np.nanmax(r)),
            "recipe_r_worst": float(np.nanmin(r)),
        })
        # For W3: per-well |GRACE trend - well trend| IQR across recipes
        residual = grace_slope_at_well - well_slope
        per_well_trend_records.append({
            "centre": centre,
            "StnID": stn,
            "trend_residual_iqr": float(
                np.nanpercentile(residual, 75) - np.nanpercentile(residual, 25)),
            "trend_residual_mean": float(np.nanmean(residual)),
            "well_trend_neg_dtw_per_yr": float(well_slope),
        })

        # For per-recipe global stats accumulate at the end -- save the r vector
        # via a temporary list per recipe
        well_meta_records.append((stn, r, grace_slope_at_well, well_slope))

        if (k + 1) % 1000 == 0:
            el = (_time.time() - t_loop_start) / 60
            print(f"  [{centre}] processed {k+1}/{len(wjoin)} wells "
                   f"({el:.1f} min)", flush=True)

    print(f"[well-level {centre}] well loop done in "
           f"{(_time.time() - t_loop_start)/60:.1f} min", flush=True)

    # Aggregate per-recipe global stats
    n_wells = len(well_meta_records)
    if n_wells == 0:
        return pd.DataFrame(), pd.DataFrame()
    r_matrix = np.full((n_wells, n_recipes), np.nan, dtype=np.float32)
    grace_slope_per_well_recipe = np.full((n_wells, n_recipes), np.nan,
                                            dtype=np.float32)
    well_slopes = np.full(n_wells, np.nan, dtype=np.float32)
    for i, (stn, r, gs, ws) in enumerate(well_meta_records):
        r_matrix[i] = r
        grace_slope_per_well_recipe[i] = gs
        well_slopes[i] = ws

    # W1: per-recipe mean Pearson r
    recipe_mean_r = np.nanmean(r_matrix, axis=0)
    recipe_median_r = np.nanmedian(r_matrix, axis=0)
    recipe_n_valid = np.sum(np.isfinite(r_matrix), axis=0)

    # W2: per-recipe Spearman rho across wells between GRACE trend and well trend
    # Compute via rankdata, vectorized
    from scipy.stats import rankdata
    rho_per_recipe = np.full(n_recipes, np.nan, dtype=np.float32)
    well_rank = rankdata(well_slopes, nan_policy="omit")
    for ri in range(n_recipes):
        gs = grace_slope_per_well_recipe[:, ri]
        ok = np.isfinite(gs) & np.isfinite(well_slopes)
        if ok.sum() < 5:
            continue
        with np.errstate(invalid="ignore"):
            rho_per_recipe[ri] = np.corrcoef(
                rankdata(gs[ok]), rankdata(well_slopes[ok]))[0, 1]

    # Recipe key table
    trunc_vals = [str(s) for s in nc4.Dataset(cube_path, "r").variables["truncation"][:]]
    nc_h = nc4.Dataset(cube_path, "r")
    filt_vals = [str(s) for s in nc_h.variables["filter"][:]]
    gia_vals = [str(s) for s in nc_h.variables["gia_model"][:]]
    c20_vals = [str(s) for s in nc_h.variables["c20_treatment"][:]]
    c30_vals = [str(s) for s in nc_h.variables["c30_treatment"][:]]
    gc_vals = [str(s) for s in nc_h.variables["geocenter"][:]]
    nc_h.close()
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
    recipe_keys["recipe_idx_per_centre"] = np.arange(n_recipes)
    recipe_keys["centre"] = centre
    recipe_keys["recipe_mean_r"] = recipe_mean_r
    recipe_keys["recipe_median_r"] = recipe_median_r
    recipe_keys["n_wells_valid"] = recipe_n_valid
    recipe_keys["spearman_trend_rho"] = rho_per_recipe

    return recipe_keys, pd.concat([
        pd.DataFrame(per_well_records),
        pd.DataFrame(per_well_trend_records)], axis=1
        # actually keep them separate
    ) if False else pd.DataFrame(per_well_records)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cube_dir", type=Path, default=DEFAULT_CUBE_DIR)
    ap.add_argument("--wells_csv", type=Path, default=DEFAULT_WELLS_CSV)
    ap.add_argument("--dtw_csv", type=Path, default=DEFAULT_DTW_CSV)
    ap.add_argument("--centres", nargs="+", default=["csr", "jpl", "gfz"])
    ap.add_argument("--out_dir", type=Path, default=None,
                     help="Default depends on --window: tier_a (poc), "
                          "tier_a_full (appendix), tier_a_grace_only (w_grace).")
    ap.add_argument("--window", default="poc",
                     choices=list(WINDOW_REGISTRY.keys()),
                     help="poc (default, unchanged) | appendix / w_full | "
                          "w_grace (<=2017-06-30) | w_fo (>=2018-06-01)")
    ap.add_argument("--min_overlap_years", type=int, default=5)
    args = ap.parse_args(argv)

    win = WINDOW_REGISTRY[args.window]
    cube_suffix = win["suffix"]
    tmin = win["tmin"]
    tmax = win["tmax"]
    if args.out_dir is None:
        args.out_dir = DEFAULT_OUT_DIR.parent / win["out_subdir"]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[tier_a_well_level] window={args.window!r}: "
          f"cube_suffix={cube_suffix!r} t in [{tmin}, {tmax}] "
          f"-> {args.out_dir}")
    wells = pd.read_csv(args.wells_csv)
    print(f"[well-level] loading annual DTW...")
    dtw = pd.read_csv(args.dtw_csv,
                        usecols=["StnID", "IntegerYear", "DepthToWater_m"])
    cohort_ids = set(wells.StnID.values)
    dtw = dtw[dtw.StnID.isin(cohort_ids)]
    print(f"[well-level] {len(dtw):,} DTW rows for cohort wells")

    all_recipe_summaries = []
    all_per_well = []
    for centre in args.centres:
        cube_path = args.cube_dir / f"grace_ewh_global_{centre}{cube_suffix}.nc"
        if not cube_path.exists():
            print(f"SKIP {centre}: cube not found ({cube_path})")
            continue
        recipe_summary, per_well = reduce_centre_well_level(
            cube_path, wells, dtw, args.min_overlap_years,
            tmin=tmin, tmax=tmax)
        all_recipe_summaries.append(recipe_summary)
        all_per_well.append(per_well)

    if not all_recipe_summaries:
        print("Nothing produced")
        return 1

    df_recipe = pd.concat(all_recipe_summaries, ignore_index=True)
    df_well = pd.concat(all_per_well, ignore_index=True)
    df_recipe.to_csv(args.out_dir / "well_level_per_recipe_skill.csv",
                      index=False)
    df_well.to_csv(args.out_dir / "well_level_per_well_residual.csv",
                    index=False)

    # Best/worst rankings (sign-flip not needed; positive recipe_mean_r =
    # positive r already since we negated well DTW upstream)
    df_recipe_sorted = df_recipe.sort_values("recipe_mean_r", ascending=False)
    top5 = df_recipe_sorted.head(5).copy(); top5["kind"] = "best 5"
    bot5 = df_recipe_sorted.tail(5).copy(); bot5["kind"] = "worst 5"
    top_bot = pd.concat([top5, bot5], ignore_index=True)
    top_bot.to_csv(args.out_dir / "well_level_global_best_worst.csv",
                    index=False)

    print("\n=== Well-level Tier-A test summary ===")
    print(f"n wells with >={args.min_overlap_years}-yr {args.window} overlap: "
          f"{df_well.StnID.nunique() if 'StnID' in df_well.columns else 'N/A'}")
    print(f"n recipes x centres scored: {len(df_recipe)}")
    print(f"\nW1: per-recipe mean Pearson r across wells (sign-flipped):")
    print(f"  range: [{df_recipe.recipe_mean_r.min():+.3f}, "
          f"{df_recipe.recipe_mean_r.max():+.3f}]")
    print(f"  median across recipes: {df_recipe.recipe_mean_r.median():+.3f}")
    print(f"\nW2: per-recipe Spearman rho(GRACE trend at well, well trend):")
    print(f"  range: [{df_recipe.spearman_trend_rho.min():+.3f}, "
          f"{df_recipe.spearman_trend_rho.max():+.3f}]")
    print(f"  median across recipes: {df_recipe.spearman_trend_rho.median():+.3f}")

    print(f"\nBest-5 recipes (by mean r):")
    print(top5[["centre", "filter", "gia_model", "recipe_mean_r",
                "spearman_trend_rho"]].to_string(index=False))
    print(f"\nWorst-5 recipes:")
    print(bot5[["centre", "filter", "gia_model", "recipe_mean_r",
                "spearman_trend_rho"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
