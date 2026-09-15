"""Tier-A mascon comparator (global).

Mirrors the IGP-scoped mascon.py but produces a GLOBAL 0.5 deg cube of the
three official mascon products (CSR-M / JPL-M / GSFC-M) on the same grid
as the SHM cubes built by pipeline_global.py, then runs the analogues of
tier_a_per_aquifer.py and tier_a_well_level.py against the 73-aquifer
Jasechko cohort and the 5,636-well W1 cohort.

The mascon products are pre-gridded GRACE solutions in which the filter /
truncation / GIA choices are baked in by the processing centre. They are
treated here as **three additional fixed "recipes"** so the Figure 5
overlay can place them against the 1,920-recipe SHM cloud produced by
pipeline_global.py.

Outputs
-------
  $GRACE_DATA_ROOT/processed/grace/global/grace_mascon_global.nc
      Optional. (time, product, lat, lon) cube on the global 0.5 deg grid.
      Written only when --write_cube is passed (~260 MB).

  analysis/datasets/output/tier_a/mascon_aquifer_metrics.csv
      One row per aquifer; per-aquifer median trend across the 3 mascons,
      across-mascon IQR / dominance / sign-agreement.

  analysis/datasets/output/tier_a/mascon_per_aquifer_trends.csv
      One row per (aquifer, product); aquifer-mean trend per mascon.

  analysis/datasets/output/tier_a/mascon_well_level_skill.csv
      One row per mascon product; mean and median Pearson r against the
      5,636 W1 wells (annual, sign-flipped to match the SHM convention).

Usage
-----
    cd <repository>/code
    python -m grace_pipeline.tier_a_mascon [--write_cube]
"""
from __future__ import annotations

import argparse
import datetime as dt
import sys
import time as _time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc4

from . import config as cfg
from . import mascon as mascon_module


# --------------------------------------------------------------------------
# Paths and grid
# --------------------------------------------------------------------------
_CODE_REPO = cfg.REPO_ROOT

GLOBAL_AOI = cfg.AOI(lat_min=-85.0, lat_max=85.0,
                     lon_min=-180.0, lon_max=180.0, res=0.5)

DEFAULT_GLOBAL_CUBE_DIR = cfg.OUT_DIR / "global"
DEFAULT_CUBE_PATH = DEFAULT_GLOBAL_CUBE_DIR / "grace_mascon_global.nc"
DEFAULT_MASKS_CSV = (_CODE_REPO / "datasets" / "jasechko_2024"
                     / "cell_masks_0p5deg.csv")
DEFAULT_WELLS_CSV = (_CODE_REPO / "datasets" / "output" / "tier_a"
                     / "wells_in_cohort.csv")
DEFAULT_DTW_CSV = (_CODE_REPO / "datasets" / "jasechko_2024"
                   / "groundwater_levels"
                   / "AnnualDepthToGroundwater.csv")
DEFAULT_SHM_CUBE_PATH = DEFAULT_GLOBAL_CUBE_DIR / "grace_ewh_global_csr.nc"
DEFAULT_OUT_DIR = _CODE_REPO / "datasets" / "output" / "tier_a"


# --------------------------------------------------------------------------
# Re-grid + temporal harmonisation
# --------------------------------------------------------------------------
def _to_global_monthly(da: xr.DataArray, aoi: cfg.AOI = GLOBAL_AOI
                       ) -> xr.DataArray:
    """Re-grid a raw mascon DataArray (0..360 longitude convention) to the
    global 0.5 deg target grid, monthly mean, snapped to month-15.

    Mascon raw files store longitude on 0..360, so the interp target is
    expressed in 0..360 via `aoi.lon360` and the result is then re-assigned
    to the target's native -180..180 convention so the cube aligns
    byte-for-byte with the SHM global cubes."""
    da = da.sortby("lat").sortby("lon")
    da = da.interp(lat=aoi.lat, lon=aoi.lon360, method="linear")
    da = da.assign_coords(lon=("lon", aoi.lon))
    monthly = da.resample(time="MS").mean()
    new_time = pd.to_datetime(
        [dt.date(t.year, t.month, 15)
         for t in pd.to_datetime(monthly.time.values)])
    monthly = monthly.assign_coords(time=("time", new_time.values))
    return monthly


def _reindex_to_shm_time(monthly: xr.DataArray, shm_cube_path: Path
                         ) -> tuple[xr.DataArray, np.ndarray]:
    """Place the monthly series onto the exact SHM cube month axis. Returns
    (windowed_da, time_days) where time_days is in 'days since 1970-01-01'
    to match the SHM cubes byte-for-byte."""
    nc = nc4.Dataset(shm_cube_path, "r")
    time_days = nc.variables["time"][:]
    nc.close()
    target_dt = pd.to_datetime(time_days, unit="D", origin="1970-01-01")
    target = pd.to_datetime(
        [dt.date(t.year, t.month, 15) for t in target_dt])
    return monthly.reindex(time=target), np.asarray(time_days)


def build_global_cube(shm_cube_path: Path = DEFAULT_SHM_CUBE_PATH
                      ) -> tuple[xr.Dataset, np.ndarray]:
    """Build the (time, product, lat, lon) global mascon cube and return it
    along with the SHM-aligned time array (days since 1970-01-01)."""
    per_product = []
    time_days = None
    for name in mascon_module.PRODUCT_ORDER:
        print(f"[mascon-global] loading + regridding {name} ...", flush=True)
        raw = mascon_module._LOADERS[name]()
        monthly = _to_global_monthly(raw)
        windowed, time_days = _reindex_to_shm_time(monthly, shm_cube_path)
        per_product.append(windowed.expand_dims(product=[name]))
        n_ok = int(np.isfinite(windowed).any(dim=("lat", "lon")).sum())
        print(f"  [{name}] {n_ok}/{windowed.sizes['time']} months with data",
              flush=True)

    cube = xr.concat(per_product, dim="product")
    cube = cube.transpose("time", "product", "lat", "lon")
    cube.name = "ewh_anom"
    cube.attrs.update({
        "units": "cm",
        "long_name": ("GRACE mascon liquid-water-equivalent thickness "
                      "anomaly, regridded to the global 0.5 deg Tier-A grid"),
    })
    ds = xr.Dataset({"ewh_anom": cube})
    ds.attrs.update({
        "title": "GRACE mascon comparator cube -- global, POC window",
        "grid": "global 0.5 deg (340 lat x 720 lon, cell centred)",
        "products": ", ".join(
            mascon_module.MASCON_LABEL[p]
            for p in mascon_module.PRODUCT_ORDER),
        "gia_models": "; ".join(
            f"{mascon_module.MASCON_LABEL[p]}={mascon_module.MASCON_GIA[p]}"
            for p in mascon_module.PRODUCT_ORDER),
        "regrid_method": "bilinear (xarray.interp) to global 0.5 deg grid",
        "consumer": ("analysis Section 5 mascon comparator "
                     "(Figure 5)"),
        "processing_date": dt.datetime.utcnow().isoformat(timespec="seconds")
                            + "Z",
    })
    return ds, time_days


# --------------------------------------------------------------------------
# Reductions
# --------------------------------------------------------------------------
def _ols_slope_along_time(y: np.ndarray, t: np.ndarray) -> np.ndarray:
    """OLS slope along time axis of y (n_t, n_cells). Returns (n_cells,)
    in y-units per t-unit. NaN-aware."""
    finite = np.isfinite(y)
    n_ok = finite.sum(axis=0)
    y_f = np.where(finite, y, 0.0)
    t_b = t[:, None]
    t_f = np.where(finite, t_b, 0.0)
    s_y = y_f.sum(axis=0)
    s_t = t_f.sum(axis=0)
    mean_y = np.where(n_ok > 0, s_y / np.maximum(n_ok, 1), np.nan)
    mean_t = np.where(n_ok > 0, s_t / np.maximum(n_ok, 1), np.nan)
    yc = np.where(finite, y - mean_y, 0.0)
    tc = np.where(finite, t_b - mean_t, 0.0)
    num = (yc * tc).sum(axis=0)
    den = (tc * tc).sum(axis=0)
    return np.where((n_ok >= 5) & (den > 0),
                    num / np.maximum(den, 1e-30), np.nan)


def reduce_per_aquifer(cube_arr: np.ndarray, time_days: np.ndarray,
                       masks: pd.DataFrame
                       ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-aquifer aggregate of the global mascon cube.

    cube_arr has shape (n_t, n_product=3, n_lat, n_lon). For each
    aquifer we extract its cells, compute an OLS slope per (product,
    cell), then average the per-cell slopes within each product to get
    one aquifer-mean trend per mascon. The across-mascon spread
    (IQR / |median|) is the analogue of dominance R for the 3-mascon
    cohort; with only 3 products it is statistically thin but reported
    as a sanity check."""
    t_years = (time_days - time_days[0]) / 365.25
    products = list(mascon_module.PRODUCT_ORDER)
    n_t, n_prod, _, _ = cube_arr.shape

    aquifer_rows: list[dict] = []
    per_trend_rows: list[dict] = []
    for aq_idx, group in masks.groupby("aquifer_idx", sort=True):
        lat_idx = group["lat_idx"].to_numpy()
        lon_idx = group["lon_idx"].to_numpy()
        study = group["Study_area"].iloc[0]
        sub = cube_arr[:, :, lat_idx, lon_idx]
        slopes = np.empty((n_prod, sub.shape[2]), dtype=np.float32)
        for p in range(n_prod):
            slopes[p] = _ols_slope_along_time(sub[:, p, :], t_years)
        prod_trends = np.nanmean(slopes, axis=1)
        med = float(np.nanmedian(prod_trends))
        iqr = float(np.nanpercentile(prod_trends, 75)
                    - np.nanpercentile(prod_trends, 25))
        dom = (iqr / abs(med) if abs(med) > 1e-6 else np.inf)
        sign_agree = float(np.mean(np.sign(prod_trends) == np.sign(med)))

        aquifer_rows.append({
            "aquifer_idx": int(aq_idx),
            "Study_area": study,
            "n_cells": int(len(lat_idx)),
            "n_products": int(n_prod),
            "aquifer_median_trend_cmyr": med,
            "aquifer_iqr_trend_cmyr": iqr,
            "dominance_R_mascon": dom,
            "sign_agree_frac": sign_agree,
        })
        for p, prod_name in enumerate(products):
            per_trend_rows.append({
                "aquifer_idx": int(aq_idx),
                "Study_area": study,
                "product": prod_name,
                "trend_cm_yr": (float(prod_trends[p])
                                if np.isfinite(prod_trends[p]) else np.nan),
            })
    return pd.DataFrame(aquifer_rows), pd.DataFrame(per_trend_rows)


def reduce_per_well(cube_arr: np.ndarray, time_days: np.ndarray,
                    wells_df: pd.DataFrame, dtw_df: pd.DataFrame,
                    min_overlap_years: int = 5
                    ) -> pd.DataFrame:
    """Per-product W1: mean Pearson r between annual mascon EWH at each
    well's cell and the well's sign-flipped annual DTW."""
    times_dt = pd.to_datetime(time_days, unit="D", origin="1970-01-01")
    years = pd.DatetimeIndex(times_dt).year.values
    years_unique = np.unique(years)
    n_y = len(years_unique)
    n_t, n_prod, n_lat, n_lon = cube_arr.shape

    annual = np.full((n_y, n_prod, n_lat, n_lon), np.nan, dtype=np.float32)
    for i, yr in enumerate(years_unique):
        sel = years == yr
        annual[i] = np.nanmean(cube_arr[sel], axis=0)

    products = list(mascon_module.PRODUCT_ORDER)

    dtw_in = dtw_df[(dtw_df.IntegerYear >= years_unique[0])
                    & (dtw_df.IntegerYear <= years_unique[-1])]
    yrs_per_well = dtw_in.groupby("StnID").IntegerYear.nunique()
    keep_ids = set(yrs_per_well[yrs_per_well >= min_overlap_years].index)
    wells_filtered = wells_df[wells_df.StnID.isin(keep_ids)].copy()
    print(f"[mascon-well] {len(wells_filtered):,} wells with "
          f">={min_overlap_years}-yr overlap (of {len(wells_df):,})",
          flush=True)

    dtw_groups = dict(iter(dtw_in.groupby("StnID")))
    per_well_r = np.full((len(wells_filtered), n_prod), np.nan,
                         dtype=np.float32)

    for k, row in enumerate(wells_filtered.itertuples(index=False)):
        well_yrs = dtw_groups[row.StnID].sort_values("IntegerYear")
        well_y = well_yrs.IntegerYear.to_numpy()
        well_v = -well_yrs.DepthToWater_m.to_numpy()
        idx_in_grid = np.searchsorted(years_unique, well_y)
        valid = ((idx_in_grid >= 0) & (idx_in_grid < n_y) &
                 (years_unique[np.clip(idx_in_grid, 0, n_y - 1)] == well_y))
        idx_in_grid = idx_in_grid[valid]
        well_v = well_v[valid]
        if len(well_v) < 3:
            continue
        grace_at_well = annual[idx_in_grid, :,
                               int(row.lat_idx), int(row.lon_idx)]
        for p in range(n_prod):
            g = grace_at_well[:, p]
            ok = np.isfinite(g)
            if ok.sum() < 3:
                continue
            gc = g[ok] - np.mean(g[ok])
            wc = well_v[ok] - np.mean(well_v[ok])
            den = np.sqrt(np.sum(gc * gc) * np.sum(wc * wc))
            if den < 1e-30:
                continue
            per_well_r[k, p] = float(np.sum(gc * wc) / den)
        if (k + 1) % 1000 == 0:
            print(f"  [mascon-well] processed {k + 1}/{len(wells_filtered)} "
                  "wells", flush=True)

    rows = []
    for p, prod_name in enumerate(products):
        r = per_well_r[:, p]
        rows.append({
            "product": prod_name,
            "mean_r": float(np.nanmean(r)),
            "median_r": float(np.nanmedian(r)),
            "n_wells_valid": int(np.sum(np.isfinite(r))),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cube_path", type=Path, default=DEFAULT_CUBE_PATH)
    ap.add_argument("--masks_csv", type=Path, default=DEFAULT_MASKS_CSV)
    ap.add_argument("--wells_csv", type=Path, default=DEFAULT_WELLS_CSV)
    ap.add_argument("--dtw_csv", type=Path, default=DEFAULT_DTW_CSV)
    ap.add_argument("--shm_cube", type=Path, default=DEFAULT_SHM_CUBE_PATH)
    ap.add_argument("--out_dir", type=Path, default=DEFAULT_OUT_DIR)
    ap.add_argument("--write_cube", action="store_true",
                    help="Persist the global mascon cube to NetCDF "
                         "(~260 MB; optional, reductions work in-memory).")
    args = ap.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.cube_path.parent.mkdir(parents=True, exist_ok=True)

    t0 = _time.time()

    print("[mascon-global] build_global_cube ...", flush=True)
    ds, time_days = build_global_cube(args.shm_cube)
    print(f"[mascon-global] cube built in "
          f"{(_time.time() - t0) / 60:.2f} min; dims = {dict(ds.dims)}",
          flush=True)

    if args.write_cube:
        ds.to_netcdf(args.cube_path)
        sz_mb = args.cube_path.stat().st_size / 1e6
        print(f"[mascon-global] wrote {args.cube_path} ({sz_mb:.1f} MB)",
              flush=True)

    cube_arr = ds["ewh_anom"].values

    print("[mascon-global] reduce_per_aquifer ...", flush=True)
    masks = pd.read_csv(args.masks_csv)
    df_aq, df_pr = reduce_per_aquifer(cube_arr, time_days, masks)
    aq_path = args.out_dir / "mascon_aquifer_metrics.csv"
    pr_path = args.out_dir / "mascon_per_aquifer_trends.csv"
    df_aq.to_csv(aq_path, index=False)
    df_pr.to_csv(pr_path, index=False)
    print(f"  wrote {aq_path}  ({len(df_aq)} rows)", flush=True)
    print(f"  wrote {pr_path}  ({len(df_pr)} rows)", flush=True)

    print("[mascon-global] reduce_per_well ...", flush=True)
    wells_df = pd.read_csv(args.wells_csv)
    dtw_df = pd.read_csv(args.dtw_csv,
                         usecols=["StnID", "IntegerYear",
                                  "DepthToWater_m"])
    cohort_ids = set(wells_df.StnID.values)
    dtw_df = dtw_df[dtw_df.StnID.isin(cohort_ids)]
    print(f"  {len(dtw_df):,} DTW rows for cohort wells", flush=True)

    df_w = reduce_per_well(cube_arr, time_days, wells_df, dtw_df)
    w_path = args.out_dir / "mascon_well_level_skill.csv"
    df_w.to_csv(w_path, index=False)
    print(f"  wrote {w_path}", flush=True)
    print(df_w.to_string(index=False))

    print(f"[mascon-global] all done in "
          f"{(_time.time() - t0) / 60:.2f} min", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
