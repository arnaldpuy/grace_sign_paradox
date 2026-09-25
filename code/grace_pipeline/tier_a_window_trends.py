"""Per-aquifer per-recipe MONTHLY series + arbitrary-window trend slicing.

The Section 6/7 literature analysis must score every cohort paper against the
recipe ensemble sliced to ITS OWN study window. Recomputing a bespoke window
per paper by re-reading the ~110 GB appendix cube is infeasible (~30 min/centre
each). Instead we read each centre's appendix cube ONCE and persist a compact
per-aquifer, per-recipe MONTHLY-MEAN EWH series; any window's per-recipe trends
are then a months-slice + OLS in seconds.

Equivalence to the existing per-recipe trend tables
---------------------------------------------------
tier_a_per_aquifer.py computes a per-recipe aquifer trend as the MEAN over
in-polygon cells of each cell's OLS slope (`np.nanmean(slopes, axis=1)`). OLS
slope is a linear functional of the response vector for a FIXED time-design,
so when every in-polygon cell shares the same set of valid months (true here:
the synthesised global EWH grid has a value at every land cell in every
available month; gaps are whole-month mission gaps, identical across cells),

    mean_over_cells( slope(y_cell) )  ==  slope( mean_over_cells(y_cell) ) .

Hence the trend of the aquifer-MEAN monthly series equals the existing
mean-of-per-cell-slopes, to float tolerance. `--validate` proves this by
reproducing the committed W-GRACE and W-FULL tables from the monthly series.

Artefacts (datasets/output/tier_a_windows/)
  monthly_<centre>.npz   series[n_aq, n_recipes, n_t] float32 (aquifer-mean
                         EWH per recipe per month), dates (jd, days since
                         1970-01-01), aquifer_idx, Study_area, recipe_keys.
  jasechko_per_recipe_trends__<start>_<end>.csv.gz   per-window trend table,
                         SAME schema as tier_a_full/jasechko_per_recipe_trends.csv,
                         emitted for every distinct window in
                         paper_study_windows.csv (--windows mode).

Usage
-----
  cd <repository>/code            # so `-m grace_pipeline...` resolves
  # 1) one ~1.5 h pass over the 3 appendix cubes -> monthly_<centre>.npz
  python -m grace_pipeline.tier_a_window_trends --emit-monthly
  # 2) prove the slicer reproduces the committed W-GRACE / W-FULL tables
  python -m grace_pipeline.tier_a_window_trends --validate
  # 3) emit one trend CSV per distinct paper window
  python -m grace_pipeline.tier_a_window_trends --windows
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
# The committed W-FULL / W-GRACE tables were built on the 276-aquifer cohort.
DEFAULT_MASKS_CSV = (_CODE_REPO / "datasets" / "jasechko_2024"
                     / "cell_masks_0p5deg_276.csv")
OUT_DIR = _CODE_REPO / "datasets" / "output" / "tier_a_windows"
EXISTING_FULL = cfg.existing_table(_CODE_REPO / "datasets" / "output" / "tier_a_full"
                                   / "jasechko_per_recipe_trends.csv")
EXISTING_GRACE = cfg.existing_table(_CODE_REPO / "datasets" / "output" / "tier_a_grace_only"
                                    / "jasechko_per_recipe_trends.csv")
WINDOWS_CSV = (_CODE_REPO / "datasets" / "bibliometric"
               / "paper_study_windows.csv")
# The literature analysis only tests the in-cohort candidate sub-basins; the
# per-window trend CSVs are filtered to these so they stay small + committable.
BASINS_FILE = (cfg.REPO_ROOT / "datasets" / "bibliometric"
               / "literature_candidate_basins.txt")

CUBE_SUFFIX = "_appendix"          # the full GRACE+GRACE-FO record
CENTRES = ("csr", "jpl", "gfz")
EPOCH = np.datetime64("1970-01-01")
MIN_MONTHS = 5                     # match tier_a_per_aquifer._ols_slopes guard


# --------------------------------------------------------------------------
# Pass 1: emit the monthly series
# --------------------------------------------------------------------------
def emit_monthly_for_centre(cube_path: Path, masks: pd.DataFrame,
                            out_npz: Path) -> None:
    """Read one appendix cube and write the per-aquifer, per-recipe
    aquifer-MEAN monthly EWH series. Memory-light: the cross-cell mean is
    taken per chunk, so only the (n_aq, n_recipes, n_t) result is held
    (~190 MB float32), never the per-cell buffers."""
    nc = nc4.Dataset(cube_path, "r")
    times_jd = np.asarray(nc.variables["time"][:]).astype("int64")
    n_t = len(times_jd)
    n_tr = len(nc.dimensions["truncation"])
    n_f = len(nc.dimensions["filter"])
    n_g = len(nc.dimensions["gia_model"])
    n_c20 = len(nc.dimensions["c20_treatment"])
    n_c30 = len(nc.dimensions["c30_treatment"])
    n_gc = len(nc.dimensions["geocenter"])
    n_recipes_per_tr = n_f * n_g * n_c20 * n_c30 * n_gc
    n_recipes = n_tr * n_recipes_per_tr

    recipe_keys = pd.DataFrame(
        [(tr, f, g, c20, c30, gc)
         for tr in (str(s) for s in nc.variables["truncation"][:])
         for f in (str(s) for s in nc.variables["filter"][:])
         for g in (str(s) for s in nc.variables["gia_model"][:])
         for c20 in (str(s) for s in nc.variables["c20_treatment"][:])
         for c30 in (str(s) for s in nc.variables["c30_treatment"][:])
         for gc in (str(s) for s in nc.variables["geocenter"][:])],
        columns=["truncation", "filter", "gia_model", "c20_treatment",
                 "c30_treatment", "geocenter"])

    _known = {"csr", "jpl", "gfz"}
    hits = [p for p in cube_path.stem.lower().split("_") if p in _known]
    if len(hits) != 1:
        raise ValueError(f"cannot determine centre from {cube_path.name}")
    centre = hits[0].upper()

    nlat = len(nc.dimensions["lat"])
    nlon = len(nc.dimensions["lon"])
    aq_groups = []
    for aq_idx, group in masks.groupby("aquifer_idx", sort=True):
        aq_groups.append({
            "aquifer_idx": int(aq_idx),
            "Study_area": group["Study_area"].iloc[0],
            "lat_idx": group["lat_idx"].to_numpy(),
            "lon_idx": group["lon_idx"].to_numpy(),
        })
    n_aq = len(aq_groups)
    series = np.full((n_aq, n_recipes, n_t), np.nan, dtype=np.float32)

    ewh = nc.variables["ewh_anom"]
    print(f"[{centre}] {n_t} months x {n_tr} truncations, {n_aq} aquifers, "
          f"{n_recipes} recipes -> series {series.nbytes/1e6:.0f} MB",
          flush=True)
    t0 = _time.time()
    with np.errstate(invalid="ignore"):                 # all-NaN month -> nan
        for ti in range(n_t):
            for tr_idx in range(n_tr):
                chunk = np.asarray(
                    ewh[ti, tr_idx, :, :, :, :, :, :, :]
                ).reshape(n_recipes_per_tr, nlat, nlon)
                r_lo = tr_idx * n_recipes_per_tr
                for ai, a in enumerate(aq_groups):
                    cells = chunk[:, a["lat_idx"], a["lon_idx"]]    # (Rtr, ncell)
                    series[ai, r_lo:r_lo + n_recipes_per_tr, ti] = \
                        np.nanmean(cells, axis=1)
            if (ti + 1) % 20 == 0 or ti == n_t - 1:
                print(f"  [{centre}] {ti + 1}/{n_t} months "
                      f"({(_time.time()-t0)/60:.1f} min)", flush=True)
    nc.close()

    out_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_npz,
        series=series,
        dates_jd=times_jd,
        aquifer_idx=np.array([a["aquifer_idx"] for a in aq_groups]),
        study_area=np.array([a["Study_area"] for a in aq_groups], dtype=object),
        recipe_keys=recipe_keys.to_numpy().astype(object),
        recipe_cols=np.array(list(recipe_keys.columns), dtype=object),
        centre=np.array(centre))
    print(f"[{centre}] wrote {out_npz.name} "
          f"({out_npz.stat().st_size/1e6:.0f} MB, {(_time.time()-t0)/60:.1f} min)",
          flush=True)


# --------------------------------------------------------------------------
# Slicing: window -> per-recipe trend table
# --------------------------------------------------------------------------
def _month_index(dates_jd: np.ndarray) -> np.ndarray:
    """year*12 + (month-1) for each julian day, for inclusive month windows."""
    d = (EPOCH + dates_jd.astype("timedelta64[D]")).astype("datetime64[M]")
    years = d.astype("datetime64[Y]").astype(int) + 1970
    months = d.astype(int) - (years - 1970) * 12          # 0..11
    return years * 12 + months


def _ym_to_idx(ym: str) -> int:
    y, m = ym.split("-")
    return int(y) * 12 + (int(m) - 1)


def _ols_slope_lastaxis(y: np.ndarray, t: np.ndarray) -> np.ndarray:
    """NaN-aware closed-form OLS slope along the last axis.

    y: (..., n_t), t: (n_t,). Returns (...,) slopes; NaN where < MIN_MONTHS
    finite points (matches tier_a_per_aquifer._ols_slopes' n_ok>=5 guard)."""
    finite = np.isfinite(y)
    n_ok = finite.sum(axis=-1)
    yf = np.where(finite, y, 0.0)
    tb = np.broadcast_to(t, y.shape)
    tf = np.where(finite, tb, 0.0)
    n_safe = np.maximum(n_ok, 1)
    mean_y = yf.sum(axis=-1) / n_safe
    mean_t = tf.sum(axis=-1) / n_safe
    y_c = np.where(finite, y - mean_y[..., None], 0.0)
    t_c = np.where(finite, tb - mean_t[..., None], 0.0)
    num = (y_c * t_c).sum(axis=-1)
    den = (t_c * t_c).sum(axis=-1)
    return np.where((n_ok >= MIN_MONTHS) & (den > 0),
                    num / np.where(den > 0, den, 1.0), np.nan)


def slice_window(npz: dict, start: str, end: str) -> pd.DataFrame:
    """Per-recipe aquifer trend (cm/yr) for the inclusive month window
    [start, end] (YYYY-MM). Schema matches the committed trend tables."""
    dates_jd = npz["dates_jd"]
    midx = _month_index(dates_jd)
    keep = (midx >= _ym_to_idx(start)) & (midx <= _ym_to_idx(end))
    if keep.sum() < MIN_MONTHS:
        raise ValueError(f"window [{start},{end}] selects {keep.sum()} months")
    jd = dates_jd[keep]
    t_years = (jd - jd[0]) / 365.25                      # offset-free slope
    series = npz["series"][:, :, keep]                   # (n_aq, n_recipes, n_w)
    slopes = _ols_slope_lastaxis(series, t_years)        # (n_aq, n_recipes)

    centre = str(npz["centre"])
    aq_idx = npz["aquifer_idx"]
    study = npz["study_area"]
    n_aq, n_recipes = slopes.shape
    df = pd.DataFrame({
        "centre": centre,
        "aquifer_idx": np.repeat(aq_idx, n_recipes),
        "Study_area": np.repeat(study, n_recipes),
        "recipe_idx_per_centre": np.tile(np.arange(n_recipes), n_aq),
        "trend_cm_yr": slopes.reshape(-1),
    })
    rk = pd.DataFrame(npz["recipe_keys"],
                      columns=[str(c) for c in npz["recipe_cols"]])
    rk["recipe_idx_per_centre"] = np.arange(len(rk))
    return df.merge(rk, on="recipe_idx_per_centre", how="left")


def load_npz(centre: str) -> dict:
    p = OUT_DIR / f"monthly_{centre}.npz"
    if not p.exists():
        raise FileNotFoundError(f"{p} missing; run --emit-monthly first")
    return dict(np.load(p, allow_pickle=True))


def trends_for_window(start: str, end: str,
                      centres=CENTRES) -> pd.DataFrame:
    return pd.concat([slice_window(load_npz(c), start, end) for c in centres],
                     ignore_index=True)


# --------------------------------------------------------------------------
# Modes
# --------------------------------------------------------------------------
def do_emit_monthly(cube_dir: Path, masks_csv: Path, centres) -> int:
    masks = pd.read_csv(masks_csv)
    print(f"[emit-monthly] masks={masks_csv.name} "
          f"({masks['aquifer_idx'].nunique()} aquifers)")
    for c in centres:
        cube = cube_dir / f"grace_ewh_global_{c}{CUBE_SUFFIX}.nc"
        if not cube.exists():
            print(f"SKIP {c}: {cube} not found")
            continue
        emit_monthly_for_centre(cube, masks, OUT_DIR / f"monthly_{c}.npz")
    return 0


def _compare(derived: pd.DataFrame, existing_path: Path, label: str) -> bool:
    ex = pd.read_csv(existing_path)
    key = ["centre", "aquifer_idx", "recipe_idx_per_centre"]
    j = ex.merge(derived[key + ["trend_cm_yr"]], on=key,
                 suffixes=("_ex", "_new"), how="inner")
    if len(j) != len(ex):
        print(f"  [{label}] ROW MISMATCH: existing {len(ex)} vs joined {len(j)}")
    a, b = j["trend_cm_yr_ex"].to_numpy(), j["trend_cm_yr_new"].to_numpy()
    both = np.isfinite(a) & np.isfinite(b)
    nan_mismatch = int((np.isfinite(a) != np.isfinite(b)).sum())
    max_abs = float(np.nanmax(np.abs(a[both] - b[both]))) if both.any() else 0.0
    rel = np.abs(a[both] - b[both]) / np.maximum(np.abs(a[both]), 1e-6)
    max_rel = float(np.nanmax(rel)) if both.any() else 0.0
    ok = (max_abs < 1e-4) and (nan_mismatch == 0)
    print(f"  [{label}] n={len(j)} finite={both.sum()} "
          f"max|Δ|={max_abs:.2e} max relΔ={max_rel:.2e} "
          f"NaN-mismatch={nan_mismatch}  -> {'PASS' if ok else 'FAIL'}")
    return ok


def do_validate(centres) -> int:
    print("[validate] reproducing committed W-FULL and W-GRACE tables")
    full = trends_for_window("2002-04", "2025-12", centres)
    grace = trends_for_window("2002-04", "2017-06", centres)
    ok_full = _compare(full, EXISTING_FULL, "W-FULL 2002-04..2025-12")
    if EXISTING_GRACE.exists():
        ok_grace = _compare(grace, EXISTING_GRACE, "W-GRACE 2002-04..2017-06")
    else:
        print(f"  [W-GRACE] {EXISTING_GRACE} not present (not distributed); skipped")
        ok_grace = True
    if ok_full and ok_grace:
        print("[validate] PASS — slicer reproduces both tables.")
        return 0
    print("[validate] FAIL — see diffs above.")
    return 1


def do_windows(centres, basins_file: Path = None) -> int:
    if not WINDOWS_CSV.exists():
        print(f"missing {WINDOWS_CSV}; run Phase 1 first")
        return 1
    basins = None
    bf = basins_file or BASINS_FILE
    if bf and Path(bf).exists():
        basins = set(l.strip() for l in Path(bf).read_text(
            encoding="utf-8").splitlines() if l.strip())
        print(f"[windows] restricting to {len(basins)} candidate basins "
              f"from {Path(bf).name}")
    w = pd.read_csv(WINDOWS_CSV)
    distinct = (w[["win_start", "win_end"]].dropna()
                .drop_duplicates().sort_values(["win_start", "win_end"]))
    print(f"[windows] {len(distinct)} distinct windows from {len(w)} rows")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    npzs = {c: load_npz(c) for c in centres}
    for _, row in distinct.iterrows():
        s, e = row["win_start"], row["win_end"]
        df = pd.concat([slice_window(npzs[c], s, e) for c in centres],
                       ignore_index=True)
        if basins is not None:
            df = df[df["Study_area"].isin(basins)]
        out = OUT_DIR / f"jasechko_per_recipe_trends__{s}_{e}.csv.gz"
        df.to_csv(out, index=False, compression="gzip")
        print(f"  {s}..{e}: {len(df)} rows -> {out.name}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cube_dir", type=Path, default=DEFAULT_CUBE_DIR)
    ap.add_argument("--masks_csv", type=Path, default=DEFAULT_MASKS_CSV)
    ap.add_argument("--centres", nargs="+", default=list(CENTRES))
    ap.add_argument("--emit-monthly", action="store_true")
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--windows", action="store_true")
    ap.add_argument("--basins-file", type=Path, default=None,
                    help="restrict --windows trend CSVs to these Study_area "
                         "names (default: literature_candidate_basins.txt)")
    args = ap.parse_args(argv)
    centres = [c.lower() for c in args.centres]
    if args.emit_monthly:
        return do_emit_monthly(args.cube_dir, args.masks_csv, centres)
    if args.validate:
        return do_validate(centres)
    if args.windows:
        return do_windows(centres, args.basins_file)
    ap.error("choose one of --emit-monthly / --validate / --windows")


if __name__ == "__main__":
    sys.exit(main())
