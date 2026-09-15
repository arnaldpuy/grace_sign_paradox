"""Tier-B sensitivity: 1-degree re-synthesis + cell-level gain-factor branch.

One pass per centre over the full-mission appendix cubes that, while
reducing each aquifer's cells, emits THREE per-recipe aquifer-mean trend
sets from the same per-cell OLS slopes:

  baseline  unweighted mean over the 0.5 deg centre-in-polygon cells --
            must reproduce datasets/output/tier_a_full_73/
            jasechko_per_recipe_trends.csv (PARITY GATE, reported per
            centre; the pass aborts nothing but flags loudly).
  deg1      re-synthesis on a 1 deg grid: each 1 deg cell is the mean of
            its four 0.5 deg constituent cells, aquifer membership is
            centre-in-polygon on the 1 deg grid
            (datasets/jasechko_2024/cell_masks_1p0deg.csv). Tests whether
            a coarser, mascon-like synthesis changes the verdicts.
  gain      cell-level multiplicative scale factors (published
            TELLUS/CLM gain grid) applied to the 0.5 deg per-cell slopes
            before aquifer averaging: trend = mean_c(g_c * slope_c).
            Valid because the gain is time-independent, so
            trend(g * y) = g * trend(y) cell by cell. Cells without a
            gain value (ocean-masked coastal cells in the gain grid)
            are dropped from the mean and counted.

All three are linear functionals of the per-cell slopes, so the cube is
read ONCE per centre; per-cell slopes are also persisted
(sensitivity_cell_slopes_<centre>.npz, ~10 MB) so the gain branch can be
recomputed against any gain grid in seconds (`--gain-only`).

Outputs (datasets/output/tier_a/):
    sensitivity_1deg_per_recipe_trends.csv.gz
    sensitivity_gain_per_recipe_trends.csv.gz
    sensitivity_cell_slopes_<centre>.npz
    sensitivity_1deg_gain_report.txt

Usage
-----
    cd analysis_repo/code
    python -m grace_pipeline.tier_a_1deg_gain                 # full pass
    python -m grace_pipeline.tier_a_1deg_gain --smoke 6       # plumbing
    python -m grace_pipeline.tier_a_1deg_gain --gain-only \
        --gain_nc /path/to/CLM4.SCALE_FACTOR....nc            # post hoc
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

REPO = cfg.REPO_ROOT
CUBE_DIR = cfg.OUT_DIR / "global"
MASKS_05 = REPO / "datasets/jasechko_2024/cell_masks_0p5deg.csv"
MASKS_10 = REPO / "datasets/jasechko_2024/cell_masks_1p0deg.csv"
PARITY_CSV = REPO / "datasets/output/tier_a_full_73/jasechko_per_recipe_trends.csv"
OUT_DIR = REPO / "datasets/output/tier_a"
CENTRES = ("csr", "jpl", "gfz")
PARITY_TOL = 1e-6  # cm/yr


def _ols_slopes(y: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Closed-form OLS slope across the first (time) axis of y.

    Identical to tier_a_per_aquifer._ols_slopes (NaN-aware, >=5 valid
    months). y: (n_t, n_recipes, n_cells); t: (n_t,). Returns
    (n_recipes, n_cells) slopes in cm/yr.
    """
    finite = np.isfinite(y)
    n_ok = finite.sum(axis=0)
    y_filled = np.where(finite, y, 0.0)
    t_b = t[:, None, None]
    t_filled = np.where(finite, t_b, 0.0)
    sum_y = y_filled.sum(axis=0)
    sum_t = t_filled.sum(axis=0)
    mean_y = np.where(n_ok > 0, sum_y / np.maximum(n_ok, 1), np.nan)
    mean_t = np.where(n_ok > 0, sum_t / np.maximum(n_ok, 1), np.nan)
    y_c = np.where(finite, y - mean_y[None, :, :], 0.0)
    t_c = np.where(finite, t_b - mean_t[None, :, :], 0.0)
    num = (y_c * t_c).sum(axis=0)
    den = (t_c * t_c).sum(axis=0)
    return np.where((n_ok >= 5) & (den > 0), num / np.maximum(den, 1e-30),
                    np.nan)


def recipe_key_frame(nc: nc4.Dataset) -> pd.DataFrame:
    """Recipe enumeration identical to tier_a_per_aquifer (truncation outer)."""
    trunc_vals = [str(s) for s in nc.variables["truncation"][:]]
    filt_vals = [str(s) for s in nc.variables["filter"][:]]
    gia_vals = [str(s) for s in nc.variables["gia_model"][:]]
    c20_vals = [str(s) for s in nc.variables["c20_treatment"][:]]
    c30_vals = [str(s) for s in nc.variables["c30_treatment"][:]]
    gc_vals = [str(s) for s in nc.variables["geocenter"][:]]
    keys = pd.DataFrame(
        [(tr, f, g, c20, c30, gc)
         for tr in trunc_vals
         for f in filt_vals
         for g in gia_vals
         for c20 in c20_vals
         for c30 in c30_vals
         for gc in gc_vals],
        columns=["truncation", "filter", "gia_model", "c20_treatment",
                 "c30_treatment", "geocenter"])
    keys["recipe_idx_per_centre"] = np.arange(len(keys))
    return keys


def build_aquifer_cells(masks05: pd.DataFrame,
                        masks10: pd.DataFrame) -> list[dict]:
    """Per-aquifer union cell set + index maps for the three reductions.

    The 0.5 deg grid is lat = -84.75..84.75 (340), lon = -179.75..179.75
    (720); the 1 deg grid is lat = -84.5..84.5 (170), lon = -179.5..179.5
    (360). The 1 deg cell (I, J) covers exactly the four 0.5 deg cells
    (2I, 2J), (2I, 2J+1), (2I+1, 2J), (2I+1, 2J+1).
    """
    aq10 = {int(i): g for i, g in masks10.groupby("aquifer_idx")}
    out = []
    for aq_idx, g05 in masks05.groupby("aquifer_idx", sort=True):
        aq_idx = int(aq_idx)
        cells05 = list(zip(g05["lat_idx"].astype(int),
                           g05["lon_idx"].astype(int)))
        g10 = aq10.get(aq_idx)
        cells10 = ([] if g10 is None else
                   list(zip(g10["lat_idx"].astype(int),
                            g10["lon_idx"].astype(int))))
        const10 = [(2 * i + di, 2 * j + dj)
                   for (i, j) in cells10 for di in (0, 1) for dj in (0, 1)]
        union = sorted(set(cells05) | set(const10))
        pos = {c: k for k, c in enumerate(union)}
        out.append({
            "aquifer_idx": aq_idx,
            "Study_area": g05["Study_area"].iloc[0],
            "lat_idx": np.array([c[0] for c in union]),
            "lon_idx": np.array([c[1] for c in union]),
            "n_union": len(union),
            # baseline: positions of the 0.5 deg mask cells in the union
            "pos05": np.array([pos[c] for c in cells05]),
            # deg1: (n_cells10, 4) constituent positions
            "pos10": np.array([[pos[(2 * i + di, 2 * j + dj)]
                                for di in (0, 1) for dj in (0, 1)]
                               for (i, j) in cells10],
                              dtype=int).reshape(len(cells10), 4),
        })
    return out


def reduce_centre(cube_path: Path, aq_cells: list[dict],
                  smoke_n: int | None = None):
    """One pass: read every (time, truncation) chunk once, buffer the
    union cells per aquifer, compute per-cell slopes."""
    nc = nc4.Dataset(cube_path, "r")
    times_jd = np.asarray(nc.variables["time"][:], dtype=float)
    if smoke_n is not None:
        times_jd = times_jd[:smoke_n]
    t_years = (times_jd - times_jd[0]) / 365.25
    n_t = len(t_years)
    n_tr = len(nc.dimensions["truncation"])
    n_f = len(nc.dimensions["filter"])
    n_g = len(nc.dimensions["gia_model"])
    n_c20 = len(nc.dimensions["c20_treatment"])
    n_c30 = len(nc.dimensions["c30_treatment"])
    n_gc = len(nc.dimensions["geocenter"])
    n_recipes = n_tr * n_f * n_g * n_c20 * n_c30 * n_gc
    n_per_tr = n_recipes // n_tr
    nlat = len(nc.dimensions["lat"])
    nlon = len(nc.dimensions["lon"])
    recipe_keys = recipe_key_frame(nc)

    _known = {"csr", "jpl", "gfz"}
    centre = [p for p in cube_path.stem.lower().split("_") if p in _known]
    assert len(centre) == 1, cube_path.name
    centre = centre[0].upper()

    for a in aq_cells:
        a["buffer"] = np.full((n_t, n_recipes, a["n_union"]), np.nan,
                              dtype=np.float32)

    ewh = nc.variables["ewh_anom"]
    n_cells_total = sum(a["n_union"] for a in aq_cells)
    print(f"[{centre}] {n_t}x{n_tr} chunks, {len(aq_cells)} aquifers, "
          f"{n_cells_total} union cells "
          f"(buffers {n_t * n_recipes * n_cells_total * 4 / 1e9:.1f} GB)",
          flush=True)
    t0 = _time.time()
    for t_idx in range(n_t):
        for tr_idx in range(n_tr):
            chunk = ewh[t_idx, tr_idx, :, :, :, :, :, :, :]
            chunk = np.asarray(chunk).reshape(n_per_tr, nlat, nlon)
            r_lo = tr_idx * n_per_tr
            r_hi = r_lo + n_per_tr
            for a in aq_cells:
                a["buffer"][t_idx, r_lo:r_hi, :] = chunk[
                    :, a["lat_idx"], a["lon_idx"]]
        if (t_idx + 1) % 10 == 0 or t_idx == n_t - 1:
            el = (_time.time() - t0) / 60
            eta = el / (t_idx + 1) * (n_t - t_idx - 1)
            print(f"  [{centre}] {t_idx + 1}/{n_t} read "
                  f"({el:.1f} min, ~{eta:.0f} min left)", flush=True)
    nc.close()

    print(f"[{centre}] computing per-cell slopes", flush=True)
    for a in aq_cells:
        a["slopes"] = _ols_slopes(a["buffer"], t_years)  # (R, n_union)
        a["buffer"] = None
    return centre, recipe_keys


def reductions_from_slopes(centre: str, aq_cells: list[dict],
                           recipe_keys: pd.DataFrame,
                           gain: "np.ndarray | None" = None):
    """baseline / deg1 / gain per-recipe aquifer trends from cell slopes."""
    base_rows, deg1_rows, gain_rows = [], [], []
    for a in aq_cells:
        s = a["slopes"]                                     # (R, n_union)
        base = np.nanmean(s[:, a["pos05"]], axis=1)         # (R,)
        if a["pos10"].size:
            cell10 = np.nanmean(s[:, a["pos10"].ravel()].reshape(
                s.shape[0], -1, 4), axis=2)                 # (R, n10)
            deg1 = np.nanmean(cell10, axis=1)
        else:
            deg1 = np.full(s.shape[0], np.nan)
        n10 = a["pos10"].shape[0]
        if gain is not None:
            g = gain[a["lat_idx"][a["pos05"]], a["lon_idx"][a["pos05"]]]
            ok = np.isfinite(g)
            gtr = (np.nanmean(s[:, a["pos05"][ok]] * g[ok][None, :], axis=1)
                   if ok.any() else np.full(s.shape[0], np.nan))
            n_g_ok = int(ok.sum())
        for ridx in range(s.shape[0]):
            common = {"centre": centre, "aquifer_idx": a["aquifer_idx"],
                      "Study_area": a["Study_area"],
                      "recipe_idx_per_centre": ridx}
            base_rows.append({**common, "trend_cm_yr": float(base[ridx]),
                              "n_cells": len(a["pos05"])})
            deg1_rows.append({**common, "trend_cm_yr": float(deg1[ridx]),
                              "n_cells_1deg": n10})
            if gain is not None:
                gain_rows.append({**common, "trend_cm_yr": float(gtr[ridx]),
                                  "n_cells_gain": n_g_ok,
                                  "n_cells": len(a["pos05"])})
    mk = lambda rows: (pd.DataFrame(rows)
                       .merge(recipe_keys, on="recipe_idx_per_centre",
                              how="left"))
    return (mk(base_rows), mk(deg1_rows),
            mk(gain_rows) if gain is not None else None)


def parity_check(df_base: pd.DataFrame, parity: pd.DataFrame,
                 centre: str) -> float:
    ref = parity[parity["centre"] == centre]
    m = df_base.merge(
        ref[["Study_area", "recipe_idx_per_centre", "trend_cm_yr"]],
        on=["Study_area", "recipe_idx_per_centre"],
        suffixes=("", "_ref"))
    assert len(m) == len(df_base), (len(m), len(df_base))
    return float((m["trend_cm_yr"] - m["trend_cm_yr_ref"]).abs().max())


def load_gain_grid(gain_nc: Path) -> np.ndarray:
    """Published gain grid -> (340, 720) array on the cube's 0.5 deg grid.

    Accepts any lat/lon-gridded NetCDF with a scale-factor variable
    (name containing 'scale' or 'factor', case-insensitive); regrids by
    nearest-neighbour (each 0.5 deg cell centre takes the value of the
    containing source-grid cell). Missing values -> NaN.
    """
    nc = nc4.Dataset(gain_nc, "r")
    var_name = None
    for v in nc.variables:
        if v.lower() in ("lat", "latitude", "lon", "longitude", "time"):
            continue
        if "scale" in v.lower() or "factor" in v.lower():
            var_name = v
            break
    if var_name is None:  # fall back: first 2-D non-coordinate variable
        for v in nc.variables:
            if (nc.variables[v].ndim == 2 and
                    v.lower() not in ("lat", "latitude", "lon", "longitude")):
                var_name = v
                break
    if var_name is None:
        raise SystemExit(f"no scale-factor variable found in {gain_nc.name}: "
                         f"{sorted(nc.variables)}")
    lat_name = "lat" if "lat" in nc.variables else "latitude"
    lon_name = "lon" if "lon" in nc.variables else "longitude"
    g_lat = np.asarray(nc.variables[lat_name][:], dtype=float)
    g_lon = np.asarray(nc.variables[lon_name][:], dtype=float)
    g = nc.variables[var_name][:]
    g = np.ma.filled(g, np.nan).astype(float)
    if g.ndim == 3:
        g = g[0]
    nc.close()
    print(f"[gain] {gain_nc.name}: variable {var_name!r}, "
          f"grid {g.shape}, lat {g_lat[0]}..{g_lat[-1]}, "
          f"lon {g_lon[0]}..{g_lon[-1]}, "
          f"finite {np.isfinite(g).mean():.1%}", flush=True)

    cube_lat = np.arange(-84.75, 85.0, 0.5)
    cube_lon = np.arange(-179.75, 180.0, 0.5)
    # nearest source index per cube coordinate (handle 0..360 longitudes)
    src_lon = g_lon.copy()
    tgt_lon = cube_lon.copy()
    if src_lon.max() > 180.0:
        tgt_lon = np.where(tgt_lon < 0, tgt_lon + 360.0, tgt_lon)
    li = np.abs(g_lat[None, :] - cube_lat[:, None]).argmin(axis=1)
    lj = np.abs(src_lon[None, :] - tgt_lon[:, None]).argmin(axis=1)
    out = g[np.ix_(li, lj)]
    print(f"[gain] regridded to (340, 720); finite "
          f"{np.isfinite(out).mean():.1%}; "
          f"range [{np.nanmin(out):.2f}, {np.nanmax(out):.2f}]", flush=True)
    return out


def save_slopes_npz(path: Path, centre: str, aq_cells: list[dict]) -> None:
    payload = {"centre": np.array(centre)}
    for a in aq_cells:
        k = f"aq{a['aquifer_idx']:03d}"
        # float64: the parity gate must hold on reload (--gain-only), not
        # only on the live pass; float32 storage degrades it to ~1.5e-6.
        payload[f"{k}_slopes"] = a["slopes"].astype(np.float64)
        payload[f"{k}_lat_idx"] = a["lat_idx"]
        payload[f"{k}_lon_idx"] = a["lon_idx"]
        payload[f"{k}_pos05"] = a["pos05"]
        payload[f"{k}_pos10"] = a["pos10"]
        payload[f"{k}_name"] = np.array(a["Study_area"])
    np.savez_compressed(path, **payload)
    print(f"wrote {path} ({path.stat().st_size / 1e6:.1f} MB)", flush=True)


def load_slopes_npz(path: Path, aq_cells: list[dict]) -> str:
    z = np.load(path, allow_pickle=False)
    for a in aq_cells:
        k = f"aq{a['aquifer_idx']:03d}"
        a["slopes"] = z[f"{k}_slopes"]
        assert (z[f"{k}_lat_idx"] == a["lat_idx"]).all()
        assert (z[f"{k}_lon_idx"] == a["lon_idx"]).all()
    return str(z["centre"])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--centres", nargs="+", default=list(CENTRES))
    ap.add_argument("--smoke", type=int, default=None,
                    help="read only the first N time steps (plumbing test; "
                         "parity will NOT hold and outputs are not written)")
    ap.add_argument("--gain_nc", type=Path, default=None,
                    help="published gain grid NetCDF; if absent the gain "
                         "branch is skipped (recompute later with "
                         "--gain-only)")
    ap.add_argument("--gain-only", action="store_true",
                    help="skip the cube pass; rebuild all three trend sets "
                         "from the persisted sensitivity_cell_slopes_*.npz")
    args = ap.parse_args(argv)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    masks05 = pd.read_csv(MASKS_05)
    masks10 = pd.read_csv(MASKS_10)
    parity = pd.read_csv(PARITY_CSV)
    gain = load_gain_grid(args.gain_nc) if args.gain_nc else None
    if args.gain_only and gain is None:
        raise SystemExit("--gain-only requires --gain_nc")

    all_base, all_deg1, all_gain, report = [], [], [], []
    for c in args.centres:
        aq_cells = build_aquifer_cells(masks05, masks10)
        if args.gain_only:
            npz = OUT_DIR / f"sensitivity_cell_slopes_{c}.npz"
            centre = load_slopes_npz(npz, aq_cells)
            recipe_keys = parity[parity["centre"] == centre][
                ["recipe_idx_per_centre", "truncation", "filter",
                 "gia_model", "c20_treatment", "c30_treatment",
                 "geocenter"]].drop_duplicates().sort_values(
                "recipe_idx_per_centre").reset_index(drop=True)
        else:
            cube = CUBE_DIR / f"grace_ewh_global_{c}_appendix.nc"
            if not cube.exists():
                print(f"SKIP {c}: {cube} not found")
                continue
            centre, recipe_keys = reduce_centre(cube, aq_cells,
                                                smoke_n=args.smoke)
            if args.smoke is None:
                save_slopes_npz(
                    OUT_DIR / f"sensitivity_cell_slopes_{c}.npz",
                    centre, aq_cells)
        df_b, df_1, df_g = reductions_from_slopes(centre, aq_cells,
                                                  recipe_keys, gain=gain)
        if args.smoke is not None:
            print(f"[{centre}] SMOKE OK: base/deg1 shapes "
                  f"{df_b.shape}/{df_1.shape}, "
                  f"base NaN {df_b.trend_cm_yr.isna().mean():.1%}")
            continue
        diff = parity_check(df_b, parity, centre)
        flag = "PASS" if diff < PARITY_TOL else "FAIL <<<<<<<<<<<<"
        line = (f"[{centre}] parity vs tier_a_full_73: max|diff| = "
                f"{diff:.3e} cm/yr [{flag}]")
        print(line, flush=True)
        report.append(line)
        all_base.append(df_b)
        all_deg1.append(df_1)
        if df_g is not None:
            all_gain.append(df_g)

    if args.smoke is not None or not all_deg1:
        return 0
    deg1 = pd.concat(all_deg1, ignore_index=True)
    p1 = OUT_DIR / "sensitivity_1deg_per_recipe_trends.csv.gz"
    deg1.to_csv(p1, index=False, float_format="%.8g")
    print(f"wrote {p1} ({len(deg1)} rows)")
    if all_gain:
        gn = pd.concat(all_gain, ignore_index=True)
        p2 = OUT_DIR / "sensitivity_gain_per_recipe_trends.csv.gz"
        gn.to_csv(p2, index=False, float_format="%.8g")
        print(f"wrote {p2} ({len(gn)} rows)")
    rep = OUT_DIR / "sensitivity_1deg_gain_report.txt"
    mode = "a" if args.gain_only else "w"
    with open(rep, mode) as fh:
        fh.write("\n".join(report) + "\n")
        if args.gain_nc:
            fh.write(f"gain grid: {args.gain_nc}\n")
    print(f"wrote {rep}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
