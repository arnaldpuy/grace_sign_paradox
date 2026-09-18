"""Per-aquifer per-recipe monthly SERIES extraction (D018).

The synthetic cubes are deleted after reduction (disk budget D006), but the
cube's information content at aquifer scale is fully captured by the
mask-mean EWH series per (aquifer, recipe, month): 73 x 640 x ~250 float32
= ~47 MB per centre. Saving the series makes every later window analysis
(rolling 10-yr null, GRACE-only / GRACE-FO decomposition, seasonal fits,
noise tests) a cheap R job instead of a 9-hour cube re-synthesis.

Parity gate: OLS trends computed from the extracted series must equal the
canonical `tier_a_per_aquifer` per-recipe trends (mean-of-cell-slopes ==
slope-of-cell-mean under identical time sampling) to < 1e-4 cm/yr; the
recovery driver refuses to delete a cube otherwise.

Usage:
    04b_extract_series.py --cube <path.nc> --out <dir>
"""
from __future__ import annotations

import argparse
import pathlib
import sys
import time as _time

import numpy as np
import pandas as pd
import netCDF4 as nc4
import xarray as xr

CONTROL_ROOT = pathlib.Path(__file__).resolve().parents[2]
REPO_ROOT = CONTROL_ROOT.parent
MASKS_CSV = pathlib.Path(REPO_ROOT / "datasets/jasechko_2024/cell_masks_0p5deg.csv")


def extract(cube_path: pathlib.Path, out_dir: pathlib.Path) -> pathlib.Path:
    masks = pd.read_csv(MASKS_CSV)
    nc = nc4.Dataset(cube_path, "r")
    n_t = len(nc.dimensions["time"])
    dims = ["truncation", "filter", "gia_model", "c20_treatment",
            "c30_treatment", "geocenter"]
    sizes = [len(nc.dimensions[d]) for d in dims]
    n_tr, rest = sizes[0], int(np.prod(sizes[1:]))
    n_recipes = n_tr * rest
    labels = {d: [str(s) for s in nc.variables[d][:]] for d in dims}
    centre = [p for p in cube_path.stem.lower().split("_")
              if p in {"csr", "jpl", "gfz"}][0]

    aq = []
    for aq_idx, grp in masks.groupby("aquifer_idx", sort=True):
        aq.append({"idx": int(aq_idx),
                   "name": grp["Study_area"].iloc[0],
                   "la": grp["lat_idx"].to_numpy(),
                   "lo": grp["lon_idx"].to_numpy()})
    series = np.full((len(aq), n_recipes, n_t), np.nan, dtype=np.float32)

    ewh = nc.variables["ewh_anom"]
    nlat, nlon = len(nc.dimensions["lat"]), len(nc.dimensions["lon"])
    t0 = _time.time()
    for t in range(n_t):
        for tr in range(n_tr):
            chunk = np.ma.filled(ewh[t, tr], np.nan).reshape(rest, nlat, nlon)
            r0 = tr * rest
            for k, a in enumerate(aq):
                series[k, r0:r0 + rest, t] = \
                    chunk[:, a["la"], a["lo"]].mean(axis=1)
        if (t + 1) % 25 == 0:
            print(f"  [series {centre}] {t+1}/{n_t} months "
                  f"({(_time.time()-t0)/60:.1f} min)", flush=True)
    times = nc.variables["time"][:]
    nc.close()

    recipe_keys = pd.MultiIndex.from_product(
        [labels[d] for d in dims], names=dims).to_frame(index=False)
    out_dir.mkdir(parents=True, exist_ok=True)
    ds = xr.Dataset(
        {"ewh_anom": (("aquifer", "recipe", "time"), series,
                      {"units": "cm",
                       "long_name": "mask-mean EWH anomaly per recipe",
                       "decision_register": "D018"})},
        coords={"aquifer": [a["idx"] for a in aq],
                "recipe": np.arange(n_recipes),
                "time": ("time", times,
                         {"units": "days since 1970-01-01"})},
        attrs={"centre": centre.upper(), "source_cube": cube_path.name,
               "recipe_key_order": ",".join(dims)})
    ds["aquifer_name"] = ("aquifer", [a["name"] for a in aq])
    path = out_dir / f"aquifer_series_{centre}.nc"
    ds.to_netcdf(path, encoding={"ewh_anom": {"zlib": True, "complevel": 4}})
    recipe_keys.to_csv(out_dir / "recipe_keys.csv", index=False)
    print(f"[series {centre}] wrote {path.name} "
          f"({path.stat().st_size/1e6:.0f} MB, "
          f"{(_time.time()-t0)/60:.1f} min)", flush=True)
    return path


def parity_check(series_nc: pathlib.Path, tier_a_csv: pathlib.Path,
                 tol: float = 1e-4) -> float:
    """Max |series-derived trend - tier_a trend| across (aquifer, recipe)."""
    ds = xr.open_dataset(series_nc)
    y = ds["ewh_anom"].to_numpy()                  # (aq, recipe, t)
    # xarray decodes `time` to datetime64[ns]; .astype(float) would yield
    # NANOSECONDS, not days. Convert through datetime64[D] to get days since
    # epoch, then to years. (The stored series VALUES are unaffected; this
    # only fixes the parity fit's x-axis units.)
    t_yr = ds["time"].to_numpy().astype("datetime64[D]").astype(float) / 365.25
    ta = pd.read_csv(tier_a_csv).sort_values(
        ["aquifer_idx", "recipe_idx_per_centre"])
    n_aq, n_r, _ = y.shape
    ref = ta["trend_cm_yr"].to_numpy().reshape(n_aq, n_r)
    finite = np.isfinite(y)
    tt = np.where(finite, t_yr[None, None, :], np.nan)
    ym = np.nanmean(y, axis=2, keepdims=True)
    tm = np.nanmean(tt, axis=2, keepdims=True)
    num = np.nansum((y - ym) * (tt - tm), axis=2)
    den = np.nansum((tt - tm) ** 2, axis=2)
    slope = num / den
    diff = float(np.nanmax(np.abs(slope - ref)))
    print(f"[series parity] max |series-trend - tier_a-trend| = "
          f"{diff:.2e} cm/yr (tol {tol:g})", flush=True)
    if not (diff < tol):
        raise RuntimeError(f"series/tier_a parity FAILED: {diff:.2e} >= {tol}")
    return diff


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cube", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    ap.add_argument("--parity_csv", type=pathlib.Path, default=None)
    args = ap.parse_args()
    path = extract(args.cube, args.out)
    if args.parity_csv is not None:
        parity_check(path, args.parity_csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
