"""Phase 1 -- truth builder (METHODOLOGY.md section 4.1).

Builds the pilot planet's layered truth in the SPATIAL domain (fairness
principle section 2.2: the truth never shares the recipes' SH basis):

  layer 1a  GLDAS-2.1 Noah rest-of-world hydrology anomaly (soil moisture
            4 layers + SWE + canopy, kg/m^2 -> cm; 0.25 -> 0.5 deg block
            mean; snow-tower cells zeroed, D013)
  layer 1b  stipulated aquifer groundwater: trend + seasonal cycle on the
            73 Tier-1 aquifer masks (truth_spec JSON, blinding-ready D004)

The GIA layer is handled in the SH domain by the forward operator
(02_forward_operator.py) directly from the Bagge et al. (2023) Stokes rates
(D011) -- it never touches the hydrology truth, because GIA is not water and
the scoring target is the WATER trend (D012).

Outputs (datasets/truth/)
  truth_tws_<planet>.nc           float32 (time, lat, lon) cm EWH anomaly
  truth_aquifer_trends_<planet>.csv   per-aquifer OLS trend of the
                                      mask-mean truth series (the tau_i
                                      scoring reference, D012)
  truth_manifest_<planet>.json    sha256 checksums + provenance

Run
---
    python \
        code/forward/01_truth_builder.py [--spec datasets/truth/truth_spec_pilot.json]
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
import pandas as pd
import xarray as xr

PIPELINE_CODE = pathlib.Path(__file__).resolve().parents[3] / "code"
sys.path.insert(0, str(PIPELINE_CODE))

from grace_pipeline import config as cfg  # noqa: E402
from grace_pipeline.io_shm import enumerate_months  # noqa: E402
from grace_pipeline.pipeline_global import GLOBAL_AOI  # noqa: E402

CONTROL_ROOT = pathlib.Path(__file__).resolve().parents[2]
REPO_ROOT = CONTROL_ROOT.parent
TRUTH_DIR = CONTROL_ROOT / "datasets" / "truth"
MASKS_CSV = pathlib.Path(REPO_ROOT / "datasets/jasechko_2024/cell_masks_0p5deg.csv")

ERA5_DIR = cfg.DATA_ROOT / "raw" / "era5_land"
# ERA5-Land background (crossed-truth planets, D039): volumetric soil water
# layers 1-4 x layer thickness (m) -> cm, snow depth and skin reservoir in
# m of water equivalent -> cm. Same composition and factors as the ERA5-Land
# auxiliary of grace_pipeline/gws_aux_aquifer.py.
ERA5_VARS = {"swvl1": 0.07 * 100, "swvl2": 0.21 * 100,
             "swvl3": 0.72 * 100, "swvl4": 1.89 * 100,
             "sd": 100.0, "src": 100.0}

GLDAS_DIR = cfg.RAW_GLDAS
GLDAS_SOIL_VARS = ("SoilMoi0_10cm_inst", "SoilMoi10_40cm_inst",
                   "SoilMoi40_100cm_inst", "SoilMoi100_200cm_inst")
GLDAS_OTHER_VARS = ("SWE_inst", "CanopInt_inst")
_MM_TO_CM = 0.1


def real_sampled_months() -> list:
    """Real CSR months common to BA01 and BB01 over the appendix window
    (D010). The truth exists only where the mission observed."""
    ba = {(m.year, m.month): m
          for m in enumerate_months("CSR", 60, months=cfg.APPENDIX_MONTHS)}
    bb = {(m.year, m.month): m
          for m in enumerate_months("CSR", 96, months=cfg.APPENDIX_MONTHS)}
    return [ba[ym] for ym in sorted(set(ba) & set(bb))]


def load_gldas_month(year: int, month: int) -> xr.DataArray:
    """One GLDAS month -> cm of water on the GLDAS 0.25 deg grid."""
    path = GLDAS_DIR / f"GLDAS_NOAH025_M.A{year:04d}{month:02d}.021.nc4"
    with xr.open_dataset(path) as ds:
        total = sum(ds[v].squeeze(drop=True)
                    for v in GLDAS_SOIL_VARS + GLDAS_OTHER_VARS)
        return (total * _MM_TO_CM).load()


_ERA5_CACHE: dict = {}


def _era5_dataset() -> xr.Dataset:
    """Both ERA5-Land monthly files, concatenated on time (cached)."""
    if "ds" not in _ERA5_CACHE:
        parts = sorted(ERA5_DIR.glob("era5_land_aux_*.nc"))
        if not parts:
            raise FileNotFoundError(f"no ERA5-Land files in {ERA5_DIR}")
        _ERA5_CACHE["ds"] = xr.concat(
            [xr.open_dataset(p) for p in parts], dim="valid_time")
    return _ERA5_CACHE["ds"]


def load_era5_month(year: int, month: int) -> np.ndarray:
    """One ERA5-Land month -> cm of water on the GLOBAL_AOI 0.5 deg grid.

    ERA5 is node-registered (lat 90..-90, lon 0..360) whereas the AOI is
    cell-centred at .25/.75 on -180..180, so a block mean would land on the
    wrong centres; the field is therefore interpolated onto the AOI grid.
    Ocean and out-of-coverage cells carry no land-water anomaly and become 0,
    exactly as for the GLDAS background.
    """
    ds = _era5_dataset()
    sel = ds.sel(valid_time=f"{year:04d}-{month:02d}")
    if "valid_time" in sel.dims:
        sel = sel.isel(valid_time=0)
    total = sum(sel[v].squeeze(drop=True) * f for v, f in ERA5_VARS.items())
    total = total.assign_coords(
        longitude=((total.longitude + 180) % 360) - 180).sortby("longitude")
    total = total.sortby("latitude")
    aoi = GLOBAL_AOI
    out = total.interp(latitude=aoi.lat, longitude=aoi.lon)
    arr = out.to_numpy().astype(np.float64)
    return np.where(np.isfinite(arr), arr, 0.0)


def to_global_half_degree(da: xr.DataArray) -> np.ndarray:
    """0.25 deg GLDAS -> 0.5 deg block mean on the GLOBAL_AOI grid.

    GLDAS is land-only from 60S; ocean and out-of-coverage cells become 0
    (no anomaly) -- the truth planet has no ocean mass variability (D012:
    GSM inputs are AOD-dealiased in reality, so omitting ocean/atmosphere
    signal is the conservative choice).
    """
    coarse = da.coarsen(lat=2, lon=2, boundary="exact").mean()
    aoi = GLOBAL_AOI
    out = coarse.reindex(lat=aoi.lat, lon=aoi.lon, method="nearest",
                         tolerance=0.01)
    arr = out.to_numpy()
    return np.where(np.isfinite(arr), arr, 0.0)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", type=pathlib.Path,
                    default=TRUTH_DIR / "truth_spec_pilot.json")
    args = ap.parse_args()
    spec = json.loads(args.spec.read_text())
    planet = spec["planet_id"]
    aoi = GLOBAL_AOI
    months = real_sampled_months()
    n_t = len(months)
    print(f"[truth {planet}] {n_t} real-sampled months "
          f"({months[0].year}-{months[0].month:02d} .. "
          f"{months[-1].year}-{months[-1].month:02d})")

    # ---- layer 1a: background hydrology (rest of world) -------------------
    # GLDAS-2.1 Noah by default. A crossed-truth planet (D039) builds the
    # background from ERA5-Land instead, so that subtracting the GLDAS-Noah
    # auxiliary in the analysis is no longer an inverse crime and the
    # groundwater subtraction can be scored against a known answer.
    bg = str(spec.get("background_model", "gldas")).lower()
    if bg not in ("gldas", "era5_land"):
        raise ValueError(f"unknown background_model {bg!r}")
    t0 = _time.time()
    stack = np.zeros((n_t, aoi.lat.size, aoi.lon.size), dtype=np.float64)
    for i, m in enumerate(months):
        stack[i] = (to_global_half_degree(load_gldas_month(m.year, m.month))
                    if bg == "gldas" else load_era5_month(m.year, m.month))
        if (i + 1) % 50 == 0:
            print(f"  {bg} {i+1}/{n_t} ({(_time.time()-t0)/60:.1f} min)")
    stack -= stack.mean(axis=0, keepdims=True)          # anomaly vs own mean

    # Snow-tower zeroing (D013): GLDAS-Noah accumulates SWE without bound in
    # permanent-snow cells (same artifact + threshold as the pipeline's
    # gws_correction.py). Those drifts are unphysical mass the real Earth
    # does not produce; left in, they would contaminate neighbouring
    # aquifers through the smoothing filters with a signal that exists in
    # no real gravity record.
    thr = float(spec.get("snow_tower_threshold_cm", 100.0))
    tower = (np.abs(stack) > thr).any(axis=0)
    stack[:, tower] = 0.0
    print(f"[truth {planet}] {bg} done ({(_time.time()-t0)/60:.1f} min); "
          f"{int(tower.sum())} snow-tower cells zeroed (D013)")

    # ---- layer 1b: stipulated aquifer groundwater --------------------------
    masks = pd.read_csv(MASKS_CSV)
    t_ref = np.datetime64(spec["trend_ref_epoch"])
    t_yr = np.array([(m.time_mid - t_ref) / np.timedelta64(1, "D")
                     for m in months]) / 365.25
    month_frac = np.array([m.month for m in months], dtype=float)

    def per_aq(spec_map: dict, aq_idx: int) -> float:
        return float(spec_map.get(str(aq_idx), spec_map["__default__"]))

    gw = np.zeros_like(stack)
    for aq_idx, grp in masks.groupby("aquifer_idx"):
        tau = per_aq(spec["aquifer_trends_cm_yr"], aq_idx)
        amp = per_aq(spec["seasonal_amp_cm"], aq_idx)
        phase = per_aq(spec["seasonal_phase_month"], aq_idx)
        series = tau * t_yr + amp * np.sin(
            2 * np.pi * (month_frac - phase) / 12.0)
        gw[:, grp["lat_idx"].to_numpy(), grp["lon_idx"].to_numpy()] = \
            series[:, None]

    # Edge taper (D027). A hard polygon boundary is unphysical -- real
    # depletion fields decay smoothly into their surroundings -- and the
    # sharp edge injects high-degree signal whose filter-to-filter
    # disagreement does not average down with aquifer size (the +0.17
    # filter-axis area exponent of calibration iterations 1-3). A small
    # Gaussian taper of the groundwater layer removes the artifact without
    # advantaging any recipe; tau_true is measured from the tapered field
    # (D012), so the answer key stays exact.
    taper = float(spec.get("gw_edge_taper_cells", 0.0))
    if taper > 0:
        from scipy.ndimage import gaussian_filter
        for i in range(gw.shape[0]):
            gw[i] = gaussian_filter(gw[i], sigma=taper, mode="nearest")
        print(f"[truth {planet}] groundwater layer edge-tapered "
              f"(Gaussian, {taper} cells; D027)")
    truth = (stack + gw).astype(np.float32)

    # ---- per-aquifer truth trends (the scoring reference, D012) ------------
    rows = []
    tm = t_yr - t_yr.mean()
    for aq_idx, grp in masks.groupby("aquifer_idx"):
        series = truth[:, grp["lat_idx"].to_numpy(),
                       grp["lon_idx"].to_numpy()].mean(axis=1)
        slope = float((tm * (series - series.mean())).sum() / (tm ** 2).sum())
        rows.append({"aquifer_idx": int(aq_idx),
                     "Study_area": grp["Study_area"].iloc[0],
                     "n_cells": len(grp),
                     "tau_true_cm_yr": slope,
                     "tau_stipulated_cm_yr":
                         per_aq(spec["aquifer_trends_cm_yr"], aq_idx)})
    trends = pd.DataFrame(rows)

    # ---- write + freeze -----------------------------------------------------
    TRUTH_DIR.mkdir(parents=True, exist_ok=True)
    nc_path = TRUTH_DIR / f"truth_tws_{planet}.nc"
    times = np.array([m.time_mid for m in months], dtype="datetime64[D]")
    ds = xr.Dataset(
        {"tws_anom": (("time", "lat", "lon"), truth,
                      {"units": "cm", "long_name":
                       "true TWS anomaly (hydrology layers only, no GIA)",
                       "planet_id": planet,
                       "decision_register":
                       "decision register (pre-registration deposit) "
                       "D003 D008 D010 D012 D013"})},
        coords={"time": times, "lat": aoi.lat, "lon": aoi.lon},
        attrs={"title": f"Control Earth truth hydrology -- {planet}",
               "spec_file": args.spec.name,
               "created": dt.datetime.utcnow().isoformat() + "Z",
               "snow_tower_cells_zeroed": int(tower.sum())})
    ds.to_netcdf(nc_path, encoding={"tws_anom": {"zlib": True,
                                                 "complevel": 4}})
    csv_path = TRUTH_DIR / f"truth_aquifer_trends_{planet}.csv"
    trends.to_csv(csv_path, index=False)

    manifest = {
        "planet_id": planet,
        "spec_sha256": hashlib.sha256(
            args.spec.read_bytes()).hexdigest(),
        "truth_nc_sha256": hashlib.sha256(
            nc_path.read_bytes()).hexdigest(),
        "trends_csv_sha256": hashlib.sha256(
            csv_path.read_bytes()).hexdigest(),
        "n_months": n_t, "n_aquifers": int(trends.shape[0]),
        "tau_true_range": [float(trends.tau_true_cm_yr.min()),
                           float(trends.tau_true_cm_yr.max())],
        "created": dt.datetime.utcnow().isoformat() + "Z"}
    (TRUTH_DIR / f"truth_manifest_{planet}.json").write_text(
        json.dumps(manifest, indent=2))
    print(f"[truth {planet}] frozen: {nc_path.name} "
          f"({nc_path.stat().st_size/1e6:.0f} MB); tau_true in "
          f"[{manifest['tau_true_range'][0]:.3f}, "
          f"{manifest['tau_true_range'][1]:.3f}] cm/yr")
    return 0


if __name__ == "__main__":
    sys.exit(main())
