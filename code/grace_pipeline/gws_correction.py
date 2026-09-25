"""Build the GRACE-GWS correction cube: the non-groundwater storage that
must be subtracted from GRACE TWS to isolate groundwater storage (GWS).

GRACE TWS = groundwater + soil moisture + snow water equivalent + canopy
water + surface water. This module builds the recipe-INDEPENDENT correction
field

    gws_correction = soil moisture + SWE + canopy  [ + surface water ]

on the AOI grid and the project month-15 time axis, as a cm-water-equivalent
anomaly. The R analysis subtracts it from every recipe
(`apply_gws_correction_fun`); because the field is the same for all 480
recipes, the per-cell inter-recipe spread is preserved exactly — the
per-axis variance decomposition is identical to TWS mode cell-for-cell
(the aggregate median differs only in the 4th decimal, solely from the
few snow-tower cells this cube NaN-masks; see _MAX_PHYSICAL_ANOMALY_CM).
Baseline- and threshold-derived metrics (dominance R, hotspot Jaccard,
flip rates) shift modestly. Either way the subtraction re-bases the
analysis from TWS to GWS, which is what makes the in-situ well comparison
physically valid.

Sources:
  * soil moisture (4 layers, 0-200 cm), SWE, canopy — GLDAS-2.1 Noah
    monthly 0.25° (`download_gldas.py`). GLDAS is on a -180..180 longitude
    axis, so it is regridded against `aoi.lon` (NATIVE) — the OPPOSITE of
    the 0..360 mascon files. Units kg/m^2 == mm; divided by 10 -> cm.
  * surface water — CDEC major-reservoir storage (`download_cdec_reservoirs`),
    Central Valley only and OPTIONAL: if cdec_reservoirs.csv is absent the
    cube is built soil+snow+canopy only and the omission is recorded in the
    attributes. The IGP has no surface-water term by design (see
    config.SURFACE_WATER) — it is documented as a stated limitation.

Anomaly baseline: the correction is expressed as an anomaly relative to the
mean over the analysis window being built (POC or appendix). Because the
trap metrics and OLS trends are invariant to a constant offset, the choice
of baseline only shifts the GWS *level* (and dominance R's median
denominator), not the spread or the trends — but it is recorded in the
`anomaly_baseline` attribute regardless.

Output: gws_correction.nc / gws_correction_appendix.nc / gws_correction_cv.nc
        / gws_correction_cv_appendix.nc  (dims time x lat x lon, var
        `gws_correction`, units cm).
"""
from __future__ import annotations

import datetime as dt
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from . import config as cfg
from .download_gldas import GLDAS_SHORT_NAME, GLDAS_VERSION

# GLDAS-2.1 Noah variables summed into the non-groundwater storage term.
# All are kg/m^2 (== mm of water).
GLDAS_SOIL_VARS = ("SoilMoi0_10cm_inst", "SoilMoi10_40cm_inst",
                   "SoilMoi40_100cm_inst", "SoilMoi100_200cm_inst")
GLDAS_OTHER_VARS = ("SWE_inst", "CanopInt_inst")
_MM_TO_CM = 0.1
_ACRE_FOOT_M3 = 1233.4818375
_EARTH_R_M = 6_371_000.0

# GLDAS-Noah "snow tower" artifact: permanent-snow cells accumulate SWE without
# bound (high-Himalaya cells drift to +/-400 cm vs <=50 cm elsewhere). Cells whose
# correction anomaly ever exceeds the bound are masked NaN; count + locations are
# recorded in the `snow_tower_cells_masked` attribute.
_MAX_PHYSICAL_ANOMALY_CM = 100.0


# --------------------------------------------------------------------------
# GLDAS soil moisture + SWE + canopy
# --------------------------------------------------------------------------
def _load_gldas_non_gw() -> xr.DataArray:
    """Open all GLDAS monthly granules, sum soil moisture + SWE + canopy,
    return a global (time, lat, lon) DataArray in cm."""
    files = sorted(cfg.RAW_GLDAS.glob("GLDAS_NOAH025_M*.nc*"))
    if not files:
        raise FileNotFoundError(
            f"no GLDAS granules under {cfg.RAW_GLDAS} — run download_gldas.py")
    ds = xr.open_mfdataset(files, combine="by_coords")
    want = GLDAS_SOIL_VARS + GLDAS_OTHER_VARS
    missing = [v for v in want if v not in ds.data_vars]
    if missing:
        raise KeyError(f"GLDAS granules missing variables: {missing}")
    total_mm = sum(ds[v].fillna(0.0) for v in want)
    da = (total_mm * _MM_TO_CM).rename("non_gw_storage_cm")
    da.attrs["units"] = "cm"
    da.attrs["components"] = ", ".join(want)
    return da


def _to_aoi_monthly(da: xr.DataArray, aoi: cfg.AOI) -> xr.DataArray:
    """Clip to a generous bbox around the AOI, bilinearly regrid to the AOI
    0.25 deg grid, resample to monthly mean, snap each stamp to month-15.

    GLDAS is stored on a -180..180 longitude axis, which already matches the
    AOI's native `aoi.lon` convention (IGP 73..88E, CV -122..-118E). So we
    regrid against `aoi.lon` (NATIVE), NOT `aoi.lon360` — the OPPOSITE of
    mascon.py, where the products were 0..360."""
    da = da.sortby("lat").sortby("lon")
    da = da.sel(lat=slice(aoi.lat_min - 1.0, aoi.lat_max + 1.0),
                lon=slice(aoi.lon_min - 1.0, aoi.lon_max + 1.0))
    da = da.interp(lat=aoi.lat, lon=aoi.lon, method="linear")
    monthly = da.resample(time="MS").mean()
    new_time = pd.to_datetime(
        [dt.date(t.year, t.month, 15)
         for t in pd.to_datetime(monthly.time.values)])
    monthly = monthly.assign_coords(time=("time", new_time.values))
    return monthly


# --------------------------------------------------------------------------
# Surface water (CDEC reservoirs) — Central Valley, optional
# --------------------------------------------------------------------------
def _aoi_area_m2(aoi: cfg.AOI) -> float:
    """Total AOI surface area (m^2), cos(lat)-weighted 0.25 deg cells."""
    cell_deg = aoi.res
    dlat_m = np.deg2rad(cell_deg) * _EARTH_R_M
    lat = aoi.lat
    dlon_m = np.deg2rad(cell_deg) * _EARTH_R_M * np.cos(np.deg2rad(lat))
    # one row of cells per latitude, len(aoi.lon) columns
    return float(np.sum(dlat_m * dlon_m) * len(aoi.lon))


def _surface_water_cv(aoi: cfg.AOI) -> xr.DataArray | None:
    """Monthly CV reservoir storage as a spatially-uniform cm series, or
    None if the CDEC file is absent. Reservoirs are sub-grid points, so the
    honest representation is an AOI-mean term applied uniformly."""
    csv_path = cfg.RAW_CDEC / "cdec_reservoirs.csv"
    if not csv_path.exists():
        print("  surface water: cdec_reservoirs.csv absent — "
              "building soil+snow+canopy only")
        return None
    df = pd.read_csv(csv_path)
    # CDEC columns vary by endpoint; normalise the date + value columns.
    date_col = next((c for c in ("date", "DATE TIME", "obsDate", "datetime")
                     if c in df.columns), None)
    val_col = next((c for c in ("value", "VALUE", "obsValue")
                    if c in df.columns), None)
    if date_col is None or val_col is None:
        raise KeyError(f"unrecognised CDEC columns: {list(df.columns)}")
    df["date"] = pd.to_datetime(df[date_col], errors="coerce")
    df["value"] = pd.to_numeric(df[val_col], errors="coerce")
    df = df.dropna(subset=["date", "value"])
    df = df[df["value"] >= 0]                       # CDEC missing flags < 0
    # total storage across reservoirs per month (acre-feet)
    df["ym"] = df["date"].dt.to_period("M")
    total_af = df.groupby("ym")["value"].sum()
    total_m3 = total_af * _ACRE_FOOT_M3
    storage_cm = (total_m3 / _aoi_area_m2(aoi)) * 100.0   # m -> cm
    times = pd.to_datetime([dt.date(p.year, p.month, 15)
                            for p in total_af.index])
    return xr.DataArray(storage_cm.values, coords={"time": times},
                        dims="time", name="surface_water_cm")


# --------------------------------------------------------------------------
# Window assembly
# --------------------------------------------------------------------------
def _reindex_to_window(da: xr.DataArray, window: str) -> xr.DataArray:
    """Place a monthly series onto the exact POC / appendix month-15 axis."""
    months = cfg.months_for_window(window)
    target = pd.to_datetime([dt.date(y, m, 15) for y, m in months])
    return da.reindex(time=target)


def _git_commit() -> str:
    try:
        out = subprocess.check_output(
            ["git", "-C", str(cfg.REPO_ROOT),
             "rev-parse", "HEAD"], stderr=subprocess.DEVNULL)
        return out.decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return "unknown"


def build_window(window: str, region: str = "igp") -> xr.Dataset:
    """Build the (time, lat, lon) GWS-correction cube for one region/window."""
    aoi = cfg.aoi_for_region(region)
    region_label = cfg.REGION_LABEL[region]
    print(f"[gws] building region={region!r} window={window!r}")

    # soil moisture + SWE + canopy, regridded + monthly
    non_gw = _to_aoi_monthly(_load_gldas_non_gw(), aoi)
    non_gw = _reindex_to_window(non_gw, window)
    n_ok = int(np.isfinite(non_gw).any(dim=("lat", "lon")).sum())
    print(f"  GLDAS soil+snow+canopy: {n_ok}/{non_gw.sizes['time']} months")

    components = ["GLDAS-2.1 soil moisture (0-200 cm)", "SWE", "canopy"]
    surface = None
    if cfg.SURFACE_WATER.get(region, False):
        sw = _surface_water_cv(aoi)
        if sw is not None:
            sw = _reindex_to_window(sw, window)
            surface = sw
            components.append("CDEC reservoir surface water")
            print(f"  + CDEC surface water: "
                  f"{int(np.isfinite(sw).sum())}/{sw.sizes['time']} months")

    # total non-groundwater storage (cm)
    total = non_gw
    if surface is not None:
        total = total + surface           # broadcasts (time,) over (lat,lon)

    # express as an anomaly vs the analysis-window mean
    baseline = total.mean(dim="time", skipna=True)
    corr = (total - baseline).rename("gws_correction")

    # mask GLDAS-Noah "snow tower" cells (see _MAX_PHYSICAL_ANOMALY_CM)
    snow_tower = (np.abs(corr) > _MAX_PHYSICAL_ANOMALY_CM).any(dim="time")
    n_masked = int(snow_tower.sum())
    masked_locs = ""
    if n_masked:
        li, loi = np.where(snow_tower.values)
        masked_locs = "; ".join(
            f"{float(corr.lat[a]):.3f}N,{float(corr.lon[b]):.3f}E"
            for a, b in zip(li, loi))
        print(f"  masked {n_masked} GLDAS snow-tower cell(s): {masked_locs}")
        corr = corr.where(~snow_tower)

    corr.attrs["units"] = "cm"
    corr.attrs["long_name"] = (
        "non-groundwater storage anomaly (soil moisture + SWE + canopy"
        + (" + surface water" if surface is not None else "")
        + ") — subtract from GRACE TWS to obtain GRACE-GWS")

    ds = xr.Dataset({"gws_correction": corr})
    ds.attrs.update({
        "title": f"GRACE-GWS correction cube — {region_label} — {window}",
        "region": region,
        "window": window,
        "components": "; ".join(components),
        "surface_water_included": str(surface is not None),
        "surface_water_note": (
            "CDEC major-reservoir storage, applied as a spatially-uniform "
            "AOI-mean term" if surface is not None else
            "no surface-water term — stated limitation (IGP by design; CV "
            "only if CDEC was unavailable at build time)"),
        "gldas_product": f"{GLDAS_SHORT_NAME} v{GLDAS_VERSION}",
        "gldas_lon_convention": "-180..180 (regridded against aoi.lon native)",
        "units_conversion": "GLDAS kg/m^2 == mm; divided by 10 -> cm",
        "snow_tower_cells_masked": (
            f"{n_masked} cell(s) with |anomaly| > "
            f"{_MAX_PHYSICAL_ANOMALY_CM:.0f} cm masked (GLDAS-Noah unbounded-"
            f"SWE artifact): {masked_locs}" if n_masked else "none"),
        "anomaly_baseline": (
            f"mean over the {window} analysis window; trap metrics and OLS "
            "trends are invariant to this constant, only the GWS level shifts"),
        "regrid_method": "bilinear (xarray.interp) to AOI 0.25 deg grid",
        "aoi_bbox": (f"{aoi.lat_min}N-{aoi.lat_max}N, "
                     f"{aoi.lon_min}E-{aoi.lon_max}E"),
        "aoi_resolution_deg": aoi.res,
        "code_git_commit": _git_commit(),
        "processing_date": dt.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "consumer": ("analysis R analysis — apply_gws_correction_fun "
                     "subtracts this from every recipe (GWS mode)"),
    })
    return ds


def write_window(window: str, region: str = "igp",
                 out_dir: Path = cfg.OUT_DIR) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    ds = build_window(window, region=region)
    region_suffix = cfg.REGION_SUFFIX[region]
    window_suffix = "_appendix" if window == "appendix" else ""
    fname = f"gws_correction{region_suffix}{window_suffix}.nc"
    path = out_dir / fname
    ds.to_netcdf(path)
    print(f"[gws] wrote {path}  "
          f"(time={ds.sizes['time']}, lat={ds.sizes['lat']}, "
          f"lon={ds.sizes['lon']})")
    return path


def main(argv: list[str] | None = None) -> None:
    import argparse
    ap = argparse.ArgumentParser(
        description="Build the GRACE-GWS correction cube(s).")
    ap.add_argument(
        "--region", default="igp", choices=tuple(cfg.AOI_REGISTRY),
        help="Study aquifer: 'igp' (Indo-Gangetic Plain) or 'cv' "
             "(California Central Valley).")
    args = ap.parse_args(argv)
    for window in ("poc", "appendix"):
        write_window(window, region=args.region)


if __name__ == "__main__":
    main()
