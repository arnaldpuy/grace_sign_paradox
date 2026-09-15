"""Build the mascon comparator cube for the IGP recipe analysis.

The 480-recipe SHM factorial (centre × truncation × filter × GIA × C20 × C30)
deliberately excludes *mascon* products — pre-gridded GRACE solutions in which
the filter / truncation / GIA choices are baked in by the processing centre.
A bibliometric survey found ~46 of 130 IGP-GRACE papers use mascons; this
module tests whether mascons escape the filter-driven spread the SHM recipes
show, or inherit it.

This module treats the three official mascon products — CSR RL0603,
JPL TELLUS RL06.3 v04 CRI, GSFC RL06v2.0 — as **three extra fixed "recipes"**
and writes them onto the SAME IGP AOI grid and the SAME POC / appendix month
axes as the SHM recipe cubes (`grace_ewh_variants_<centre>.nc`). The R
analysis layer (`load_mascon_recipes_fun.R` → `mascon_well_skill_fun.R`) then
runs them through the identical well-skill + trend machinery.

Decisions, all explicit:
  * pre-gridded mascon products instead of the SH-domain factorial — that
    is the whole point of the comparator (the SHM cube is the sensitivity
    arm; this is the "what does the literature actually use" arm).
  * THREE SEPARATE products, not an ensemble mean — we need to place each
    mascon individually against the SHM cloud, not collapse them.
  * bilinear interpolation from native grids (CSR 0.25°, JPL/GSFC 0.5°) to
    the AOI 0.25° grid. Mascons retain their native ~300 km support; the
    0.25° spacing is cosmetic and recorded honestly in the attributes.
  * monthly, snapped to month-15 — matches the SHM cubes; NO monthly→daily
    interpolation (the cube_pipeline harmoniser does that for a different
    consumer; see code/cube_pipeline/harmonise_grace_mascon.py).
  * NO re-baselining. Each mascon is an anomaly vs its own centre baseline;
    a constant offset changes neither an OLS trend nor a deseasonalised
    correlation, and the SHM pipeline does not cross-baseline recipes
    either. The native baseline is recorded as a NetCDF attribute.

Output: grace_mascon_variants.nc          (POC, 89 months)
        grace_mascon_variants_appendix.nc (appendix, 272 months)
both with dims (time, product, lat, lon), product = [CSR, JPL, GSFC].
"""
from __future__ import annotations

import datetime as dt
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from . import config as cfg

# --------------------------------------------------------------------------
# Raw mascon product locations (downloaded once by the cube pipeline's
# download_grace_mascon.py — the same three files serve both consumers).
# --------------------------------------------------------------------------
RAW_MASCON = cfg.DATA_ROOT / "raw" / "grace_mascon"
PRODUCT_ORDER = ("CSR", "JPL", "GSFC")

# Each mascon centre's GIA model is baked in by the processing centre. All
# three current products use ICE6G-D — this is what makes them comparable to
# the GIA-CORRECTED sub-cloud of the SHM factorial, not the full factorial.
MASCON_GIA = {"CSR": "ICE6G-D", "JPL": "ICE6G-D", "GSFC": "ICE6G-D"}
MASCON_LABEL = {"CSR": "CSR-M", "JPL": "JPL-M", "GSFC": "GSFC-M"}
# Native anomaly baseline each centre subtracts (for the provenance record).
MASCON_BASELINE = {
    "CSR": "2004.000-2009.999",
    "JPL": "2004.000-2009.999",
    "GSFC": "2004.000-2009.999",
}


def _csr_path() -> Path:
    return (RAW_MASCON / "csr" /
            "CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc")


def _jpl_path() -> Path:
    # Glob rather than hard-code: the JPL filename carries the data span and
    # changes with each release/refresh.
    hits = sorted((RAW_MASCON / "jpl").glob("GRCTellus.JPL.*.nc"))
    if not hits:
        raise FileNotFoundError(f"no JPL mascon file under {RAW_MASCON/'jpl'}")
    return hits[-1]


def _gsfc_path() -> Path:
    hits = sorted((RAW_MASCON / "gsfc").glob("gsfc.glb_*.nc"))
    if not hits:
        raise FileNotFoundError(f"no GSFC mascon file under {RAW_MASCON/'gsfc'}")
    return hits[-1]


# --------------------------------------------------------------------------
# Per-product loaders duplicated here so grace_pipeline stays self-contained
# (not coupled to the cube pipeline's config/AOI).
# --------------------------------------------------------------------------
def _load_csr() -> xr.DataArray:
    """CSR mascons. The time attribute uses a capital-U `Units`, which breaks
    xarray's default CF decoding; decode the days-since-2002-01-01 manually."""
    ds = xr.open_dataset(_csr_path(), decode_cf=False)
    epoch = np.datetime64("2002-01-01")
    time = epoch + (ds.time.values * np.timedelta64(1, "D")).astype(
        "timedelta64[D]")
    da = ds["lwe_thickness"].assign_coords(time=("time", time))
    da.attrs["centre"] = "CSR"
    return da


def _load_jpl() -> xr.DataArray:
    ds = xr.open_dataset(_jpl_path())
    da = ds["lwe_thickness"]
    # JPL ships a per-cell scale_factor restoring grid-cell equivalent depth
    # from the mascon-block average; without it the values are block means.
    if "scale_factor" in ds.data_vars:
        da = da * ds["scale_factor"]
    da.attrs["centre"] = "JPL"
    return da


def _load_gsfc() -> xr.DataArray:
    ds = xr.open_dataset(_gsfc_path())
    da = ds["lwe_thickness"]
    da.attrs["centre"] = "GSFC"
    return da


_LOADERS = {"CSR": _load_csr, "JPL": _load_jpl, "GSFC": _load_gsfc}


# --------------------------------------------------------------------------
# Regrid + temporal harmonisation
# --------------------------------------------------------------------------
def _to_aoi_monthly(da: xr.DataArray, aoi: cfg.AOI) -> xr.DataArray:
    """Clip to a generous bbox around the AOI, bilinearly regrid to the AOI
    0.25° grid, resample to monthly mean and snap each stamp to month-15.

    Mascon files are stored on a 0–360 longitude axis, so the bbox subset
    and the interpolation target are both expressed in 0–360 (`aoi.lon360`,
    a no-op for the IGP at 73–88E). The result's lon coordinate is then
    re-assigned to the AOI's native convention so the mascon cube aligns
    byte-for-byte with the SHM recipe cubes."""
    da = da.sortby("lat").sortby("lon")
    # Subset a bbox strictly containing the AOI BEFORE interpolation so the
    # source spans the target interval. CV's box (−122..−118E → 237..242E)
    # does not cross the 0/360 seam, so the wrapped bounds stay ascending.
    lon_lo = (aoi.lon_min - 1.0) % 360.0
    lon_hi = (aoi.lon_max + 1.0) % 360.0
    da = da.sel(lat=slice(aoi.lat_min - 1.0, aoi.lat_max + 1.0),
                lon=slice(lon_lo, lon_hi))
    da = da.interp(lat=aoi.lat, lon=aoi.lon360, method="linear")
    # Re-express lon in the AOI's native convention (matches the SHM cubes).
    da = da.assign_coords(lon=("lon", aoi.lon))

    # Monthly mean (mascon stamps are already ~monthly; this absorbs the
    # slight per-product stamping offsets), then snap to month-15.
    monthly = da.resample(time="MS").mean()
    new_time = pd.to_datetime(
        [dt.date(t.year, t.month, 15) for t in pd.to_datetime(monthly.time.values)])
    monthly = monthly.assign_coords(time=("time", new_time.values))
    return monthly


def _reindex_to_window(monthly: xr.DataArray, window: str) -> xr.DataArray:
    """Place the monthly series onto the exact POC / appendix month axis,
    NaN-filling any month the product does not cover (e.g. GSFC stops at
    2025-11, so it is NaN for the final POC month)."""
    months = cfg.months_for_window(window)
    target = pd.to_datetime([dt.date(y, m, 15) for y, m in months])
    return monthly.reindex(time=target)


# --------------------------------------------------------------------------
# Build + write
# --------------------------------------------------------------------------
def _git_commit() -> str:
    try:
        out = subprocess.check_output(
            ["git", "-C", str(cfg.REPO_ROOT),
             "rev-parse", "HEAD"], stderr=subprocess.DEVNULL)
        return out.decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return "unknown"


def build_window(window: str, region: str = "igp") -> xr.Dataset:
    """Build the (time, product, lat, lon) mascon cube for one window."""
    aoi = cfg.aoi_for_region(region)
    region_label = cfg.REGION_LABEL[region]
    print(f"[mascon] building region={region!r} window={window!r}")
    per_product = []
    for name in PRODUCT_ORDER:
        raw = _LOADERS[name]()
        monthly = _to_aoi_monthly(raw, aoi)
        windowed = _reindex_to_window(monthly, window)
        n_ok = int(np.isfinite(windowed).any(dim=("lat", "lon")).sum())
        print(f"  {name:5s}: {n_ok}/{windowed.sizes['time']} months with data")
        per_product.append(windowed.expand_dims(product=[name]))

    cube = xr.concat(per_product, dim="product")
    cube.name = "ewh_anom"
    cube.attrs["units"] = "cm"
    cube.attrs["long_name"] = (
        "GRACE mascon liquid-water-equivalent thickness anomaly, "
        f"regridded to the {region_label} AOI")
    ds = xr.Dataset({"ewh_anom": cube})
    ds.attrs.update({
        "title": f"GRACE mascon comparator cube — {region_label} — {window}",
        "region": region,
        "window": window,
        "products": ", ".join(MASCON_LABEL[p] for p in PRODUCT_ORDER),
        "csr_file": _csr_path().name,
        "jpl_file": _jpl_path().name,
        "gsfc_file": _gsfc_path().name,
        "gia_models": "; ".join(f"{MASCON_LABEL[p]}={MASCON_GIA[p]}"
                                for p in PRODUCT_ORDER),
        "native_baselines": "; ".join(f"{MASCON_LABEL[p]}={MASCON_BASELINE[p]}"
                                      for p in PRODUCT_ORDER),
        "regrid_method": "bilinear (xarray.interp) to AOI 0.25 deg grid",
        "support_kernel_note": (
            "mascons retain their native ~300 km effective support; the "
            "0.25 deg grid spacing is cosmetic, not resolved information"),
        "baseline_note": (
            "each product is an anomaly vs its centre's own baseline; NOT "
            "re-baselined here — a constant offset changes neither an OLS "
            "trend nor a deseasonalised correlation"),
        "aoi_bbox": (f"{aoi.lat_min}N-{aoi.lat_max}N, "
                     f"{aoi.lon_min}E-{aoi.lon_max}E"),
        "aoi_resolution_deg": aoi.res,
        "code_git_commit": _git_commit(),
        "processing_date": dt.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "consumer": "analysis R analysis — mascon comparator (Result 14)",
    })
    return ds


def write_window(window: str, region: str = "igp",
                 out_dir: Path = cfg.OUT_DIR) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    ds = build_window(window, region=region)
    region_suffix = cfg.REGION_SUFFIX[region]
    window_suffix = "_appendix" if window == "appendix" else ""
    fname = f"grace_mascon_variants{region_suffix}{window_suffix}.nc"
    path = out_dir / fname
    ds.to_netcdf(path)
    print(f"[mascon] wrote {path}  "
          f"(time={ds.sizes['time']}, product={ds.sizes['product']}, "
          f"lat={ds.sizes['lat']}, lon={ds.sizes['lon']})")
    return path


def main(argv: list[str] | None = None) -> None:
    import argparse
    ap = argparse.ArgumentParser(
        description="Build the GRACE mascon comparator cube(s).")
    ap.add_argument(
        "--region", default="igp", choices=tuple(cfg.AOI_REGISTRY),
        help="Study aquifer: 'igp' (Indo-Gangetic Plain) or 'cv' "
             "(California Central Valley).")
    args = ap.parse_args(argv)
    for window in ("poc", "appendix"):
        write_window(window, region=args.region)


if __name__ == "__main__":
    main()
