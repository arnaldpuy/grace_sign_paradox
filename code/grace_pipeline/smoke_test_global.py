"""Tier-A smoke test: per-field byte cost for one global EWH synthesis.

Validates the on-disk storage projection for Strategy B (full 1,920-recipe
materialisation everywhere) before launching the multi-day full synthesis.

What it does
------------
1. Synthesises ONE EWH field for ONE recipe x ONE month x ONE centre
   on a global lat/lon grid at user-chosen resolution (default 0.5 deg).
2. Writes it to NetCDF four ways:
     A. float64 no compression  (current pipeline default)
     B. float64 zlib deflate L4
     C. float32 no compression
     D. float32 zlib deflate L4
3. Reports each on-disk size + the projection to the full Tier-A cube:
     n_fields = 640 variants/centre * 89 months * 3 centres = 170_880.

Usage
-----
    cd <repository>/code
    python -m grace_pipeline.smoke_test_global \\
        --centre CSR --year 2024 --month 6 --grid_res 0.5 \\
        --out_dir /tmp/tier_a_smoke

Outputs
-------
    {out_dir}/smoke_f64_raw.nc       (variant A)
    {out_dir}/smoke_f64_deflate.nc   (variant B)
    {out_dir}/smoke_f32_raw.nc       (variant C)
    {out_dir}/smoke_f32_deflate.nc   (variant D)
    {out_dir}/storage_smoke_test.csv  one row per variant
"""
from __future__ import annotations

import argparse
import sys
import time as _time
import csv
from dataclasses import replace
from pathlib import Path

import numpy as np
import xarray as xr

from . import config as cfg
from .corrections_sh import replace_low_degrees
from .filters_sh import apply_filter
from .geocenter import apply_geocenter, load_geocenter
from .io_aux import load_tn14
from .io_shm import enumerate_months, load_shm
from .synthesis import synthesise_to_aoi


# Total field count for the planned Tier-A cube:
N_FIELDS_TOTAL = 640 * 89 * 3   # variants_per_centre x POC_months x centres
                                # = 170,880


def make_global_aoi(res: float) -> cfg.AOI:
    """Global AOI dataclass. Latitudes capped at +-85 to avoid polar
    pathologies in SH evaluation; aquifers are well within [-60, +60]."""
    return cfg.AOI(lat_min=-85.0, lat_max=85.0,
                   lon_min=-180.0, lon_max=180.0, res=res)


def synthesise_one_global(centre: str, year: int, month: int,
                          grid_res: float, filt: str = "G300",
                          truncation: int = 60) -> np.ndarray:
    """One-recipe global synthesis. Reuses the production pipeline modules."""
    aoi = make_global_aoi(res=grid_res)
    nlat, nlon = aoi.lat.size, aoi.lon.size
    print(f"[smoke] global AOI: {nlat}x{nlon} = {nlat*nlon:,} cells "
          f"({grid_res} deg)")

    months = enumerate_months(centre, truncation, months=[(year, month)])
    if not months:
        raise SystemExit(
            f"No SHM file for {centre} lmax={truncation} {year}-{month:02d}; "
            "pick a different month.")
    m = months[0]
    print(f"[smoke] {centre} {year}-{month:02d} ({m.file_path.name})")

    tn14 = load_tn14()
    geocenter_series = load_geocenter("none")

    t0 = _time.time()
    harm_raw  = load_shm(m.file_path, lmax=truncation)
    t1 = _time.time()
    harm_corr = replace_low_degrees(harm_raw, tn14, m, "replaced", "replaced")
    harm_gc   = apply_geocenter(harm_corr, geocenter_series, m)
    harm_filt = apply_filter(harm_gc, filt, truncation)
    t2 = _time.time()
    ewh = synthesise_to_aoi(harm_filt, aoi, truncation)
    t3 = _time.time()
    print(f"[smoke] load_shm:        {t1-t0:6.2f}s")
    print(f"[smoke] SH-domain corrs: {t2-t1:6.2f}s")
    print(f"[smoke] SH->grid synth:  {t3-t2:6.2f}s   (this is the dominant per-field cost)")
    return ewh, aoi


def write_and_measure(ewh: np.ndarray, aoi: cfg.AOI, dtype: str,
                       zlib: bool, complevel: int, path: Path) -> int:
    ds = xr.Dataset(
        data_vars={"ewh": (("lat", "lon"), ewh.astype(dtype))},
        coords={"lat": aoi.lat, "lon": aoi.lon},
        attrs={"smoke_test": "true",
                "dtype": dtype,
                "compression": f"zlib_complevel_{complevel}" if zlib else "none"})
    enc = {"ewh": {"zlib": zlib, "complevel": complevel}} if zlib else None
    ds.to_netcdf(path, encoding=enc)
    return path.stat().st_size


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--centre", required=True, choices=cfg.CENTRES)
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--month", type=int, required=True)
    ap.add_argument("--grid_res", type=float, default=0.5,
                     help="Grid resolution in degrees (default 0.5)")
    ap.add_argument("--filter", default="G300",
                     help="Filter to apply (default G300; light "
                          "smoother is the fastest test).")
    ap.add_argument("--trunc", type=int, default=60,
                     help="SH truncation (default 60)")
    ap.add_argument("--out_dir", type=Path, default=Path("/tmp/tier_a_smoke"))
    args = ap.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ewh, aoi = synthesise_one_global(args.centre, args.year, args.month,
                                       args.grid_res, args.filter, args.trunc)

    print("\n[smoke] writing 4 variants:")
    variants = [
        ("f64_raw",     "float64", False, 0),
        ("f64_deflate", "float64", True,  4),
        ("f32_raw",     "float32", False, 0),
        ("f32_deflate", "float32", True,  4),
    ]
    rows = []
    for label, dt, zlib, lvl in variants:
        path = args.out_dir / f"smoke_{label}.nc"
        size = write_and_measure(ewh, aoi, dt, zlib, lvl, path)
        proj_bytes = size * N_FIELDS_TOTAL
        proj_gb = proj_bytes / 1e9
        print(f"  {label:14s}: per-field {size/1e3:8.1f} kB   "
               f"x {N_FIELDS_TOTAL:,} -> {proj_gb:7.1f} GB total")
        rows.append({
            "variant": label,
            "dtype": dt,
            "compression": f"zlib_L{lvl}" if zlib else "none",
            "per_field_bytes": size,
            "per_field_kB": round(size / 1e3, 2),
            "n_fields_total": N_FIELDS_TOTAL,
            "projected_GB_total": round(proj_gb, 2),
            "grid_res_deg": args.grid_res,
            "global_cells_nlat": aoi.lat.size,
            "global_cells_nlon": aoi.lon.size,
            "test_field_centre": args.centre,
            "test_field_year": args.year,
            "test_field_month": args.month,
            "test_field_filter": args.filter,
            "test_field_trunc": args.trunc,
        })

    csv_path = args.out_dir / "storage_smoke_test.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n[smoke] CSV: {csv_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
