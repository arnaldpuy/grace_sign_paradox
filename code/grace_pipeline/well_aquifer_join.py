"""Spatial join: Jasechko monitoring wells -> 73-cohort aquifer polygons.

For each well in `AnnualDepthToGroundwater.csv` (Jasechko's redistributable
subset), determine which (if any) cohort polygon contains it. Output a CSV
that downstream well-skill machinery uses to filter wells per aquifer.

Output: datasets/output/tier_a/wells_in_cohort.csv
   columns: StnID, Lat, Lon, aquifer_idx, Study_area, lat_idx, lon_idx
   (lat_idx, lon_idx are the global 0.5deg grid indices for the cell
    containing the well; used by the well-vs-GRACE correlation step.)

Also outputs per-aquifer well counts so we can confirm the 56-aquifer
'>=10 wells' subset from the scouting day matches what the actual
spatial join produces.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd

from . import config as cfg


_CODE_REPO = cfg.REPO_ROOT
_JAS_DIR = _CODE_REPO / "datasets" / "jasechko_2024"

DEFAULT_SHAPEFILE = _JAS_DIR / "aquifers_shp" / "jasechko_et_al_2024_aquifers.shp"
DEFAULT_WELLS_CSV = _JAS_DIR / "groundwater_levels" / "AnnualDepthToGroundwater.csv"
DEFAULT_COHORT_CSV = _JAS_DIR / "cohort_50K_n73.csv"
DEFAULT_OUT_CSV = _CODE_REPO / "datasets" / "output" / "tier_a" / "wells_in_cohort.csv"


def join_wells_to_cohort(wells_csv: Path = DEFAULT_WELLS_CSV,
                          shapefile: Path = DEFAULT_SHAPEFILE,
                          cohort_csv: Path = DEFAULT_COHORT_CSV,
                          grid_res: float = 0.5) -> pd.DataFrame:
    # 1. Load wells (~170k records) -> deduplicate to unique (StnID, Lat, Lon)
    print(f"[wells] loading {wells_csv.name} ...")
    wells = pd.read_csv(wells_csv, usecols=["StnID", "Lat", "Lon"])
    wells = wells.drop_duplicates("StnID").reset_index(drop=True)
    print(f"[wells] {len(wells):,} unique wells "
          f"(lat {wells.Lat.min():+.1f}..{wells.Lat.max():+.1f}, "
          f"lon {wells.Lon.min():+.1f}..{wells.Lon.max():+.1f})")

    # 2. Filter cohort polygons
    shp = gpd.read_file(shapefile)
    cohort = pd.read_csv(cohort_csv)
    shp["combo_key"] = [
        a if (b is None or (isinstance(b, float) and np.isnan(b))
              or str(b).strip() in ("", "-"))
        else f"{a} ({b})"
        for a, b in zip(shp["Aquifer"], shp["Broader"])
    ]
    shp = shp[shp["combo_key"].str.lower().isin(
        cohort["Study_area"].str.lower())].copy()
    shp = shp.merge(
        cohort[["Study_area"]].assign(_k=cohort["Study_area"].str.lower()),
        left_on=shp["combo_key"].str.lower(),
        right_on="_k", how="left")
    shp = shp.sort_values(
        "Study_area",
        key=lambda s: s.map({n: i for i, n in enumerate(
            cohort["Study_area"])})).reset_index(drop=True)
    shp["aquifer_idx"] = range(len(shp))
    print(f"[wells] cohort polygons matched: {len(shp)} of 73")

    # 3. Build well GeoDataFrame
    wells_gdf = gpd.GeoDataFrame(
        wells, geometry=gpd.points_from_xy(wells.Lon, wells.Lat),
        crs="EPSG:4326")

    # 4. Spatial join
    print(f"[wells] joining {len(wells_gdf):,} wells against "
          f"{len(shp)} polygons ...")
    joined = gpd.sjoin(
        wells_gdf, shp[["aquifer_idx", "Study_area", "geometry"]],
        how="inner", predicate="within")
    print(f"[wells] {len(joined):,} wells inside cohort polygons "
          f"({100 * len(joined) / len(wells_gdf):.1f}% of unique wells)")

    # 5. Compute global-grid cell index for each well
    lats = np.arange(-85.0 + grid_res / 2, 85.0, grid_res)
    lons = np.arange(-180.0 + grid_res / 2, 180.0, grid_res)
    # Closest lat/lon idx for each well
    lat_idx = np.searchsorted(lats - grid_res / 2, joined.Lat.values) - 1
    lat_idx = np.clip(lat_idx, 0, len(lats) - 1)
    lon_idx = np.searchsorted(lons - grid_res / 2, joined.Lon.values) - 1
    lon_idx = np.clip(lon_idx, 0, len(lons) - 1)

    out = joined[["StnID", "Lat", "Lon", "aquifer_idx",
                   "Study_area"]].copy()
    out["lat_idx"] = lat_idx
    out["lon_idx"] = lon_idx
    out = out.sort_values(["aquifer_idx", "StnID"]).reset_index(drop=True)
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--wells_csv", type=Path, default=DEFAULT_WELLS_CSV)
    ap.add_argument("--shapefile", type=Path, default=DEFAULT_SHAPEFILE)
    ap.add_argument("--cohort_csv", type=Path, default=DEFAULT_COHORT_CSV)
    ap.add_argument("--out_csv", type=Path, default=DEFAULT_OUT_CSV)
    args = ap.parse_args(argv)

    df = join_wells_to_cohort(args.wells_csv, args.shapefile, args.cohort_csv)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    print(f"\n[wells] wrote {args.out_csv} ({args.out_csv.stat().st_size / 1e6:.2f} MB)")

    # Per-aquifer summary
    counts = df.groupby("Study_area").size().sort_values(ascending=False)
    print(f"\n[wells] well counts per aquifer (top 5 + bottom 5 of "
          f"{len(counts)} aquifers with >=1 well):")
    print(counts.head(5).to_string())
    print("  ...")
    print(counts.tail(5).to_string())
    print(f"\n[wells] aquifers with >=10 wells:  "
          f"{(counts >= 10).sum()}")
    print(f"[wells] aquifers with >=50 wells:  "
          f"{(counts >= 50).sum()}")
    print(f"[wells] aquifers with >=100 wells: "
          f"{(counts >= 100).sum()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
