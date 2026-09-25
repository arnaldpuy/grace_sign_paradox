"""Tier-A polygon masks: which global-grid cells fall inside each cohort aquifer.

For each of the 73 Jasechko aquifers (area >= 50,000 km^2, frozen in
`datasets/jasechko_2024/cohort_50K_n73.csv`), this module computes the
set of (lat_idx, lon_idx) cells of the canonical Tier-A 0.5 deg global
grid whose CENTROID falls inside the polygon.

Output artifact (read from both R and Python downstream):
    datasets/jasechko_2024/cell_masks_0p5deg.parquet

Schema:
    aquifer_idx   int    0..72, index into the cohort table
    Study_area    str    Jasechko name (matches cohort_50K_n73.csv)
    lat_idx       int    0..339, index into the global lat axis
    lon_idx       int    0..719, index into the global lon axis
    lat           float  cell centre latitude (deg)
    lon           float  cell centre longitude (deg, -180..180)

Usage
-----
    cd <repository>/code
    python -m grace_pipeline.polygon_masks --grid_res 0.5

Validation: a centroid-based mask is the simplest, defensible choice
(used by every gridded reanalysis -> aquifer-mean pipeline). For v2 we
could move to area-weighted intersection if a small-cell-fraction
aquifer turns out to need it.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

from . import config as cfg


# Default paths: cohort CSV and shapefile live under datasets/ in this repository.
_CODE_REPO = cfg.REPO_ROOT
_JAS_DIR = _CODE_REPO / "datasets" / "jasechko_2024"

DEFAULT_SHAPEFILE = _JAS_DIR / "aquifers_shp" / "jasechko_et_al_2024_aquifers.shp"
DEFAULT_COHORT_CSV = _JAS_DIR / "cohort_50K_n73.csv"
DEFAULT_OUT_PARQUET = _JAS_DIR / "cell_masks_0p5deg.parquet"


def build_global_grid(res: float = 0.5,
                       lat_min: float = -85.0, lat_max: float = 85.0,
                       lon_min: float = -180.0, lon_max: float = 180.0
                       ) -> tuple[np.ndarray, np.ndarray]:
    """Canonical Tier-A global grid (centre-aligned cells).

    Returns (lats, lons) 1-D arrays of cell-centre coordinates,
    matching the convention used in smoke_test_global.make_global_aoi().
    """
    lats = np.arange(lat_min + res / 2, lat_max, res)
    lons = np.arange(lon_min + res / 2, lon_max, res)
    return lats, lons


def build_combo_key(aquifer: str, broader: str | float) -> str:
    """Reconstruct Jasechko Table-16 style "Aquifer (Broader)" key.

    Empty or "-" Broader yields just the Aquifer name (matching how
    Jasechko formats sub-basins-without-broader-system in Table 16).
    """
    if broader is None or (isinstance(broader, float) and np.isnan(broader)) \
            or str(broader).strip() in ("", "-"):
        return str(aquifer)
    return f"{aquifer} ({broader})"


def compute_cell_masks(shapefile: Path = DEFAULT_SHAPEFILE,
                        cohort_csv: Path = DEFAULT_COHORT_CSV,
                        grid_res: float = 0.5) -> pd.DataFrame:
    """For each cohort aquifer, return (lat_idx, lon_idx) of cells inside.

    Uses centroid-based assignment: a cell belongs to a polygon iff its
    centre point lies within the polygon geometry.
    """
    # 1. Load full 1,693 shapefile + cohort filter
    shp = gpd.read_file(shapefile)
    cohort = pd.read_csv(cohort_csv)
    cohort_names = set(cohort["Study_area"].str.lower())
    shp["combo_key"] = [
        build_combo_key(a, b) for a, b in zip(shp["Aquifer"], shp["Broader"])
    ]
    cohort_shp = shp[shp["combo_key"].str.lower().isin(cohort_names)].copy()
    cohort_shp = cohort_shp.merge(
        cohort[["Study_area"]].assign(_key=cohort["Study_area"].str.lower()),
        left_on=cohort_shp["combo_key"].str.lower(),
        right_on="_key", how="left")
    if len(cohort_shp) != len(cohort):
        print(f"[masks] WARNING: cohort has {len(cohort)} rows but only "
              f"{len(cohort_shp)} match the shapefile. Missing: "
              f"{set(cohort['Study_area']) - set(cohort_shp['Study_area'])}")
    # Order by cohort CSV order (= descending area; see scout_jasechko.R)
    cohort_shp = cohort_shp.sort_values(
        "Study_area", key=lambda s: s.map(
            {n: i for i, n in enumerate(cohort["Study_area"])})).reset_index(
        drop=True)
    cohort_shp["aquifer_idx"] = range(len(cohort_shp))

    # 2. Build global grid + a GeoDataFrame of centre points for fast STRtree
    lats, lons = build_global_grid(res=grid_res)
    lat_idx_arr, lon_idx_arr = np.meshgrid(
        np.arange(lats.size), np.arange(lons.size), indexing="ij")
    lat_grid, lon_grid = np.meshgrid(lats, lons, indexing="ij")
    centres = gpd.GeoDataFrame(
        {"lat_idx": lat_idx_arr.ravel(),
         "lon_idx": lon_idx_arr.ravel(),
         "lat": lat_grid.ravel(),
         "lon": lon_grid.ravel()},
        geometry=gpd.points_from_xy(lon_grid.ravel(), lat_grid.ravel()),
        crs="EPSG:4326")
    print(f"[masks] global grid: {lats.size}x{lons.size} = "
          f"{centres.shape[0]:,} centre points")

    # 3. Spatial join: which polygon (if any) contains each cell centre
    print(f"[masks] joining {len(cohort_shp)} polygons against grid "
          f"(this is the slow step; ~1 min)")
    cell_to_aq = gpd.sjoin(centres, cohort_shp[["aquifer_idx", "Study_area",
                                                  "combo_key", "geometry"]],
                            how="inner", predicate="within")
    print(f"[masks] {len(cell_to_aq):,} (cell, aquifer) pairs")

    out = cell_to_aq[["aquifer_idx", "Study_area", "lat_idx", "lon_idx",
                      "lat", "lon"]].sort_values(
        ["aquifer_idx", "lat_idx", "lon_idx"]).reset_index(drop=True)

    # 4. Sanity print: cell counts per aquifer
    counts = out.groupby("Study_area").size().sort_values(ascending=False)
    print(f"\n[masks] cell-count distribution across {len(counts)} aquifers:")
    print(f"  min  / 25% / med / 75% / max  =  "
          f"{counts.min()} / {int(counts.quantile(0.25))} / "
          f"{int(counts.median())} / {int(counts.quantile(0.75))} / "
          f"{counts.max()}")
    if counts.min() < 1:
        print(f"  WARNING: {(counts < 1).sum()} aquifers have 0 cells")
    if counts.min() < 5:
        n_thin = (counts < 5).sum()
        print(f"  WARNING: {n_thin} aquifers have <5 cells (sub-pixel risk)")
        print(counts[counts < 5].to_string())
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--shapefile", type=Path, default=DEFAULT_SHAPEFILE)
    ap.add_argument("--cohort_csv", type=Path, default=DEFAULT_COHORT_CSV)
    ap.add_argument("--grid_res", type=float, default=0.5)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT_PARQUET)
    args = ap.parse_args(argv)

    df = compute_cell_masks(args.shapefile, args.cohort_csv, args.grid_res)
    try:
        df.to_parquet(args.out, index=False)
        out_path = args.out
    except (ImportError, ValueError) as e:
        out_path = args.out.with_suffix(".csv")
        df.to_csv(out_path, index=False)
        print(f"[masks] parquet write failed ({e}); falling back to CSV")
    size_mb = out_path.stat().st_size / 1e6
    print(f"\n[masks] wrote {out_path} ({size_mb:.2f} MB, "
           f"{len(df):,} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
