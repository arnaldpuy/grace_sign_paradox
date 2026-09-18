"""Climate-zone classification for the 73-cohort aquifers.

Computes the area-weighted mean Aridity Index (AI = MAP / MAE) for each
of the 73 Jasechko cohort polygons using the CGIAR-CSI Global-Aridity
v3 raster (Zomer et al. 2022). Classifies each aquifer per the bins used
in Jasechko et al. 2024:

    AI < 0.03         hyper-arid
    0.03 <= AI < 0.20 arid
    0.20 <= AI < 0.50 semi-arid
    0.50 <= AI < 0.65 dry sub-humid
    AI >= 0.65        humid

The CGIAR-CSI raster stores AI x 10000 as int16; we rescale to float
per the published convention.

Output: datasets/jasechko_2024/cohort_climate_zones.csv
    columns: Study_area, ai_mean, ai_std, n_pixels, climate_zone
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterstats import zonal_stats

from . import config as cfg


_CODE_REPO = cfg.REPO_ROOT
_JAS_DIR = _CODE_REPO / "datasets" / "jasechko_2024"

DEFAULT_AI_RASTER = (_CODE_REPO / "datasets" / "cgiar_csi_v3"
                      / "Global-AI_ET0__annual_v3_1" / "ai_v31_yr.tif")
DEFAULT_SHAPEFILE = _JAS_DIR / "aquifers_shp" / "jasechko_et_al_2024_aquifers.shp"
DEFAULT_COHORT_CSV = _JAS_DIR / "cohort_50K_n73.csv"
DEFAULT_OUT_CSV = _JAS_DIR / "cohort_climate_zones.csv"

# Jasechko bins (matching their Methods section)
CLIMATE_BINS = [
    (-np.inf, 0.03,  "hyper-arid"),
    (0.03,     0.20, "arid"),
    (0.20,     0.50, "semi-arid"),
    (0.50,     0.65, "dry sub-humid"),
    (0.65,     np.inf, "humid"),
]


def classify_ai(ai: float) -> str:
    for lo, hi, label in CLIMATE_BINS:
        if lo <= ai < hi:
            return label
    return "unknown"


def compute_climate_zones(raster_path: Path = DEFAULT_AI_RASTER,
                           shapefile: Path = DEFAULT_SHAPEFILE,
                           cohort_csv: Path = DEFAULT_COHORT_CSV,
                           out_csv: Path = DEFAULT_OUT_CSV) -> pd.DataFrame:
    # Filter shapefile to cohort
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
        right_on="_k", how="left").reset_index(drop=True)
    print(f"[climate] cohort polygons: {len(shp)}")

    # rasterstats zonal_stats; the AI raster is int16 with scale = 1/10000
    print(f"[climate] reading {raster_path.name} (~400 MB; this can take ~1 min)")
    stats = zonal_stats(
        shp, str(raster_path),
        stats=["mean", "std", "count", "nodata"],
        nodata=-32768,   # CGIAR-CSI nodata
        all_touched=False,
    )
    # AI raster is stored x10000; rescale
    SCALE = 1.0 / 10000.0
    rows: list[dict] = []
    for poly, s in zip(shp.itertuples(index=False), stats):
        ai_mean = (s["mean"] or 0) * SCALE if s["mean"] is not None else np.nan
        ai_std = (s["std"] or 0) * SCALE if s["std"] is not None else np.nan
        n_pix = int(s.get("count") or 0)
        rows.append({
            "Study_area": poly.Study_area,
            "ai_mean": float(ai_mean) if np.isfinite(ai_mean) else np.nan,
            "ai_std": float(ai_std) if np.isfinite(ai_std) else np.nan,
            "n_pixels": n_pix,
            "climate_zone": classify_ai(ai_mean) if np.isfinite(ai_mean)
                              else "unknown",
        })
    df = pd.DataFrame(rows)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"\n[climate] wrote {out_csv} "
          f"({out_csv.stat().st_size / 1e3:.1f} kB)")

    print("\n[climate] zone distribution:")
    print(df["climate_zone"].value_counts().to_string())
    print(f"\n[climate] AI summary:")
    print(df["ai_mean"].describe().to_string())
    return df


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--raster", type=Path, default=DEFAULT_AI_RASTER)
    ap.add_argument("--shapefile", type=Path, default=DEFAULT_SHAPEFILE)
    ap.add_argument("--cohort_csv", type=Path, default=DEFAULT_COHORT_CSV)
    ap.add_argument("--out_csv", type=Path, default=DEFAULT_OUT_CSV)
    args = ap.parse_args(argv)

    compute_climate_zones(args.raster, args.shapefile, args.cohort_csv,
                            args.out_csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
