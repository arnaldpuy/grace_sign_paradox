"""Download the ERA5-Land monthly auxiliary storage terms from the CDS.

Second, independent auxiliary for the TWS -> GWS sensitivity analysis
(gws_aux_aquifer.py / the robustness-gws chunk): where GLDAS-2.1 Noah is
NASA's land-surface system, ERA5-Land is ECMWF's, so the difference
between the two per-aquifer auxiliary trends measures the variance the
auxiliary-model CHOICE adds on top of the 1,920-recipe preprocessing
spread - the quantity that makes the floor claim conservative.

Variables (monthly means, 2002-01 .. 2025-12, regridded server-side to
0.25 deg to match the GLDAS pathway):
  * volumetric_soil_water_layer_1..4  (m3/m3; layer thicknesses
    0.07 / 0.21 / 0.72 / 1.89 m -> column water in cm)
  * snow_depth_water_equivalent       (m of water equivalent)
  * skin_reservoir_content            (m of water equivalent; the
    canopy-interception analogue of GLDAS CanopInt)

Auth: ~/.cdsapirc. Two requests (12 years each) to stay under the CDS
field-count limit. Output: data/raw/era5_land/era5_land_aux_<years>.nc
"""
from __future__ import annotations

from pathlib import Path

import cdsapi

from . import config as cfg

OUT_DIR = cfg.DATA_ROOT / "raw" / "era5_land"
VARIABLES = ["volumetric_soil_water_layer_1", "volumetric_soil_water_layer_2",
             "volumetric_soil_water_layer_3", "volumetric_soil_water_layer_4",
             "snow_depth_water_equivalent", "skin_reservoir_content"]
BATCHES = [(2002, 2013), (2014, 2025)]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    client = cdsapi.Client()
    for y0, y1 in BATCHES:
        target = OUT_DIR / f"era5_land_aux_{y0}_{y1}.nc"
        if target.exists() and target.stat().st_size > 1e6:
            print(f"already present, skipping: {target.name}", flush=True)
            continue
        print(f"requesting {y0}-{y1} ...", flush=True)
        client.retrieve(
            "reanalysis-era5-land-monthly-means",
            {
                "product_type": ["monthly_averaged_reanalysis"],
                "variable": VARIABLES,
                "year": [str(y) for y in range(y0, y1 + 1)],
                "month": [f"{m:02d}" for m in range(1, 13)],
                "time": ["00:00"],
                "grid": [0.25, 0.25],
                "data_format": "netcdf",
                "download_format": "unarchived",
            },
            str(target))
        print(f"done: {target.name} "
              f"({target.stat().st_size / 1e6:.0f} MB)", flush=True)


if __name__ == "__main__":
    main()
