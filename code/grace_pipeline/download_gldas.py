"""Download GLDAS-2.1 Noah monthly land-surface fields from NASA GES DISC.

These are the non-groundwater storage terms used to decompose GRACE TWS
into GRACE-GWS (groundwater storage): soil moisture (four layers,
0-200 cm), snow water equivalent, and canopy-intercepted water. See
grace_pipeline/gws_correction.py for how they are combined.

Product: GLDAS_NOAH025_M v2.1 — 0.25° (matches the AOI grid), monthly,
global, 2000-01 → present. We pull full-globe granules and clip in the
builder (same as mascon.py); the monthly NetCDFs are small.

Auth: NASA Earthdata via the repo .netrc (`machine urs.earthdata.nasa.gov`).
GLDAS collections require the GES DISC application to be authorised once
in the Earthdata URS profile — a 401 from earthaccess.download means that,
not a .netrc problem.
"""
from __future__ import annotations

from pathlib import Path

import earthaccess

from . import config as cfg

GLDAS_SHORT_NAME = "GLDAS_NOAH025_M"
GLDAS_VERSION = "2.1"

# Download a generous span covering both analysis windows: the appendix
# window starts 2002-04 and the POC window ends 2025-12.
TEMPORAL_START = "2002-01-01"
TEMPORAL_END = "2025-12-31"


def download_gldas(out_dir: Path = cfg.RAW_GLDAS) -> list[Path]:
    """Fetch all GLDAS-2.1 Noah monthly granules for 2002-2025."""
    out_dir.mkdir(parents=True, exist_ok=True)
    print("=== GLDAS-2.1 Noah monthly ===")
    earthaccess.login()
    results = earthaccess.search_data(
        short_name=GLDAS_SHORT_NAME,
        version=GLDAS_VERSION,
        temporal=(TEMPORAL_START, TEMPORAL_END),
    )
    print(f"  found {len(results)} granules "
          f"({TEMPORAL_START} -> {TEMPORAL_END})")
    if not results:
        return []
    files = earthaccess.download(results, local_path=str(out_dir))
    paths = [Path(f) for f in files]
    total_mb = sum(p.stat().st_size for p in paths if p.exists()) / 1024**2
    print(f"  -> {len(paths)} files, {total_mb:.0f} MB in {out_dir}")
    return paths


if __name__ == "__main__":
    download_gldas()
