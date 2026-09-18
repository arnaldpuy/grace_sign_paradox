"""Download GRACE / GRACE-FO Level-2 Stokes-coefficient monthly files.

Fetches the monthly SHM record per centre and truncation
(CSR, JPL, GFZ × BA01, BB01).

This script downloads:
  * GSM-2 files: monthly spherical harmonic coefficients (the GRACE
    "gravity Stokes monthly" product). Two truncations: BA01 (60×60)
    and BB01 (96×96).
  * GAC files: atmospheric+oceanic dealiasing monthly products,
    required by some preprocessing chains for restoration. Match
    each GSM month.

Total volume estimate
---------------------
~120 months × 3 centres × ~5 files per month × ~3 MB each ≈ ~5 GB.
Comfortable on the current disk.

Source
------
**NASA PO.DAAC** is the unified distribution point for all three
processing centres' Level-2 products.

  Provider:        NASA PO.DAAC + Earthdata Login
  Search:          NASA CMR (Common Metadata Repository)
  Auth:            Earthdata Login via ~/.netrc (already configured)
  Collection IDs:
    CSR L2 RL06.3:     GRACEFO_L2_GRAV_CSR_RL06.3
    CSR L2 RL06 (pre-FO): GRACE_L2_GRAV_CSR_RL06
    JPL L2 RL06.3:     GRACEFO_L2_GRAV_JPL_RL06.3
    JPL L2 RL06 (pre-FO): GRACE_L2_GRAV_JPL_RL06
    GFZ L2 RL06.3:     GRACEFO_L2_GRAV_GFZ_RL06.3
    GFZ L2 RL06 (pre-FO): GRACE_L2_GRAV_GFZ_RL06

Time coverage notes
-------------------
  * GRACE-1 mission: 2002-04 → 2017-10
  * Mission gap:     2017-10 → 2018-06 (no data, by design)
  * GRACE-FO:        2018-06 → present
  * Our paper window 2015-04 → 2025-12 spans both missions; the
    gap is unavoidable. Existing `enumerate_months` handles missing
    months without complaint.

Filename conventions (matching what grace_pipeline expects)
-----------------------------------------------------------
The existing pipeline's `io_shm.py` parses filenames of the form:
  GSM-2_<yyyyddd>-<yyyyddd>_<MISSION>_<CENTRE>_<TRUNC>_<VER>
where MISSION ∈ {GRAC, GRFO} and CENTRE ∈ {UTCSR, JPLEM, GFZOP}.
This script preserves that naming exactly.

Alternative sources if PO.DAAC is unavailable
---------------------------------------------
  * CSR direct:   https://download.csr.utexas.edu/pub/slr/grace/
  * GFZ ISDC:     https://isdc.gfz-potsdam.de/
  * ICGEM:        https://icgem.gfz-potsdam.de/
"""
from __future__ import annotations

import netrc
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import requests

from . import config as cfg
from ._resumable_download import (
    DownloadRegistry,
    atomic_get,
    downloader_session,
    shutdown_requested,
)


# --------------------------------------------------------------------------
# Constants — collection IDs and CMR endpoint
# --------------------------------------------------------------------------

CMR_BASE = "https://cmr.earthdata.nasa.gov/search"
EDL_BASE = "https://urs.earthdata.nasa.gov"

# (centre_key, mission, collection_short_name), CMR-confirmed 2026-05-12. GRACE-FO
# bundles all products per centre; filtered by granule title prefix at parse time.
COLLECTIONS: tuple[tuple[str, str, str], ...] = (
    # GRACE-FO (2018-06 → present) — one collection per centre, multiple products
    ("CSR", "GRFO", "GRACEFO_L2_CSR_MONTHLY_0063"),
    ("JPL", "GRFO", "GRACEFO_L2_JPL_MONTHLY_0063"),
    ("GFZ", "GRFO", "GRACEFO_L2_GFZ_MONTHLY_0063"),
    # GRACE original (2002-04 → 2017-10) — per-product collections
    ("CSR", "GRAC", "GRACE_GSM_L2_GRAV_CSR_RL06"),
    ("CSR", "GRAC", "GRACE_GAC_L2_GRAV_CSR_RL06"),
    ("JPL", "GRAC", "GRACE_GSM_L2_GRAV_JPL_RL06"),
    ("JPL", "GRAC", "GRACE_GAC_L2_GRAV_JPL_RL06"),
    ("GFZ", "GRAC", "GRACE_GSM_L2_GRAV_GFZ_RL06"),
    ("GFZ", "GRAC", "GRACE_GAC_L2_GRAV_GFZ_RL06"),
)

# Truncations we need: matches cfg.TRUNCATIONS
TRUNCATIONS = (60, 96)
TRUNC_MNEM = {60: "BA01", 96: "BB01"}

# File types to download per month
FILE_TYPES = ("GSM", "GAC")

N_CONCURRENT = 4
N_RETRIES = 5
POLL_DELAY_S = 0.5


# --------------------------------------------------------------------------
# Auth (Earthdata Login)
# --------------------------------------------------------------------------

def _edl_credentials() -> tuple[str, str]:
    """Read NASA Earthdata Login from ~/.netrc."""
    nrc = netrc.netrc()
    auth = nrc.authenticators("urs.earthdata.nasa.gov")
    if auth is None:
        raise RuntimeError(
            "no NASA Earthdata credentials in ~/.netrc; expected "
            "machine urs.earthdata.nasa.gov entry"
        )
    user, _, password = auth
    return user, password


def _make_session() -> requests.Session:
    """Build a requests Session with EDL basic auth + redirect handling."""
    user, password = _edl_credentials()
    s = requests.Session()
    s.auth = (user, password)
    s.headers.update({"User-Agent": "grace-sign-paradox/0.1 grace-shm-downloader"})
    return s


# --------------------------------------------------------------------------
# CMR search
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Granule:
    """One month's worth of SHM file metadata."""
    centre: str
    mission: str
    file_type: str       # "GSM" or "GAC"
    truncation: int      # 60 or 96 — only meaningful for GSM
    yyyyddd_start: str
    yyyyddd_end: str
    granule_id: str
    download_url: str


def _cmr_search(
    session: requests.Session,
    collection_id: str,
    time_start: str,
    time_end: str,
) -> list[dict]:
    """Query CMR for granules in collection in time range. Returns raw items.

    CMR search is anonymous. We deliberately use a fresh
    `requests.get` rather than the EDL-auth-bearing session because
    CMR rejects requests that include basic-auth headers (401
    Unauthorized). The `session` parameter is retained for API
    consistency but unused for the search itself.
    """
    url = f"{CMR_BASE}/granules.json"
    out: list[dict] = []
    page = 1
    while True:
        params = {
            "short_name": collection_id,
            "temporal": f"{time_start}T00:00:00Z,{time_end}T23:59:59Z",
            "page_size": 2000,
            "page_num": page,
        }
        r = requests.get(url, params=params, timeout=120)
        r.raise_for_status()
        entries = r.json().get("feed", {}).get("entry", [])
        if not entries:
            break
        out.extend(entries)
        if len(entries) < 2000:
            break
        page += 1
    return out


def _parse_granule(
    centre: str, mission: str, item: dict,
) -> Granule | None:
    """Parse one CMR entry into a Granule. Returns None if the filename
    doesn't match an expected (GSM/GAC)-2 product or no download URL is
    discoverable.

    Title format (PODAAC convention, all RL06+ releases):
      <PRODUCT>-2_<yyyyddd>-<yyyyddd>_<MISSION>_<CENTRE>_<TRUNC>_<VER>
    e.g. GSM-2_2020001-2020031_GRFO_UTCSR_BA01_0603
    """
    title = item.get("title", "")
    parts = title.split("_")
    if len(parts) < 6:
        return None
    file_type_raw, date_range = parts[0], parts[1]
    file_type = file_type_raw.split("-")[0]  # GSM-2 → GSM
    if file_type not in FILE_TYPES:
        return None
    truncation_mnem = parts[4] if len(parts) > 4 else ""
    truncation = next(
        (t for t, m in TRUNC_MNEM.items() if m == truncation_mnem),
        0,
    )
    # For GAC products, truncation token is BC01 (not in TRUNC_MNEM).
    # We accept truncation=0 for non-GSM products since the pipeline
    # routes GAC by file_type, not by truncation.
    download_url = ""
    for link in item.get("links", []):
        rel = link.get("rel", "")
        href = link.get("href", "")
        if rel.endswith("/data#") and "archive.podaac.earthdata.nasa.gov/podaac-ops" in href:
            download_url = href
            break
    if not download_url:
        return None
    try:
        yyyyddd_start, yyyyddd_end = date_range.split("-")
    except ValueError:
        return None
    return Granule(
        centre=centre,
        mission=mission,
        file_type=file_type,
        truncation=truncation,
        yyyyddd_start=yyyyddd_start,
        yyyyddd_end=yyyyddd_end,
        granule_id=title,
        download_url=download_url,
    )


# --------------------------------------------------------------------------
# Per-granule download
# --------------------------------------------------------------------------

def _local_path(g: Granule) -> Path:
    """Where the granule lives on disk. Matches existing pipeline convention."""
    centre_dir = cfg.CENTRE_DIRS[g.centre]
    return centre_dir / g.granule_id


def _download_granule(
    session: requests.Session, g: Granule,
    registry: DownloadRegistry,
) -> Path | None:
    out_path = _local_path(g)
    if registry.is_complete(g.granule_id):
        return out_path  # registry shortcut — no stat needed
    result = atomic_get(
        g.download_url,
        out_path,
        session=session,
        progress_label=g.granule_id,
        # GRACE PO.DAAC delivers Content-Length; size-validated rename
        # is automatic.
    )
    if result is None:
        return None
    path, n_bytes = result
    registry.mark_complete(g.granule_id, n_bytes, extra={
        "centre": g.centre, "mission": g.mission,
        "file_type": g.file_type, "truncation": g.truncation,
    })
    time.sleep(POLL_DELAY_S)
    return path


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------

def download_all(
    time_start: str = "2015-04-01",
    time_end: str = "2025-12-31",
) -> list[Path]:
    print("\n=== GRACE / GRACE-FO L2 SHM download ===", flush=True)
    print(f"  window: {time_start} → {time_end}", flush=True)

    session = _make_session()
    print(f"  auth: Earthdata Login OK", flush=True)

    # Phase 1: enumerate granules across all (centre, mission, collection)
    granules: list[Granule] = []
    for centre, mission, collection_id in COLLECTIONS:
        print(f"  searching {centre} / {mission} ({collection_id})…",
              flush=True)
        items = _cmr_search(session, collection_id, time_start, time_end)
        for item in items:
            g = _parse_granule(centre, mission, item)
            if g is not None:
                granules.append(g)
        print(f"    → {len(items)} CMR entries", flush=True)

    # Filter to truncations we care about (GSM only — GAC has no truncation)
    keep: list[Granule] = []
    for g in granules:
        if g.file_type == "GAC":
            keep.append(g)
        elif g.file_type == "GSM" and g.truncation in TRUNCATIONS:
            keep.append(g)
    print(f"  total granules to download: {len(keep)}", flush=True)

    # Phase 2: download concurrently, skip-if-exists.
    # The downloader_session ctx installs SIGINT handler + opens
    # the persistent registry under data/raw/grace/.
    paths: list[Path] = []
    registry_dir = cfg.RAW_GRACE
    registry_dir.mkdir(parents=True, exist_ok=True)
    with downloader_session(registry_dir) as registry:
        with ThreadPoolExecutor(max_workers=N_CONCURRENT) as ex:
            futures = {
                ex.submit(_download_granule, session, g, registry):
                    g.granule_id
                for g in keep
            }
            for i, fut in enumerate(as_completed(futures)):
                if shutdown_requested():
                    # Drain quickly so ^C exits in seconds, not minutes
                    for f in futures:
                        f.cancel()
                    break
                p = None
                try:
                    p = fut.result()
                except Exception as exc:  # noqa: BLE001
                    gid = futures[fut]
                    registry.mark_failed(gid, str(exc))
                    print(f"  ! {gid}: "
                          f"{type(exc).__name__}: {str(exc)[:120]}",
                          flush=True)
                if p is not None:
                    paths.append(p)
                if i % 50 == 0 and i > 0:
                    print(f"  progress: {i}/{len(keep)} granules",
                          flush=True)

    return paths


# --------------------------------------------------------------------------
# Data-coverage notes
# --------------------------------------------------------------------------
# Notes: BB01 (96x96) only became routine with GRACE-FO, so truncation=96 is sparse
# for 2015-2017 (unrun slots stay NaN). Missing months at the GRACE/GRACE-FO
# interface (~2017-10..2018-06) are the expected mission gap, not download failures.


if __name__ == "__main__":
    t0 = datetime.now()
    paths = download_all()
    elapsed = (datetime.now() - t0).total_seconds() / 60
    total_mb = sum(p.stat().st_size for p in paths) / 1024**2
    print(f"\n=== summary ({elapsed:.1f} min) ===")
    print(f"  {len(paths)} files, {total_mb:.1f} MB total")
