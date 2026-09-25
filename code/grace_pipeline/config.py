"""Static configuration for the GRACE gridding pipeline.

Encodes the variant matrix (centre × truncation × filter × GIA × C20 × C30)
and the POC AOI / time window. Everything that downstream modules read as
"the project's preprocessing space" lives here.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import os
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
# REPO_ROOT is this repository. DATA_ROOT holds the raw and processed GRACE
# tree (data/raw/..., data/processed/...); set GRACE_DATA_ROOT to relocate it.
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get("GRACE_DATA_ROOT", REPO_ROOT / "data"))
RAW_GRACE = DATA_ROOT / "raw" / "grace"
AUX_DIR = RAW_GRACE / "auxiliary"
OUT_DIR = DATA_ROOT / "processed" / "grace"

CENTRE_DIRS: dict[str, Path] = {
    "CSR": RAW_GRACE / "csr",
    "JPL": RAW_GRACE / "jpl",
    "GFZ": RAW_GRACE / "gfz",
}
CENTRE_TOKEN: dict[str, str] = {"CSR": "UTCSR", "JPL": "JPLEM", "GFZ": "GFZOP"}

# GLDAS-2.1 Noah land-surface model — the non-groundwater storage terms
# (soil moisture, snow water equivalent, canopy) used to decompose GRACE
# TWS into GRACE-GWS. See grace_pipeline/gws_correction.py.
RAW_GLDAS = DATA_ROOT / "raw" / "gldas"
# California reservoir storage (CDEC) — the surface-water term of the GWS
# decomposition. Only material / tractable for the Central Valley; the IGP
# GWS correction is soil + snow + canopy only (surface water a stated
# limitation). See SURFACE_WATER below + download_cdec_reservoirs.py.
RAW_CDEC = DATA_ROOT / "raw" / "cdec"

TN14_FILE = AUX_DIR / "TN-14_C30_C20_GSFC_SLR.txt"
GIA_FILES: dict[str, Path] = {
    "Geruo_ICE5G_VM2": AUX_DIR / "Tellus_GIA_L3_A-WAHR_ICE5G-VM2_0.5-DEG_v1.0.nc",
    "Caron_2018":      AUX_DIR / "Tellus_GIA_L3_CARON-2018_0.5-DEG_v1.0.nc",
    "Peltier_ICE6G_D": AUX_DIR / "Tellus_GIA_L3_PELTIER_ICE6G-D_0.5-DEG_v1.0.nc",
}
# DDK kernels are not yet on disk — see register entry to be added when
# download_grace_aux.py is extended.
DDK_KERNEL_DIR = AUX_DIR / "ddk"

# Geocenter coefficient files. GravIS (Dahle/Murböck, GFZ) is the primary
# product because it is open-access (no Earthdata Login), single-file,
# centre-agnostic, and covers the full GRACE+GRACE-FO mission window
# (2002-04 → 2026-02 in the v0004 distribution).
GEOCENTER_DIR = AUX_DIR / "geocenter"
GRAVIS_GEOCENTER_FILE = GEOCENTER_DIR / "GRAVIS-2B_GFZOP_GEOCENTER_0004.dat"

# --------------------------------------------------------------------------
# AOI registry — analysis regions at 0.25°
# --------------------------------------------------------------------------
# Each study aquifer is an AOI in AOI_REGISTRY; REGION_SUFFIX maps region key ->
# output-filename suffix (IGP keeps the bare filename for back-compatibility).
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class AOI:
    lat_min: float = 24.0
    lat_max: float = 32.0
    lon_min: float = 73.0
    lon_max: float = 88.0
    res: float = 0.25

    @property
    def lat(self) -> np.ndarray:
        return np.arange(self.lat_min + self.res / 2, self.lat_max, self.res)

    @property
    def lon(self) -> np.ndarray:
        """Cell-centre longitudes in the AOI's native convention
        (−180..180; negative for the Central Valley)."""
        return np.arange(self.lon_min + self.res / 2, self.lon_max, self.res)

    @property
    def lon360(self) -> np.ndarray:
        """Cell-centre longitudes wrapped to 0..360, in the SAME order as
        `lon`. Auxiliary grids (GIA TELLUS fields, mascon products) are
        stored on a 0..360 longitude axis, so any `.interp`/`.sel` against
        them must use this convention. For the IGP (73..88 E) this is a
        no-op; for the Central Valley (−122..−118 E) it maps to 237..242."""
        return self.lon % 360.0


# Indo-Gangetic Plain — the original POC aquifer.
AOI_IGP = AOI()

# Central Valley box includes the western Sierra front (lon to -118.5): snowpack
# leakage into the valley is the signal the recipes disagree about.
AOI_CV = AOI(lat_min=35.0, lat_max=40.5, lon_min=-122.5, lon_max=-118.5,
             res=0.25)

# High Plains box: central+southern depletion core; excludes the recharge-dominated
# Nebraska Sand Hills north of ~41N (hydrologic outlier).
AOI_HP = AOI(lat_min=34.0, lat_max=41.0, lon_min=-103.0, lon_max=-98.0,
             res=0.25)

AOI_REGISTRY: dict[str, AOI] = {"igp": AOI_IGP, "cv": AOI_CV, "hp": AOI_HP}
REGION_SUFFIX: dict[str, str] = {"igp": "", "cv": "_cv", "hp": "_hp"}
REGION_LABEL: dict[str, str] = {
    "igp": "IGP", "cv": "Central Valley", "hp": "High Plains"}

# Surface-water term: ON for the Central Valley (CDEC reservoirs material), OFF for
# the IGP (no clean basin-wide series; documented as a limitation).
SURFACE_WATER: dict[str, bool] = {"igp": False, "cv": True, "hp": False}


def aoi_for_region(region: str) -> AOI:
    """Resolve a region key ('igp' | 'cv' | 'hp') to its AOI dataclass."""
    try:
        return AOI_REGISTRY[region]
    except KeyError:
        raise ValueError(
            f"unknown region {region!r}; expected one of "
            f"{sorted(AOI_REGISTRY)}") from None

# --------------------------------------------------------------------------
# POC_MONTHS: 89 GRACE-FO months (2018-06..2025-12) common to all centres/truncations
# (2018-08/09 have no solutions) -- the legacy default build window.
# APPENDIX_MONTHS: full GRACE+GRACE-FO 2002-04..2025-12 minus the 2017-07..2018-05
# bridge and 2018-08/09 commissioning gaps = 272 months; the paper's headline
# (W-FULL) window. Months missing at a centre are handled as NaN by io_shm.py.
# --------------------------------------------------------------------------
_GRACE_GAP_MONTHS = {(2018, 8), (2018, 9)}

POC_MONTHS: tuple[tuple[int, int], ...] = tuple(
    (y, m)
    for y in range(2018, 2026)
    for m in range(1, 13)
    if (y, m) >= (2018, 6)
    and (y, m) <= (2025, 12)
    and (y, m) not in _GRACE_GAP_MONTHS
)

APPENDIX_MONTHS: tuple[tuple[int, int], ...] = tuple(
    (y, m)
    for y in range(2002, 2026)
    for m in range(1, 13)
    if (
        ((y, m) >= (2002, 4) and (y, m) <= (2017, 6))
        or ((y, m) >= (2018, 6) and (y, m) <= (2025, 12))
    )
    and (y, m) not in _GRACE_GAP_MONTHS
)


def months_for_window(window: str) -> tuple[tuple[int, int], ...]:
    """Resolve a window label to its month tuple."""
    if window == "poc":
        return POC_MONTHS
    if window == "appendix":
        return APPENDIX_MONTHS
    raise ValueError(f"unknown window {window!r}; expected 'poc' or 'appendix'")

# GIA reference epoch = mid-window: corrected and uncorrected fields are equal here
# and diverge symmetrically, encoding that GIA contaminates the trend, not the level.
GIA_REF_EPOCH = np.datetime64("2022-04-01")

# --------------------------------------------------------------------------
# Variant axes (ORDER PRESERVED — used as Dataset coord ordering)
# --------------------------------------------------------------------------
CENTRES: tuple[str, ...] = ("CSR", "JPL", "GFZ")

# Truncation: BA01 = 60×60, BB01 = 96×96 (verified against JPL BB01 header).
TRUNCATIONS: tuple[int, ...] = (60, 96)
TRUNC_MNEM: dict[int, str] = {60: "BA01", 96: "BB01"}

# Filter axis order indexes the cube's `filter` dim: original POC five, then the two
# Swenson+Gaussian combos, then DDK2/DDK8 bracketing the DDK3/5/7 triplet; G400
# fills the Gaussian-radius gap.
FILTERS: tuple[str, ...] = (
    "G300", "G400", "G500",
    "Swenson_G300", "Swenson_G500",
    "DDK2", "DDK3", "DDK5", "DDK7", "DDK8",
)
GAUSSIAN_RADII_KM: dict[str, float] = {
    "G300": 300.0, "G400": 400.0, "G500": 500.0,
}
DDK_KERNEL_IDS: dict[str, int] = {
    "DDK2": 2, "DDK3": 3, "DDK5": 5, "DDK7": 7, "DDK8": 8,
}
# Second-stage radius after Swenson destriping; kept separate from GAUSSIAN_RADII_KM
# so Swenson_G* is never dispatched as a plain Gaussian by name lookup.
SWENSON_GAUSSIAN_KM: dict[str, float] = {
    "Swenson_G300": 300.0, "Swenson_G500": 500.0,
}

GIA_MODELS: tuple[str, ...] = ("none",) + tuple(GIA_FILES.keys())
C20_TREATMENTS: tuple[str, ...] = ("replaced", "original")
C30_TREATMENTS: tuple[str, ...] = ("replaced", "original")
# Geocenter axis. "none" must be index 0 so the existing 960-recipe
# (filter-only) cubes archived in _tier1_filters_only/ byte-match the
# geocenter="none" slice of the new 1,920-recipe cubes.
GEOCENTER_TREATMENTS: tuple[str, ...] = ("none", "GravIS")

# --------------------------------------------------------------------------
# Spherical-harmonic synthesis parameters
# --------------------------------------------------------------------------
# Love numbers: gravity_toolkit ships PREM-based defaults; we use them.
# Surfaced as a config knob so the choice is explicit and citable.
LOVE_NUMBERS_FILE: Path | None = None  # None ⇒ built-in PREM

# --------------------------------------------------------------------------
# Variant tuples
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Variant:
    centre: str
    truncation: int
    filter: str
    gia: str
    c20: str
    c30: str
    geocenter: str


def variant_matrix(centre: str | None = None) -> list[Variant]:
    """All non-time variant combinations. If `centre` is given, restrict to it.

    Geocenter is the fastest-varying axis in the Cartesian product (matches
    the R loader's column-major flattening; see load_grace_recipes_fun.R).
    """
    centres = (centre,) if centre is not None else CENTRES
    return [
        Variant(c, t, f, g, c20, c30, gc)
        for c, t, f, g, c20, c30, gc in product(
            centres, TRUNCATIONS, FILTERS, GIA_MODELS, C20_TREATMENTS, C30_TREATMENTS,
            GEOCENTER_TREATMENTS,
        )
    ]


def variant_matrix_size_per_month() -> int:
    """Sanity check. Tier-1 with geocenter: 2·10·4·2·2·2 = 640 per centre per
    month (was 2·10·4·2·2 = 320 in the filter-only Tier-1 intermediate, and
    2·5·4·2·2 = 160 in the original 5-filter set).
    """
    return (
        len(TRUNCATIONS)
        * len(FILTERS)
        * len(GIA_MODELS)
        * len(C20_TREATMENTS)
        * len(C30_TREATMENTS)
        * len(GEOCENTER_TREATMENTS)
    )


def existing_table(path: Path) -> Path:
    """Return `path` if it exists, else its gzipped twin when that exists.

    The release ships the large per-recipe trend tables as `.csv.gz`; pandas
    and data.table read either form transparently.
    """
    if path.exists():
        return path
    gz = path.with_name(path.name + ".gz")
    return gz if gz.exists() else path
