"""Degree-1 (geocenter) coefficient injection.

GRACE/GRACE-FO GSM Level-2 files have ``C10 = C11 = S11 = 0`` by convention
(degree-1 is unobservable from satellite tracking alone). The geocenter
treatment axis covers whether and how those three coefficients are restored
before SH-domain filtering and synthesis.

Two levels are wired:

  * ``none``    — pass-through; the GSM degree-1 zeros are kept. This slice
                  of the cube byte-matches the pre-geocenter (filter-only,
                  960-recipe) cubes archived in
                  ``data/processed/grace/_tier1_filters_only/``.
  * ``GravIS``  — Dahle & Murböck (GFZ) combined degree-1 series, computed
                  via the Swenson & Wahr 2008 method with self-attraction
                  and loading (SAL) correction. Single file, centre-agnostic,
                  covers April 2002 → February 2026 (i.e. full POC and full
                  appendix windows). Source:
                  https://isdc-data.gfz.de/grace/GravIS/GFZ/Level-2B/aux_data/
                  Read with ``gravity_toolkit.geocenter.from_gravis``.

Why GravIS and not TN-13: TN-13 is the JPL/PO.DAAC per-centre geocenter file
historically cited in the GRACE community, but the PO.DAAC distribution
requires NASA Earthdata authentication and ships per-centre. GravIS uses the
same Swenson 2008 + SAL methodology, is open access without auth, covers the
full mission window, and is single-file (centre-agnostic). For the
preprocessing variance question, the substantive comparison
is "geocenter on" vs "geocenter off" — the per-centre vs centre-agnostic
distinction is a second-order refinement that future tiers can address.

Operator order: ``apply_geocenter`` is called between ``replace_low_degrees``
(C20/C30 SLR replacement) and ``apply_filter`` — it only mutates the
``clm[1, 0]``, ``clm[1, 1]``, ``slm[1, 1]`` triplet, which is orthogonal to
C20/C30 replacement (degrees 2-3 only) but feeds into both filtering and
synthesis.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any

import numpy as np

from . import config as cfg
from .io_shm import Month


@dataclass(frozen=True)
class GeocenterSeries:
    """Time series of degree-1 Stokes coefficients.

    ``time`` is in decimal years (e.g. 2018.456 = 2018-06-15-ish). ``c10``,
    ``c11``, ``s11`` are dimensionless Stokes coefficients (same units as
    ``harm.clm``/``harm.slm``).
    """
    time: np.ndarray  # decimal year
    c10: np.ndarray
    c11: np.ndarray
    s11: np.ndarray
    source: str


def _month_to_decimal_year(month: Month) -> float:
    """Mid-month date → decimal year matching GravIS time convention."""
    t_mid = month.time_mid.astype("datetime64[D]")
    year = int(t_mid.astype("datetime64[Y]").astype(int)) + 1970
    year_start = np.datetime64(f"{year:04d}-01-01")
    year_end = np.datetime64(f"{year + 1:04d}-01-01")
    year_len = (year_end - year_start) / np.timedelta64(1, "D")
    day_of_year = (t_mid - year_start) / np.timedelta64(1, "D")
    return year + day_of_year / year_len


@lru_cache(maxsize=4)
def load_geocenter(name: str) -> GeocenterSeries:
    """Load a geocenter series by treatment name.

    ``"none"`` returns an all-zeros series (one entry; nearest-time lookup
    in ``apply_geocenter`` will always return zeros). ``"GravIS"`` reads the
    Dahle/Murböck file via ``gravity_toolkit.geocenter.from_gravis``.
    """
    if name == "none":
        return GeocenterSeries(
            time=np.array([0.0]),
            c10=np.array([0.0]),
            c11=np.array([0.0]),
            s11=np.array([0.0]),
            source="zeros (GSM convention)",
        )
    if name == "GravIS":
        # Lazy import — gravity_toolkit is already a runtime dep but keep
        # this module importable even if the package is broken.
        from gravity_toolkit.geocenter import geocenter
        path = cfg.GRAVIS_GEOCENTER_FILE
        if not path.exists():
            raise FileNotFoundError(
                f"GravIS geocenter file expected at {path}. Fetch via "
                f"download_grace_aux.py (geocenter section) from "
                f"isdc-data.gfz.de/grace/GravIS/GFZ/Level-2B/aux_data/"
            )
        gc = geocenter().from_gravis(str(path))
        return GeocenterSeries(
            time=np.asarray(gc.time, dtype=np.float64),
            c10=np.asarray(gc.C10, dtype=np.float64),
            c11=np.asarray(gc.C11, dtype=np.float64),
            s11=np.asarray(gc.S11, dtype=np.float64),
            source=path.name,
        )
    raise ValueError(
        f"Unknown geocenter treatment {name!r}; expected one of "
        f"{cfg.GEOCENTER_TREATMENTS}"
    )


def apply_geocenter(harm: Any, series: GeocenterSeries, month: Month,
                     tol_years: float = 0.05) -> Any:
    """Inject degree-1 coefficients into a copy of ``harm`` for one month.

    Picks the series entry whose decimal year is closest to the month's
    mid-time. Raises if the closest entry is further than ``tol_years``
    (~18 days) away — that protects against silently using a stale
    coefficient when the GravIS series has not yet been updated for a
    recent GRACE-FO month.

    ``"none"`` (zeros series of length 1) bypasses the distance check —
    it is always the correct value to substitute.
    """
    h = harm.copy()
    if series.source.startswith("zeros"):
        # GSM convention: C10 = C11 = S11 = 0. No-op (already zero in GSM).
        return h

    t_target = _month_to_decimal_year(month)
    idx = int(np.argmin(np.abs(series.time - t_target)))
    dt = float(np.abs(series.time[idx] - t_target))
    if dt > tol_years:
        raise ValueError(
            f"Geocenter series {series.source!r} has no entry within "
            f"{tol_years} yr of month {month.year}-{month.month:02d} "
            f"(decimal {t_target:.4f}); closest is {series.time[idx]:.4f} "
            f"({dt:.4f} yr away). Update the file or expand the tolerance."
        )

    h.clm[1, 0] = series.c10[idx]
    h.clm[1, 1] = series.c11[idx]
    h.slm[1, 1] = series.s11[idx]
    return h
