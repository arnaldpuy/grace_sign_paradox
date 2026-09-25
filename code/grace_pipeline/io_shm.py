"""Discover and read GRACE/GRACE-FO L2 SHM files."""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import gravity_toolkit as gt

from . import config as cfg


@dataclass(frozen=True)
class Month:
    year: int
    month: int
    file_path: Path
    time_start: np.datetime64
    time_end: np.datetime64
    time_mid: np.datetime64

    @property
    def mjd_start(self) -> float:
        # MJD = (datetime - 1858-11-17) in days
        return float(
            (self.time_start - np.datetime64("1858-11-17")) / np.timedelta64(1, "D")
        )


# Accepts both GRACE (GRAC, RL06 "0600") and GRACE-FO (GRFO, RL06.3 "0603")
# monthly SH releases. The mission token is captured so downstream code can
# tell them apart if needed; in practice we treat them identically once
# loaded.
_GSM_RE = re.compile(
    r"^GSM-2_(?P<ys>\d{4})(?P<ds>\d{3})-(?P<ye>\d{4})(?P<de>\d{3})_"
    r"(?P<mission>GRAC|GRFO)_"
    r"(?P<token>[A-Z]+)_(?P<mnem>B[AB]01)_(?P<release>060\d)$"
)


def _doy_to_date(year: int, doy: int) -> np.datetime64:
    return np.datetime64(f"{year:04d}-01-01") + np.timedelta64(doy - 1, "D")


def enumerate_months(centre: str, truncation: int,
                      months: tuple[tuple[int, int], ...] | None = None
                      ) -> list[Month]:
    """Find every GSM file matching `centre` / `truncation` whose midpoint
    lands inside the requested month set (defaults to POC_MONTHS)."""
    centre_dir = cfg.CENTRE_DIRS[centre]
    token = cfg.CENTRE_TOKEN[centre]
    mnem = cfg.TRUNC_MNEM[truncation]
    poc = set(months if months is not None else cfg.POC_MONTHS)

    found: dict[tuple[int, int], Month] = {}
    for path in centre_dir.glob("GSM-2_*"):
        m = _GSM_RE.match(path.name)
        if m is None or m["token"] != token or m["mnem"] != mnem:
            continue
        ys, ds = int(m["ys"]), int(m["ds"])
        ye, de = int(m["ye"]), int(m["de"])
        t_start = _doy_to_date(ys, ds)
        t_end = _doy_to_date(ye, de)
        # GSM time_coverage uses an open-ended end day-of-year (matches the
        # next month's first day for monthly solutions); midpoint computed on
        # the open interval.
        span_days = (t_end - t_start) / np.timedelta64(1, "D")
        t_mid = t_start + np.timedelta64(int(span_days // 2), "D")
        ym = (int(t_mid.astype("datetime64[Y]").astype(int)) + 1970,
              int(t_mid.astype("datetime64[M]").astype(int) % 12) + 1)
        if ym not in poc:
            continue
        candidate = Month(
            year=ym[0], month=ym[1], file_path=path,
            time_start=t_start, time_end=t_end, time_mid=t_mid,
        )
        if ym in found:
            # 2011-10 and 2015-04 ship two short-arc solutions at every
            # centre (mid-month and end-month windows) because of
            # accelerometer anomalies; keep the one whose midpoint is
            # closest to the canonical month-15 timestamp.
            target = np.datetime64(f"{ym[0]:04d}-{ym[1]:02d}-15")
            d_existing = abs(found[ym].time_mid - target)
            d_new = abs(candidate.time_mid - target)
            if d_new < d_existing:
                found[ym] = candidate
        else:
            found[ym] = candidate

    missing = poc - set(found)
    if missing and len(found) == 0:
        raise FileNotFoundError(
            f"{centre} {mnem}: no GSM files found for any POC month"
        )
    if missing:
        # GRACE-mission months frequently lack a solution at one or more
        # centres (instrument anomalies, accelerometer dropouts in 2015-17,
        # truncation BB01 not routinely produced pre-2015). Log and proceed
        # with the available subset; downstream code handles missing recipe
        # slots as NaN in the cube.
        import warnings
        warnings.warn(
            f"{centre} {mnem}: {len(missing)} POC months without a GSM "
            f"file (e.g. {sorted(missing)[:3]}…). Proceeding with "
            f"{len(found)} available months.",
            stacklevel=2,
        )
    return [found[ym] for ym in sorted(found)]


def load_shm(file_path: Path, lmax: int) -> Any:
    """Read a SHM file into a gravity_toolkit.harmonics object."""
    return gt.harmonics(lmax=lmax).from_SHM(str(file_path))
