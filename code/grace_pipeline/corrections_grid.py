"""Grid-domain corrections: GIA rate × time-offset subtraction at AOI cells."""
from __future__ import annotations

import numpy as np
import xarray as xr

from . import config as cfg
from .io_shm import Month


_DAYS_PER_YEAR = 365.25
_M_TO_CM = 100.0


def subtract_gia_at_epoch(
    ewh_cm: np.ndarray,
    gia_rate_m_per_yr: xr.DataArray | None,
    month: Month,
    ref_epoch: np.datetime64,
    aoi: cfg.AOI,
) -> np.ndarray:
    """Subtract GIA cumulative correction at AOI grid cells.

    GIA correction = rate × (month_mid − ref_epoch). The correction is zero
    at ref_epoch and grows linearly outward, encoding GIA as a trend
    contamination (not a state contamination).

    Returns a fresh (nlat, nlon) array in cm water-equivalent.
    """
    if gia_rate_m_per_yr is None:
        return ewh_cm.copy()

    delta_t_days = (month.time_mid - ref_epoch) / np.timedelta64(1, "D")
    delta_t_years = float(delta_t_days) / _DAYS_PER_YEAR

    # Bilinear interpolation from the native 0.5° grid to the AOI 0.25° grid.
    # GIA fields are global lon 0–360°, lat -90 to 90. The AOI is queried in
    # the file's 0–360 convention via `aoi.lon360` (a no-op for the IGP at
    # 73–88E; maps the Central Valley −122..−118E to 237..242E). `lon360`
    # preserves cell order, so the returned array still aligns with `aoi.lon`.
    gia_aoi = gia_rate_m_per_yr.interp(
        lat=aoi.lat, lon=aoi.lon360, method="linear", assume_sorted=True
    )
    correction_cm = np.asarray(gia_aoi.values) * delta_t_years * _M_TO_CM
    return ewh_cm - correction_cm
