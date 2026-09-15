"""SH-domain corrections: TN-14 C20 / C30 replacement."""
from __future__ import annotations

from typing import Any

import numpy as np

from .io_aux import TN14
from .io_shm import Month


def replace_low_degrees(
    harm: Any, tn14: TN14, month: Month, c20: str, c30: str
) -> Any:
    """Apply TN-14 SLR replacement of C20 and/or C30 in spherical-harmonic domain.

    `c20` / `c30` are each either ``"replaced"`` or ``"original"``. Returns a
    fresh harmonics object so the caller can reuse the unmodified input
    across other variant combinations.
    """
    h = harm.copy()
    if c20 == "original" and c30 == "original":
        return h

    idx = tn14.at_mjd_start(month.mjd_start)
    if c20 == "replaced":
        h.clm[2, 0] = tn14.c20[idx]
    if c30 == "replaced":
        c30_val = tn14.c30[idx]
        if not np.isnan(c30_val):
            h.clm[3, 0] = c30_val
        # Pre-MJD 55987 TN-14 has no C30; in that case we silently retain the
        # GSM value. The 2020 POC window is fully post-MJD 55987 so this branch
        # is never hit here, but kept for forward compatibility.
    return h
