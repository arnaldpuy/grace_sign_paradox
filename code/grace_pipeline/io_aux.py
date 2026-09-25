"""Read auxiliary GRACE corrections (TN-14 SLR replacement, TELLUS GIA)."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr

from . import config as cfg


# --------------------------------------------------------------------------
# TN-14: GSFC SLR C20 / C30 replacement series
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class TN14:
    mjd_start: np.ndarray   # MJD of beginning of solution data span (matches GSM)
    mjd_end: np.ndarray
    time_start: np.ndarray  # decimal year, beginning of solution
    time_end: np.ndarray
    c20: np.ndarray
    c20_sigma: np.ndarray   # already in absolute units (1e-10 scaling applied)
    c30: np.ndarray         # NaN before MJD 55987
    c30_sigma: np.ndarray

    def at_mjd_start(self, mjd: float, tol_days: float = 15.0) -> int:
        """Index of the row whose MJD_start matches `mjd` within `tol_days`.

        GRACE-FO L2 file starts align with TN-14 to ~1 day; GRACE-mission
        files routinely sit 3-5 days off because the centres' start-date
        conventions diverged in 2002-2017. 15 days = half a month is the
        safe upper bound (TN-14 has exactly one row per GRACE month, so a
        ±15-day window can never match two rows).
        """
        diff = np.abs(self.mjd_start - mjd)
        idx = int(np.argmin(diff))
        if diff[idx] > tol_days:
            raise KeyError(
                f"TN-14 has no row within {tol_days} d of MJD {mjd} "
                f"(closest: {self.mjd_start[idx]}, Δ={diff[idx]:.2f} d)"
            )
        return idx


def load_tn14(path: Path = cfg.TN14_FILE) -> TN14:
    """Parse the TN-14 GSFC SLR file.

    Native columns (from the file header):
      1 MJD start, 2 year start (frac), 3 C20, 4 ΔC20×1e-10, 5 σC20×1e-10,
      6 C30, 7 ΔC30×1e-10, 8 σC30×1e-10, 9 MJD end, 10 year end (frac).

    Replaces upstream gravity_toolkit.SLR.C20 / C30, which crash on this
    distribution due to a numpy-scalar-coercion regression.
    """
    rows: list[list[str]] = []
    in_data = False
    with path.open("r") as fh:
        for line in fh:
            if not in_data:
                if line.startswith("Product:"):
                    in_data = True
                continue
            parts = line.split()
            if len(parts) >= 10 and parts[0].replace(".", "").replace("-", "").isdigit():
                rows.append(parts[:10])
    if not rows:
        raise ValueError(f"No data rows parsed from {path}")

    arr = np.array(
        [[float("nan") if p == "NaN" else float(p) for p in r] for r in rows]
    )
    return TN14(
        mjd_start=arr[:, 0],
        time_start=arr[:, 1],
        c20=arr[:, 2],
        c20_sigma=arr[:, 4] * 1e-10,
        c30=arr[:, 5],
        c30_sigma=arr[:, 7] * 1e-10,
        mjd_end=arr[:, 8],
        time_end=arr[:, 9],
    )


# --------------------------------------------------------------------------
# TELLUS L3 GIA mass-rate fields
# --------------------------------------------------------------------------
def load_gia_rate(name: str) -> xr.DataArray:
    """Load one TELLUS L3 GIA mass-rate field as an xarray.DataArray.

    Returns the native 0.5° global grid (lon 0–360°, lat -90 to 90) with
    units of m H₂O / yr. Filter-compatibility caveat: these fields are
    pre-filtered for 3° JPL Mascon compatibility, not for Gaussian/DDK SH
    solutions. They are retained because subtracting them from SH-derived
    storage is common published practice, and the smoothness of a GIA trend
    field keeps the mismatch residual small against the spread between GIA
    models (see the dataset comment in pipeline.py).
    """
    if name == "none":
        raise ValueError("load_gia_rate('none') is undefined; handle 'none' upstream")
    if name not in cfg.GIA_FILES:
        raise KeyError(f"Unknown GIA model: {name!r}")
    ds = xr.open_dataset(cfg.GIA_FILES[name])
    da = ds["GIA_mass_rate_3degJPL_MSCN"]
    da = da.rename({"lat": "lat", "lon": "lon"})
    da.attrs["gia_model"] = name
    da.attrs["source_file"] = cfg.GIA_FILES[name].name
    da.attrs["filter_compatibility"] = "JPL 3° Mascon (NOT SH-filter compatible)"
    return da
