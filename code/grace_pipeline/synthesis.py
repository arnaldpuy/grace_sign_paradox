"""Spherical-harmonic synthesis at AOI grid cells, in cm water-equivalent."""
from __future__ import annotations

import pathlib
from functools import lru_cache
from typing import Any

import gravity_toolkit as gt
import numpy as np

from . import config as cfg


_GT_DATA = pathlib.Path(gt.__file__).parent / "data" / "love_numbers"


@lru_cache(maxsize=8)
def _cmwe_factor(lmax: int) -> np.ndarray:
    """Per-degree Stokes → cm water-equivalent conversion (Wahr et al. 1998).

    Cached because the factor depends only on lmax and PREM Love numbers,
    both of which are fixed across the entire variant matrix.
    """
    hl, kl, ll = gt.read_love_numbers(str(_GT_DATA), LMAX=lmax, FORMAT="tuple")
    units = gt.units(lmax=lmax)
    # Forward Wahr (1998) SH → spatial factor: ρ_e · a · (2l+1) / (3·(1+k_l)).
    # `units.spatial` gives the inverse (spatial → SH); `units.harmonic` is
    # what synthesis needs.
    units.harmonic(hl, kl, ll)
    return units.cmwe.copy()


def synthesise_to_aoi(harm: Any, aoi: cfg.AOI, lmax: int) -> np.ndarray:
    """Spherical harmonics → EWH (cm) at AOI grid cells, shape (nlat, nlon)."""
    h = harm.copy().convolve(_cmwe_factor(lmax))
    field_lon_lat = gt.harmonic_summation(
        h.clm, h.slm, aoi.lon, aoi.lat, LMIN=0, LMAX=lmax,
    )
    # gt.harmonic_summation returns (nlon, nlat); transpose to (nlat, nlon)
    # to match the standard convention used in the output Dataset.
    return np.asarray(field_lon_lat).T
