"""SH-domain filtering for the variant matrix.

Three filter families:
  * Gaussian (G300 / G400 / G500): Jekeli (1981) per-degree weights.
  * Swenson (Swenson_G300 / Swenson_G500): Swenson & Wahr 2006 destriping
    of orders m≥5 followed by the corresponding Gaussian smoother.
  * DDK (DDK2 / DDK3 / DDK5 / DDK7 / DDK8): Kusche 2009 decorrelation
    kernels, per-order block-diagonal mixing of degrees.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Any

import gravity_toolkit as gt
import numpy as np

from . import config as cfg
from .ddk import DDK_FILENAME, apply_ddk, read_ddk_kernel
from .destripe import swenson_destripe


def _gaussian_weights(hw_km: float, lmax: int) -> np.ndarray:
    """Per-degree Jekeli (1981) Gaussian weights, normalised so weight[0] = 1.

    gravity_toolkit returns the raw 1/(2π)-normalised weights; we divide by
    weight[0] so the filter acts as identity on the global mean and the rest
    is a proper smoothing kernel.
    """
    w = gt.gauss_weights(hw_km, lmax)
    return w / w[0]


# Cache holds five DDK kernels in steady state (DDK2/3/5/7/8 each ~9.4 MB);
# size 16 leaves headroom for any future addition without thrashing.
@lru_cache(maxsize=16)
def _load_ddk(level: str):
    """Cached DDK kernel — kernels are 9.4 MB binary; load once per session."""
    path = cfg.DDK_KERNEL_DIR / DDK_FILENAME[level]
    if not path.exists():
        raise FileNotFoundError(
            f"DDK kernel {level} expected at {path}. Run download_grace_aux.py "
            f"(DDK section) to fetch it from strawpants/GRACE-filter (MIT)."
        )
    return read_ddk_kernel(path)


def apply_filter(harm: Any, filt: str, lmax: int) -> Any:
    """Apply Gaussian, Swenson-destripe-then-Gaussian, or DDK in SH domain.

    Returns a fresh harmonics copy — caller can reuse the input across
    other filter slots.
    """
    h = harm.copy()
    if filt in cfg.GAUSSIAN_RADII_KM:
        w = _gaussian_weights(cfg.GAUSSIAN_RADII_KM[filt], lmax)
        return h.convolve(w)
    if filt in cfg.SWENSON_GAUSSIAN_KM:
        # Swenson & Wahr 2006 canonical pipeline: destripe orders m≥5 first,
        # then smooth with a Gaussian of 300 or 500 km half-width. Order
        # matters — Gaussian-then-destripe attenuates the high-m signal that
        # the destriping polynomial is meant to fit against.
        h_d = swenson_destripe(h, lmax=lmax)
        w = _gaussian_weights(cfg.SWENSON_GAUSSIAN_KM[filt], lmax)
        return h_d.convolve(w)
    if filt in cfg.DDK_KERNEL_IDS:
        return apply_ddk(h, _load_ddk(filt))
    raise ValueError(f"Unknown filter: {filt!r}")
