"""Swenson & Wahr 2006 destriping of GRACE Stokes coefficients.

Removes the "correlated errors" that show up as north-south stripes in
GRACE mascons. For each order m ≥ 5, the per-order coefficient series is
split into even-l and odd-l sub-series, a quadratic is fitted to a
sliding window of adjacent degrees (window size narrows with m, following
Swenson & Wahr's exp(−m/10)·15 schedule), and the fit is subtracted from
the central coefficient. Orders m < 5 are passed through unchanged.

This module is a thin wrapper around `gravity_toolkit.destripe_harmonics`
(Sutterley 2023, MIT-licensed), which is itself a port of the original
Velicogna Fortran code. We re-export it through a project-typed function
so the destriping operator lives next to the other SH-domain filters and
so the call site is grep-able as ``swenson_destripe``.

References
----------
Swenson, S., and J. Wahr (2006), Post-processing removal of correlated
errors in GRACE data, Geophys. Res. Lett., 33, L08402,
doi:10.1029/2005GL025285.
"""
from __future__ import annotations

from typing import Any

from gravity_toolkit.destripe_harmonics import destripe_harmonics


def swenson_destripe(harm: Any, lmax: int, mmax: int | None = None) -> Any:
    """Destripe a gravity_toolkit.harmonics object (Swenson & Wahr 2006).

    Parameters
    ----------
    harm : gravity_toolkit.harmonics
        Input field. Not mutated — a fresh copy is returned.
    lmax : int
        Maximum spherical-harmonic degree to operate on. Must equal the
        field's own ``harm.lmax`` (caller's responsibility to keep them
        in sync; we don't reach into the harm object's internals here).
    mmax : int or None
        Maximum order. Defaults to ``lmax`` (square spectrum), which is
        the only case GRACE GSM files exercise.

    Returns
    -------
    gravity_toolkit.harmonics
        Copy of ``harm`` with ``clm`` and ``slm`` destriped for m ≥ 5.

    Notes
    -----
    The downstream Gaussian smoothing in ``filters_sh.apply_filter`` should
    be applied AFTER this destriping step — Swenson destriping reshapes
    the per-order spectrum but does not down-weight high degrees; the
    Gaussian still does that.
    """
    h = harm.copy()
    result = destripe_harmonics(h.clm, h.slm, LMAX=lmax, MMAX=mmax)
    h.clm = result["clm"]
    h.slm = result["slm"]
    return h
