"""Tier 1 cube-level byte-parity test against the pre-expansion archive.

Loads one of the new 10-filter cubes and the corresponding 5-filter cube in
``_pre_expansion/``, slices to the 5 shared filter names, and reports
max|Δewh_abs|, max|Δewh_anom| per filter. The expectation is exact zero
on every slice — the new cube's geocenter="none"-implicit slice on the
5 shared filters should byte-match the archived cube.

Usage
-----
    python -m code.grace_pipeline.verify_tier1_parity              # all 12 cubes
    python -m code.grace_pipeline.verify_tier1_parity --only csr/igp/poc
"""
from __future__ import annotations

import argparse
import sys
from itertools import product
from pathlib import Path

import numpy as np
import xarray as xr

from . import config as cfg

SHARED_FILTERS: tuple[str, ...] = ("G300", "G500", "DDK3", "DDK5", "DDK7")
PRE_DIR = cfg.OUT_DIR / "_pre_expansion"


def _cube_filename(centre: str, region: str, window: str) -> str:
    centre_l = centre.lower()
    region_suffix = cfg.REGION_SUFFIX[region]
    window_suffix = "_appendix" if window == "appendix" else ""
    return f"grace_ewh_variants_{centre_l}{region_suffix}{window_suffix}.nc"


def diff_one_cube(centre: str, region: str, window: str) -> dict[str, float]:
    """Return per-filter max-abs differences for one (centre, region, window).

    Reports NaN for any filter where one side is missing, or for a cube
    pair where the new cube has not yet been written.
    """
    fname = _cube_filename(centre, region, window)
    new_path = cfg.OUT_DIR / fname
    old_path = PRE_DIR / fname

    out: dict[str, float] = {}
    if not new_path.exists():
        out["__status__"] = float("nan")
        out["__error__"] = -1.0  # new cube missing
        return out
    if not old_path.exists():
        out["__status__"] = float("nan")
        out["__error__"] = -2.0  # archived cube missing
        return out

    ds_new = xr.open_dataset(new_path)
    ds_old = xr.open_dataset(old_path)
    try:
        max_overall_abs = 0.0
        max_overall_anom = 0.0
        for f in SHARED_FILTERS:
            if f not in ds_new.filter.values or f not in ds_old.filter.values:
                out[f] = float("nan")
                continue
            a_abs = ds_new.ewh_abs.sel(filter=f).values
            b_abs = ds_old.ewh_abs.sel(filter=f).values
            a_anom = ds_new.ewh_anom.sel(filter=f).values
            b_anom = ds_old.ewh_anom.sel(filter=f).values
            d_abs = float(np.nanmax(np.abs(a_abs - b_abs)))
            d_anom = float(np.nanmax(np.abs(a_anom - b_anom)))
            out[f + ":abs"] = d_abs
            out[f + ":anom"] = d_anom
            max_overall_abs = max(max_overall_abs, d_abs)
            max_overall_anom = max(max_overall_anom, d_anom)
        out["__max_abs__"] = max_overall_abs
        out["__max_anom__"] = max_overall_anom
    finally:
        ds_new.close()
        ds_old.close()
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--only", default=None,
        help="One '<centre>/<region>/<window>' triple to check (lowercased), "
             "e.g. 'csr/igp/poc'. Default: all 12 cubes.",
    )
    args = ap.parse_args(argv)

    cubes = []
    if args.only is not None:
        c, r, w = args.only.split("/")
        cubes.append((c.upper(), r, w))
    else:
        for c, r, w in product(cfg.CENTRES, ("igp", "cv", "hp"), ("poc",)):
            cubes.append((c, r, w))
        for c in cfg.CENTRES:
            cubes.append((c, "igp", "appendix"))

    overall_fail = False
    print(f"{'cube':<28} {'status':<10} {'max|Δabs|':>14} {'max|Δanom|':>14}")
    print("-" * 70)
    for c, r, w in cubes:
        label = f"{c.lower()}/{r}/{w}"
        d = diff_one_cube(c, r, w)
        if "__error__" in d:
            tag = "MISSING-NEW" if d["__error__"] == -1 else "MISSING-OLD"
            print(f"{label:<28} {tag:<10}")
            continue
        max_abs = d.get("__max_abs__", float("nan"))
        max_anom = d.get("__max_anom__", float("nan"))
        # Strict tolerance: byte-parity means literally zero. Allow a single
        # ULP of float64 slack just in case xarray's _FillValue → NaN round
        # trip introduces it (it shouldn't, but be defensive).
        passed = max_abs == 0.0 and max_anom == 0.0
        tag = "OK" if passed else "DIVERGE"
        if not passed:
            overall_fail = True
        print(f"{label:<28} {tag:<10} {max_abs:>14.3e} {max_anom:>14.3e}")
    return 1 if overall_fail else 0


if __name__ == "__main__":
    sys.exit(main())
