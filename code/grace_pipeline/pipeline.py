"""Variant-matrix walker for the GRACE EWH pipeline.

The order of operations is:
  load SHM → SH-domain C20/C30 replacement → SH-domain geocenter injection →
  SH-domain filter → AOI synthesis (cmwe) → grid-domain GIA subtraction.

Loop nesting is chosen so each expensive step is reused as many times
as possible inside the loop:
  load_shm                  once per (month, truncation)
  replace_low_degrees       once per (month, truncation, c20, c30)
  apply_geocenter           once per (month, truncation, c20, c30, geocenter)
  apply_filter              once per (month, truncation, c20, c30, geocenter, filter)
  synthesise_to_aoi         once per (month, truncation, c20, c30, geocenter, filter)
  subtract_gia_at_epoch     once per (month, ..., gia)

Geocenter slots between replace_low_degrees and apply_filter because
apply_geocenter only mutates clm[1,*]/slm[1,1] (orthogonal to C20/C30) but
the filter operators read all degrees.

`filters` is a runtime parameter so the matrix can be partially populated
when DDK kernels are missing — the cube keeps all filter slots; the
unrun ones stay NaN, and `output.write_centre` records which were realised.
"""
from __future__ import annotations

import argparse
import sys
import time as _time
from typing import Iterable

import numpy as np
import xarray as xr

from . import config as cfg
from .corrections_grid import subtract_gia_at_epoch
from .corrections_sh import replace_low_degrees
from .filters_sh import apply_filter
from .geocenter import apply_geocenter, load_geocenter
from .io_aux import load_gia_rate, load_tn14
from .io_shm import Month, enumerate_months, load_shm
from .output import write_centre
from .synthesis import synthesise_to_aoi


def compute_centre(
    centre: str,
    filters: Iterable[str] = cfg.FILTERS,
    window: str = "poc",
    region: str = "igp",
) -> tuple[xr.Dataset, list[Month]]:
    """Compute every variant for one centre. Returns (Dataset, months_used).

    Parameters
    ----------
    centre : "CSR" | "JPL" | "GFZ"
    filters : iterable of filter names (default = all four)
    window : "poc" (89 GRACE-FO months, headline) or "appendix"
        (272 months full GRACE + GRACE-FO, supplementary panel).
    region : "igp" (Indo-Gangetic Plain, default) or "cv" (California
        Central Valley). GRACE Stokes coefficients are global, so the
        region only changes the AOI the cube is synthesised onto.
    """
    aoi = cfg.aoi_for_region(region)
    nlat, nlon = aoi.lat.size, aoi.lon.size
    filters = tuple(filters)
    if not set(filters).issubset(cfg.FILTERS):
        raise ValueError(f"Unknown filter(s) in {filters!r}")

    window_months = cfg.months_for_window(window)
    n_t = len(window_months)
    n_tr = len(cfg.TRUNCATIONS)
    n_f = len(cfg.FILTERS)
    n_g = len(cfg.GIA_MODELS)
    n_c20 = len(cfg.C20_TREATMENTS)
    n_c30 = len(cfg.C30_TREATMENTS)
    n_gc = len(cfg.GEOCENTER_TREATMENTS)

    arr = np.full(
        (n_t, n_tr, n_f, n_g, n_c20, n_c30, n_gc, nlat, nlon),
        np.nan, dtype=np.float64,
    )
    # Global month → t-index inside this window. Truncations don't share
    # calendar coverage in the appendix window (BB01 is sparse pre-2015), so
    # every variant is positioned by its month's index in the window list.
    month_to_t = {ym: i for i, ym in enumerate(window_months)}
    times = np.array(
        [np.datetime64(f"{y:04d}-{m:02d}-15") for (y, m) in window_months],
        dtype="datetime64[D]",
    )

    tn14 = load_tn14()
    gia_rates: dict[str, xr.DataArray | None] = {"none": None}
    for name in cfg.GIA_MODELS:
        if name != "none":
            gia_rates[name] = load_gia_rate(name)
    geocenter_series = {name: load_geocenter(name)
                         for name in cfg.GEOCENTER_TREATMENTS}

    months_first: list[Month] = []
    t_start = _time.time()
    for tr_idx, trunc in enumerate(cfg.TRUNCATIONS):
        months = enumerate_months(centre, trunc, months=window_months)
        if tr_idx == 0:
            months_first = months

        for m in months:
            t_idx = month_to_t[(m.year, m.month)]
            times[t_idx] = m.time_mid.astype("datetime64[D]")

            harm_raw = load_shm(m.file_path, lmax=trunc)
            elapsed_min = (_time.time() - t_start) / 60
            print(
                f"[{centre}] {cfg.TRUNC_MNEM[trunc]} {m.year}-{m.month:02d} "
                f"(t={t_idx + 1}/{n_t}) elapsed {elapsed_min:.1f} min",
                flush=True,
            )

            for c20_idx, c20 in enumerate(cfg.C20_TREATMENTS):
                for c30_idx, c30 in enumerate(cfg.C30_TREATMENTS):
                    harm_corr = replace_low_degrees(harm_raw, tn14, m, c20, c30)

                    for gc_idx, gc_name in enumerate(cfg.GEOCENTER_TREATMENTS):
                        harm_gc = apply_geocenter(
                            harm_corr, geocenter_series[gc_name], m)

                        for f_idx, filt in enumerate(cfg.FILTERS):
                            if filt not in filters:
                                continue
                            harm_filt = apply_filter(harm_gc, filt, trunc)
                            ewh = synthesise_to_aoi(harm_filt, aoi, trunc)

                            for g_idx, gia in enumerate(cfg.GIA_MODELS):
                                arr[t_idx, tr_idx, f_idx, g_idx,
                                    c20_idx, c30_idx, gc_idx] = \
                                    subtract_gia_at_epoch(
                                        ewh, gia_rates[gia], m,
                                        cfg.GIA_REF_EPOCH, aoi,
                                    )

    return _assemble_dataset(arr, times, centre, region), months_first


def _assemble_dataset(arr: np.ndarray, times: np.ndarray, centre: str,
                      region: str = "igp") -> xr.Dataset:
    """Assemble the per-centre Dataset with both absolute and anomaly fields.

    We store
    both `ewh_abs` (full synthesised EWH, ~10⁹ cm dynamic range) and
    `ewh_anom` (= ewh_abs − time-mean over the reference period, ~±25 cm),
    so downstream consumers see the textbook view immediately without losing
    the absolute-provenance trail. Reference period is the full POC window.

    Note on the GIA axis: anomaly subtraction does NOT collapse the
    gia_model axis. It removes only the time-mean of the GIA correction (≈0 by
    construction); the time-varying component is preserved. At the IGP,
    over a 12-month window, GIA spread across models is ~0.06 cm — small
    compared to ~8 cm hydrology signal — but the axis is faithful in both
    `ewh_abs` and `ewh_anom`.
    """
    aoi = cfg.aoi_for_region(region)
    region_label = cfg.REGION_LABEL[region]
    dims = ("time", "truncation", "filter", "gia_model",
            "c20_treatment", "c30_treatment", "geocenter", "lat", "lon")

    # Anomaly: subtract time-mean (NaN-aware so unrun filter slots stay NaN).
    arr_mean = np.nanmean(arr, axis=0, keepdims=True)
    arr_anom = arr - arr_mean

    common_comment = (
        "GIA correction (when applied) uses the gridded TELLUS L3 mass-rate "
        "fields, which are pre-filtered for 3° JPL Mascon compatibility "
        "rather than for Gaussian/DDK SH solutions. Subtracting these grids "
        "from SH-derived storage is common published practice, which the "
        "recipe ensemble is built to span; the GIA trend field is smooth at "
        "continental scale, so the residual effect of the filter mismatch on "
        "an aquifer mean is small against the spread between GIA models and "
        "the no-correction level."
    )
    return xr.Dataset(
        data_vars={
            "ewh_abs": (
                dims, arr,
                {
                    "long_name": "equivalent water height (absolute)",
                    "units": "cm",
                    "comment": (
                        "Full synthesised EWH including the static gravity "
                        "field. Dynamic range is dominated by the static "
                        "field (~10⁹ cm); inter-variant differences are in "
                        "the trailing digits. For visualisation and most "
                        "analysis use ewh_anom; use ewh_abs for "
                        "GIA-model comparison and for rebasing anomalies "
                        "against an alternative reference period. " + common_comment
                    ),
                },
            ),
            "ewh_anom": (
                dims, arr_anom,
                {
                    "long_name": "equivalent water height anomaly",
                    "units": "cm",
                    "comment": (
                        "ewh_abs minus its time-mean over the reference "
                        "period (= reference_period_start..reference_period_end). "
                        "GIA axis is preserved in this view; at the IGP within "
                        "the 12-month POC the GIA spread across models is ~0.06 "
                        "cm against ~8 cm hydrology σ. " + common_comment
                    ),
                },
            ),
        },
        coords={
            "time": times,
            "truncation": list(cfg.TRUNCATIONS),
            "filter": list(cfg.FILTERS),
            "gia_model": list(cfg.GIA_MODELS),
            "c20_treatment": list(cfg.C20_TREATMENTS),
            "c30_treatment": list(cfg.C30_TREATMENTS),
            "geocenter": list(cfg.GEOCENTER_TREATMENTS),
            "lat": aoi.lat,
            "lon": aoi.lon,
        },
        attrs={
            "title": f"GRACE-FO EWH variant cube — {centre} — {region_label}",
            "centre": centre,
            "region": region,
            "aoi_bbox": f"{aoi.lat_min}N–{aoi.lat_max}N, {aoi.lon_min}E–{aoi.lon_max}E",
            "aoi_resolution_deg": aoi.res,
            "reference_period_start": str(times[0]),
            "reference_period_end": str(times[-1]),
            "reference_period_note": (
                "ewh_anom = ewh_abs − ewh_abs.mean('time'); reference period "
                "is the full POC window. For consumers wanting a longer "
                "baseline (e.g. 2004–2010) use ewh_abs and recompute."
            ),
            "gia_reference_epoch": str(cfg.GIA_REF_EPOCH),
            "gia_application": (
                "rate × (month_mid - gia_reference_epoch); "
                "see config.GIA_REF_EPOCH"
            ),
            "gia_axis_note": (
                "GIA axis is preserved in both ewh_abs and ewh_anom. At the "
                "IGP over the 12-month POC the spread across gia_model is "
                "~0.06 cm vs ~8 cm hydrology σ — small, but real and faithful."
            ),
            "Conventions": "CF-1.10",
        },
    )


# --------------------------------------------------------------------------
# CLI: python -m code.grace_pipeline.pipeline --centre CSR
# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Compute the variant-matrix EWH cube for one centre."
    )
    ap.add_argument(
        "--centre", required=True, choices=cfg.CENTRES,
        help="Processing centre (CSR / JPL / GFZ).",
    )
    ap.add_argument(
        "--filters", nargs="+", default=list(cfg.FILTERS),
        choices=cfg.FILTERS,
        help=(f"Filters to realise (default: all {len(cfg.FILTERS)}). "
              "Unrealised slots stay NaN in the cube."),
    )
    ap.add_argument(
        "--window", default="poc", choices=("poc", "appendix"),
        help="Analysis window. 'poc' = 89 GRACE-FO months (headline); "
             "'appendix' = 272 months full GRACE+GRACE-FO (supplement).",
    )
    ap.add_argument(
        "--region", default="igp", choices=tuple(cfg.AOI_REGISTRY),
        help="Study aquifer. 'igp' = Indo-Gangetic Plain (headline); "
             "'cv' = California Central Valley (second-aquifer case).",
    )
    args = ap.parse_args(argv)

    tag = f"{args.centre}/{args.region}/{args.window}"
    months_set = cfg.months_for_window(args.window)
    t0 = _time.time()
    print(f"[{tag}] starting compute_centre over "
          f"{len(months_set)} months × {cfg.variant_matrix_size_per_month()} "
          f"variants = {len(months_set) * cfg.variant_matrix_size_per_month()} "
          "fields", flush=True)
    ds, months = compute_centre(args.centre, filters=args.filters,
                                 window=args.window, region=args.region)
    elapsed = _time.time() - t0
    print(f"[{tag}] compute finished in {elapsed/60:.1f} min", flush=True)

    out_path = write_centre(ds, args.centre, months,
                              ran_filters=args.filters,
                              window=args.window, region=args.region)
    size_mb = out_path.stat().st_size / 1e6
    print(f"[{tag}] wrote {out_path} ({size_mb:.1f} MB)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
