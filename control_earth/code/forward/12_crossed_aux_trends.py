"""12_crossed_aux_trends.py -- the answer key for the groundwater subtraction
on a crossed-truth planet (D045).

`truth_aquifer_trends_<planet>.csv` scores the TOTAL storage trend, which on a
crossed planet is contaminated by the background model's own trends: the field
is `truth = background + gw`, so `tau_true` mixes the planted groundwater
signal with ERA5-Land hydrology. Scoring the SUBTRACTION needs the groundwater
layer's own aquifer-mean trend, and the two auxiliary trends an analyst would
remove.

Three quantities per aquifer, all on the planet's own 272 months, its mask set
and its OLS convention (identical to 01_truth_builder.py):

  tau_gw     trend of the tapered groundwater layer alone -- THE answer key
             for the subtraction. Below the stipulated trend in magnitude
             because the 1-cell Gaussian edge taper (D027) bleeds signal out
             of the polygon.
  aux_era5   trend of the planet's own background. The planet was built from
             ERA5-Land with exactly the composition of the ERA5-Land
             auxiliary (ERA5_VARS), so this is the RIGHT model: recovered
             minus this isolates the filter/attenuation mismatch. Obtained as
             tau_true - tau_gw, which is exact because OLS is linear and
             truth = background + gw.
  aux_gldas  trend of the GLDAS-2.1 Noah auxiliary over the same months --
             the WRONG model, and what the real analysis actually subtracts.
             Loaded from the raw GLDAS files, snow-tower zeroed on its own
             field exactly as the background was (D013).

Run (about 6 min, GLDAS dominates):

  python code/forward/12_crossed_aux_trends.py \
      --planet crossed_era5_v1
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time as _time

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from importlib import import_module  # noqa: E402

_tb = import_module("01_truth_builder")

CONTROL_ROOT = _tb.CONTROL_ROOT
TRUTH_DIR = _tb.TRUTH_DIR


def ols_slope(t_yr: np.ndarray, series: np.ndarray) -> float:
    """Same estimator as the truth builder's answer key (D012)."""

    tm = t_yr - t_yr.mean()
    return float((tm * (series - series.mean())).sum() / (tm ** 2).sum())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--planet", default="crossed_era5_v1")
    args = ap.parse_args()
    planet = args.planet

    spec = pd.read_json(TRUTH_DIR / f"truth_spec_{planet}.json",
                        typ="series")
    masks = pd.read_csv(_tb.MASKS_CSV)
    months = _tb.real_sampled_months()
    n_t = len(months)
    aoi = _tb.GLOBAL_AOI

    t_ref = np.datetime64(spec["trend_ref_epoch"])
    t_yr = np.array([(m.time_mid - t_ref) / np.timedelta64(1, "D")
                     for m in months]) / 365.25
    month_frac = np.array([m.month for m in months], dtype=float)
    print(f"[aux {planet}] {n_t} months, {masks.aquifer_idx.nunique()} aquifers")

    # ---- groundwater layer alone, tapered exactly as the truth builder -----

    def per_aq(spec_map: dict, aq_idx: int) -> float:
        return float(spec_map.get(str(aq_idx), spec_map["__default__"]))

    gw = np.zeros((n_t, aoi.lat.size, aoi.lon.size), dtype=np.float64)
    for aq_idx, grp in masks.groupby("aquifer_idx"):
        tau = per_aq(spec["aquifer_trends_cm_yr"], aq_idx)
        amp = per_aq(spec["seasonal_amp_cm"], aq_idx)
        phase = per_aq(spec["seasonal_phase_month"], aq_idx)
        series = tau * t_yr + amp * np.sin(
            2 * np.pi * (month_frac - phase) / 12.0)
        gw[:, grp["lat_idx"].to_numpy(), grp["lon_idx"].to_numpy()] = \
            series[:, None]

    taper = float(spec.get("gw_edge_taper_cells", 0.0))
    if taper > 0:
        from scipy.ndimage import gaussian_filter
        for i in range(gw.shape[0]):
            gw[i] = gaussian_filter(gw[i], sigma=taper, mode="nearest")
        print(f"[aux {planet}] groundwater layer edge-tapered "
              f"(Gaussian, {taper} cells; D027)")

    # ---- total truth, to recover the background by difference --------------

    with xr.open_dataset(TRUTH_DIR / f"truth_tws_{planet}.nc") as ds:
        truth = ds["tws_anom"].to_numpy().astype(np.float64)
    assert truth.shape == gw.shape, (truth.shape, gw.shape)

    # ---- GLDAS-Noah auxiliary over the same months -------------------------

    t0 = _time.time()
    gldas = np.zeros_like(gw)
    for i, m in enumerate(months):
        gldas[i] = _tb.to_global_half_degree(
            _tb.load_gldas_month(m.year, m.month))
        if (i + 1) % 50 == 0:
            print(f"  gldas {i+1}/{n_t} ({(_time.time()-t0)/60:.1f} min)")
    gldas -= gldas.mean(axis=0, keepdims=True)
    thr = float(spec.get("snow_tower_threshold_cm", 100.0))
    tower = (np.abs(gldas) > thr).any(axis=0)
    gldas[:, tower] = 0.0
    print(f"[aux {planet}] gldas done ({(_time.time()-t0)/60:.1f} min); "
          f"{int(tower.sum())} snow-tower cells zeroed (D013)")

    # ---- per-aquifer trends -------------------------------------------------

    rows = []
    for aq_idx, grp in masks.groupby("aquifer_idx"):
        li = grp["lat_idx"].to_numpy()
        lo = grp["lon_idx"].to_numpy()
        tau_true = ols_slope(t_yr, truth[:, li, lo].mean(axis=1))
        tau_gw = ols_slope(t_yr, gw[:, li, lo].mean(axis=1))
        rows.append({
            "aquifer_idx": int(aq_idx),
            "Study_area": grp["Study_area"].iloc[0],
            "n_cells": len(grp),
            "tau_true_cm_yr": tau_true,
            "tau_gw_cm_yr": tau_gw,
            "tau_stipulated_cm_yr": per_aq(spec["aquifer_trends_cm_yr"],
                                           aq_idx),
            "aux_era5_cm_yr": tau_true - tau_gw,
            "aux_gldas_cm_yr": ols_slope(t_yr, gldas[:, li, lo].mean(axis=1)),
        })
    out = pd.DataFrame(rows)

    # ---- verification against the frozen answer key -------------------------

    key = pd.read_csv(TRUTH_DIR / f"truth_aquifer_trends_{planet}.csv")
    merged = out.merge(key[["aquifer_idx", "tau_true_cm_yr"]],
                       on="aquifer_idx", suffixes=("", "_key"))
    dmax = float(np.abs(merged.tau_true_cm_yr
                        - merged.tau_true_cm_yr_key).max())
    print(f"[aux {planet}] max |tau_true - banked key| = {dmax:.2e} cm/yr")
    assert dmax < 1e-6, "reconstruction does not reproduce the banked key"

    out_csv = TRUTH_DIR / f"truth_gw_aux_trends_{planet}.csv"
    out.to_csv(out_csv, index=False)
    print(f"[aux {planet}] wrote {out_csv.name} ({len(out)} aquifers)")
    print(f"  tau_gw       {out.tau_gw_cm_yr.min():+.3f} .. "
          f"{out.tau_gw_cm_yr.max():+.3f}  (taper retains "
          f"{(out.tau_gw_cm_yr / out.tau_stipulated_cm_yr).median():.2f} "
          f"of the stipulated trend, median)")
    print(f"  aux_era5     {out.aux_era5_cm_yr.min():+.3f} .. "
          f"{out.aux_era5_cm_yr.max():+.3f}")
    print(f"  aux_gldas    {out.aux_gldas_cm_yr.min():+.3f} .. "
          f"{out.aux_gldas_cm_yr.max():+.3f}")
    disc = out.aux_gldas_cm_yr - out.aux_era5_cm_yr
    print(f"  gldas - era5 median |discrepancy| {disc.abs().median():.3f} "
          f"cm/yr, max {disc.abs().max():.3f}")


if __name__ == "__main__":
    main()
