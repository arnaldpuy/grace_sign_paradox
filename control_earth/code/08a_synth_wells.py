"""Synthetic wells for the validation-of-validation (Phase 7, D034).

METHODOLOGY §8: sample each planet's TRUE groundwater field at the real
Jasechko well locations, so 08b_val_of_val.R can convert storage to head via
a prescribed specific yield, add calibrated measurement noise, and re-run the
paper's above/below-boundary well-agreement test against known truth.

Two products per planet:
  1. synth_wells_<planet>.csv — one row per real well inside the 73 Tier-1
     aquifers: StnID, Study_area, and the true groundwater-layer trend at the
     well's 0.5 deg cell (cm/yr EWH). The gw trend grid is rebuilt exactly as
     01_truth_builder.py builds the gw layer (planted per-aquifer trend on
     the mask cells, then the D027 Gaussian edge taper; tapering is linear,
     so tapering the trend field equals the trend of the tapered series).
     Wells measure HEAD, i.e. the groundwater layer only — not the GLDAS
     soil-moisture/snow layer that GRACE also sees. The gap between the
     gw-layer trend and the total-truth answer key tau_true (D012) is
     therefore part of what the val-of-val measures.
  2. truth_gw_aquifer_<planet>.csv — mask-averaged gw-layer trend per
     aquifer (the well-accessible truth), alongside tau_true.

Plus one calibration product (planet-independent):
  well_noise_calibration.csv — per-aquifer SD of real per-well depth-to-water
  trends (OLS on annual values, 2002–2025, wells with >= 10 annual values),
  from Jasechko 2024's AnnualDepthToGroundwater.csv. The cohort median SD is
  the sigma_well used by 08b. This treats the real within-aquifer trend
  spread (genuine heterogeneity + measurement error) as i.i.d. per-well
  noise — a registered simplification (D034) that is conservative for the
  well-aggregate sign (real heterogeneity is spatially correlated and
  averages down more slowly).

Run:
  python code/08a_synth_wells.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parent
MASKS_CSV = REPO_ROOT / "datasets/jasechko_2024/cell_masks_0p5deg.csv"
WELLS_CSV = REPO_ROOT / "datasets/output/tier_a_full/wells_in_cohort.csv"  # not distributed
DTW_CSV = (REPO_ROOT / "datasets/jasechko_2024/groundwater_levels"  # not distributed
           / "AnnualDepthToGroundwater.csv")
OUT_DIR = ROOT / "datasets" / "output" / "val_of_val"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_LAT, N_LON = 340, 720          # the truth-builder global half-degree grid
LAT0 = -84.75

SPECS = {
    "pilot_uniform_depletion_v2": "truth_spec_pilot_v2.json",
    "graded_floor_v2": "truth_spec_graded_floor_v2.json",
    "graded_floor_v3": "truth_spec_graded_floor_v3.json",
}


def per_aq(spec_map: dict, aq_idx: int) -> float:
    return float(spec_map.get(str(aq_idx), spec_map["__default__"]))


def main() -> None:
    masks = pd.read_csv(MASKS_CSV)
    wells = pd.read_csv(WELLS_CSV)

    # Grid-convention check: the well's latitude must sit inside its cell ----
    cell_lat = LAT0 + 0.5 * wells["lat_idx"].to_numpy()
    assert np.abs(cell_lat - wells["Lat"].to_numpy()).max() <= 0.25 + 1e-9, \
        "wells_in_cohort lat_idx is not on the 340-lat truth grid"

    for planet, spec_name in SPECS.items():
        spec = json.loads((ROOT / "datasets/truth" / spec_name).read_text())
        taper = float(spec.get("gw_edge_taper_cells", 0.0))

        gw_trend = np.zeros((N_LAT, N_LON))
        for aq_idx, grp in masks.groupby("aquifer_idx"):
            tau = per_aq(spec["aquifer_trends_cm_yr"], aq_idx)
            gw_trend[grp["lat_idx"].to_numpy(), grp["lon_idx"].to_numpy()] = tau
        if taper > 0:
            gw_trend = gaussian_filter(gw_trend, sigma=taper, mode="nearest")

        cohort73 = set(masks["Study_area"].unique())
        w = wells[wells["Study_area"].isin(cohort73)].copy()
        w["true_gw_trend_cm_yr"] = gw_trend[w["lat_idx"].to_numpy(),
                                            w["lon_idx"].to_numpy()]
        out = OUT_DIR / f"synth_wells_{planet}.csv"
        w[["StnID", "Study_area", "lat_idx", "lon_idx",
           "true_gw_trend_cm_yr"]].to_csv(out, index=False,
                                          float_format="%.6g")

        # Mask-averaged gw-layer truth per aquifer (the well-accessible truth)
        rows = []
        for aq_idx, grp in masks.groupby("aquifer_idx"):
            tau_gw = gw_trend[grp["lat_idx"].to_numpy(),
                              grp["lon_idx"].to_numpy()].mean()
            rows.append({"aquifer_idx": int(aq_idx),
                         "Study_area": grp["Study_area"].iloc[0],
                         "tau_gw_cm_yr": tau_gw})
        gw_aq = pd.DataFrame(rows)
        tt = pd.read_csv(ROOT / "datasets/truth" /
                         f"truth_aquifer_trends_{planet}.csv")
        gw_aq = gw_aq.merge(tt[["aquifer_idx", "tau_true_cm_yr",
                                "tau_stipulated_cm_yr"]], on="aquifer_idx")
        gw_aq.to_csv(OUT_DIR / f"truth_gw_aquifer_{planet}.csv", index=False,
                     float_format="%.6g")
        n_sign_diff = int((np.sign(gw_aq["tau_gw_cm_yr"]) !=
                           np.sign(gw_aq["tau_true_cm_yr"])).sum())
        print(f"[{planet}] {len(w):,} wells in cohort; gw-vs-total truth "
              f"sign differs for {n_sign_diff}/73 aquifers", flush=True)

    # Well-noise calibration from the real per-well records -------------------
    dtw = pd.read_csv(DTW_CSV)
    dtw = dtw[(dtw["IntegerYear"] >= 2002) & (dtw["IntegerYear"] <= 2025)]
    dtw = dtw.merge(wells[["StnID", "Study_area"]], on="StnID", how="inner")
    g = dtw.groupby("StnID")
    counts = g["IntegerYear"].count()
    keep = counts[counts >= 10].index
    dtw = dtw[dtw["StnID"].isin(keep)]

    def well_slope(d: pd.DataFrame) -> float:
        x = d["IntegerYear"].to_numpy(float)
        y = d["DepthToWater_m"].to_numpy(float)
        x = x - x.mean()
        return float((x * (y - y.mean())).sum() / (x * x).sum())

    slopes = dtw.groupby(["Study_area", "StnID"]).apply(
        well_slope, include_groups=False).rename("trend_m_yr").reset_index()
    calib = slopes.groupby("Study_area").agg(
        n_wells=("trend_m_yr", "size"),
        sd_trend_m_yr=("trend_m_yr", "std"),
        median_trend_m_yr=("trend_m_yr", "median")).reset_index()
    calib.to_csv(OUT_DIR / "well_noise_calibration.csv", index=False,
                 float_format="%.6g")
    c10 = calib[calib["n_wells"] >= 10]
    print(f"[calibration] {len(slopes):,} real wells with >=10 annual values "
          f"(2002-2025) in {len(calib)} aquifers; cohort-median within-"
          f"aquifer SD = {c10['sd_trend_m_yr'].median():.4f} m/yr "
          f"({len(c10)} aquifers with >=10 such wells)", flush=True)


if __name__ == "__main__":
    main()
