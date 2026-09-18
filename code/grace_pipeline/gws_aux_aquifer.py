"""Per-aquifer GLDAS-Noah auxiliary storage (soil moisture + SWE + canopy)
series and month-matched trends, for the TWS -> GWS sensitivity analysis.

Rationale:
the 1,920-recipe ensemble propagates pre-processing choices into TWS,
while published GRACE groundwater claims subtract soil moisture, snow and
surface water from TWS using auxiliary land-surface models. The auxiliary
field is recipe-INDEPENDENT, so for OLS trends fitted on the same months

    trend(GWS_r) = trend(TWS_r) - trend(AUX)        for every recipe r,

which means the GWS ensemble can be obtained from the existing per-recipe
per-aquifer TWS trend table by subtracting ONE auxiliary trend per
(aquifer x centre month-set) - no cube pass needed. Per-aquifer recipe
IQRs are unchanged exactly; medians (and hence sign counts, dominance R
and the detectability boundary) shift.

This script computes that auxiliary trend, mirroring the conventions of
tier_a_per_aquifer.py exactly:
  * auxiliary = GLDAS-2.1 Noah monthly 0.25 deg: 4 soil-moisture layers
    (0-200 cm) + SWE + canopy water, kg/m^2 (= mm) summed and /10 -> cm
    (same composition as grace_pipeline/gws_correction.py; surface water
    is NOT included - documented limitation, CDEC covers only the
    Central Valley);
  * spatial aggregation: plain (unweighted) mean over the aquifer's mask
    cells (tier_a_per_aquifer.py averages per-cell slopes without area
    weights; with a common month axis the mean of per-cell OLS slopes
    equals the OLS slope of the cell-mean series, so the cell-mean series
    is fitted here). Each 0.5 deg mask cell maps to its 2x2 GLDAS 0.25 deg
    subcells; ocean/no-data subcells are NaN and excluded;
  * time axis: per centre, the months actually present in that centre's
    appendix cube (matched by calendar YYYY-MM), with
    t = (jd - jd[0]) / 365.25 years - the same convention as
    tier_a_per_aquifer._ols_slopes, so the subtraction is exact.

Inputs:
  * datasets/jasechko_2024/cell_masks_0p5deg_276.csv  (covers Tier 1 + 2)
  * $GRACE_DATA_ROOT/raw/gldas/GLDAS_NOAH025_M.A{YYYYMM}.021.nc4
  * $GRACE_DATA_ROOT/processed/grace/global/
        grace_ewh_global_{csr,jpl,gfz}_appendix.nc   (time axis only)

Outputs (datasets/output/tier_a/):
  * gws_aux_gldas_monthly.csv  - aquifer x month auxiliary storage (cm)
  * gws_aux_trends.csv         - aquifer x centre month-matched OLS trend
                                 (cm/yr) of the auxiliary series

A second, independent auxiliary (--source era5land) reads the ERA5-Land
monthly means downloaded by download_era5land_aux.py: volumetric soil
water layers 1-4 (x layer thickness 0.07/0.21/0.72/1.89 m -> cm), snow
water equivalent and skin-reservoir (canopy) content (m -> cm). The
difference between the two per-aquifer auxiliary trends measures the
variance added by the auxiliary-model choice itself. Geometry note: the
regridded ERA5 nodes sit AT multiples of 0.25 deg, so each 0.5 deg mask
cell takes the 2x2 nodes at offsets (-0.25, 0) per axis (exact match,
no rounding); the GLDAS grid is offset-centred and uses (+-0.125).

Run:  .venv python, from the repo root:
  python code/grace_pipeline/gws_aux_aquifer.py [--source gldas|era5land]
"""
from __future__ import annotations

import argparse
import datetime as dt
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

from . import config as cfg

REPO = cfg.REPO_ROOT
GLDAS_DIR = cfg.DATA_ROOT / "raw" / "gldas"
ERA5_DIR = cfg.DATA_ROOT / "raw" / "era5_land"
CUBE_DIR = cfg.OUT_DIR / "global"
MASKS = REPO / "datasets/jasechko_2024/cell_masks_0p5deg_276.csv"
OUT_DIR = REPO / "datasets/output/tier_a"
CENTRES = ("csr", "jpl", "gfz")

GLDAS_VARS = ("SoilMoi0_10cm_inst", "SoilMoi10_40cm_inst",
              "SoilMoi40_100cm_inst", "SoilMoi100_200cm_inst",
              "SWE_inst", "CanopInt_inst")
# ERA5-Land: variable -> factor converting to cm of column water
ERA5_VARS = {"swvl1": 0.07 * 100, "swvl2": 0.21 * 100,
             "swvl3": 0.72 * 100, "swvl4": 1.89 * 100,
             "sd": 100.0, "src": 100.0}


def build_aquifer_subcells(masks: pd.DataFrame, glat: np.ndarray,
                           glon: np.ndarray,
                           offsets: tuple[float, float]) -> list[dict]:
    """Map each aquifer's 0.5 deg mask cells to 0.25 deg grid indices."""
    lat0, dlat = float(glat[0]), float(glat[1] - glat[0])
    lon0, dlon = float(glon[0]), float(glon[1] - glon[0])
    aquifers = []
    for aq_idx, grp in masks.groupby("aquifer_idx"):
        la, lo = [], []
        for lat_c, lon_c in zip(grp["lat"].values, grp["lon"].values):
            for dla in offsets:
                for dlo in offsets:
                    i = int(round((lat_c + dla - lat0) / dlat))
                    j = int(round((lon_c + dlo - lon0) / dlon))
                    if 0 <= i < len(glat):
                        la.append(i)
                        lo.append(j % len(glon))
        aquifers.append({"aquifer_idx": int(aq_idx),
                         "Study_area": grp["Study_area"].iloc[0],
                         "lat_idx": np.array(la), "lon_idx": np.array(lo)})
    return aquifers


def _aquifer_rows(aquifers: list[dict], total: np.ndarray,
                  month: str) -> list[dict]:
    rows = []
    for a in aquifers:
        vals = total[a["lat_idx"], a["lon_idx"]]
        ok = np.isfinite(vals)
        rows.append({"aquifer_idx": a["aquifer_idx"],
                     "Study_area": a["Study_area"], "month": month,
                     "aux_cm": float(np.nanmean(vals)) if ok.any()
                     else np.nan,
                     "n_subcells_ok": int(ok.sum())})
    return rows


def monthly_aux_series_gldas() -> pd.DataFrame:
    """One pass over the GLDAS granules -> per-aquifer cell-mean series."""
    files = sorted(GLDAS_DIR.glob("GLDAS_NOAH025_M.A*.021.nc4"))
    if not files:
        raise SystemExit(f"no GLDAS granules under {GLDAS_DIR}")
    with netCDF4.Dataset(files[0]) as nc0:
        glat = nc0.variables["lat"][:].filled(np.nan)
        glon = nc0.variables["lon"][:].filled(np.nan)
    masks = pd.read_csv(MASKS)
    aquifers = build_aquifer_subcells(masks, glat, glon, (-0.125, 0.125))
    print(f"{len(aquifers)} aquifers, "
          f"{sum(len(a['lat_idx']) for a in aquifers)} GLDAS subcells, "
          f"{len(files)} months", flush=True)

    rows = []
    for k, f in enumerate(files):
        ym = f.name.split(".A")[1][:6]
        month = f"{ym[:4]}-{ym[4:6]}"
        with netCDF4.Dataset(f) as nc:
            total = None
            for v in GLDAS_VARS:
                arr = nc.variables[v][0].filled(np.nan)  # (lat, lon), mm
                total = arr if total is None else total + arr
        total /= 10.0  # mm -> cm
        rows.extend(_aquifer_rows(aquifers, total, month))
        if (k + 1) % 48 == 0 or k == len(files) - 1:
            print(f"  {k + 1}/{len(files)} months done", flush=True)
    return pd.DataFrame(rows)


def monthly_aux_series_era5land() -> pd.DataFrame:
    """ERA5-Land monthly means -> per-aquifer cell-mean series (cm)."""
    files = sorted(ERA5_DIR.glob("era5_land_aux_*.nc"))
    if not files:
        raise SystemExit(f"no ERA5-Land files under {ERA5_DIR} -- run "
                         "download_era5land_aux.py first")
    masks = pd.read_csv(MASKS)
    rows, aquifers = [], None
    for f in files:
        with netCDF4.Dataset(f) as nc:
            glat = nc.variables["latitude"][:].filled(np.nan)
            glon = nc.variables["longitude"][:].filled(np.nan)
            if aquifers is None:
                aquifers = build_aquifer_subcells(masks, glat, glon,
                                                  (-0.25, 0.0))
                print(f"{len(aquifers)} aquifers, "
                      f"{sum(len(a['lat_idx']) for a in aquifers)} "
                      "ERA5-Land subcells", flush=True)
            tname = "valid_time" if "valid_time" in nc.variables else "time"
            tvar = nc.variables[tname]
            dates = netCDF4.num2date(tvar[:], tvar.units,
                                     getattr(tvar, "calendar", "standard"))
            missing = [v for v in ERA5_VARS if v not in nc.variables]
            if missing:
                raise SystemExit(f"{f.name}: variables missing {missing}; "
                                 f"present: {sorted(nc.variables)}")
            for ti, d in enumerate(dates):
                month = f"{d.year:04d}-{d.month:02d}"
                total = None
                for v, fac in ERA5_VARS.items():
                    arr = nc.variables[v][ti].filled(np.nan) * fac
                    total = arr if total is None else total + arr
                rows.extend(_aquifer_rows(aquifers, total, month))
        print(f"  {f.name} done", flush=True)
    return pd.DataFrame(rows)


NPZ_DIR = cfg.REPO_ROOT / "datasets/output/tier_a_windows"


def valid_month_sets(centre: str) -> list[tuple[list[int], list[str],
                                                np.ndarray]]:
    """(truncations, months, jd) per distinct valid-month set of a centre.

    The per-recipe trends were fitted NaN-aware on each recipe's valid
    months (the cube time axis is the cross-centre union; months a centre
    or a truncation lacks are NaN). The auxiliary trend must be fitted on
    the SAME months for the subtraction trend(GWS) = trend(TWS) -
    trend(AUX) to be exact. The valid-month pattern is truncation-specific
    for JPL (251 months at lmax 60, 242 at lmax 96) and common for
    CSR/GFZ (250); patterns are read from the monthly_<centre>.npz series
    emitted by tier_a_window_trends.py.
    """
    npz = np.load(NPZ_DIR / f"monthly_{centre}.npz", allow_pickle=True)
    masks = np.isnan(npz["series"][0])          # (n_rec, n_t)
    if not (np.isnan(npz["series"]) == masks[None, :, :]).all():
        raise SystemExit(f"[{centre}] NaN pattern varies across aquifers")
    keys = np.asarray(npz["recipe_keys"])
    ti = list(npz["recipe_cols"]).index("truncation")
    jd_all = npz["dates_jd"].astype(float)
    epoch = dt.date(1970, 1, 1)
    out = []
    patterns, inv = np.unique(masks, axis=0, return_inverse=True)
    for p_i, pat in enumerate(patterns):
        truncs = sorted({int(v) for v in keys[inv == p_i, ti]})
        jd = jd_all[~pat]
        months = [(epoch + dt.timedelta(days=float(d))).strftime("%Y-%m")
                  for d in jd]
        out.append((truncs, months, jd))
    return out


def month_matched_trends(series: pd.DataFrame) -> pd.DataFrame:
    """OLS trend of the auxiliary on each (centre, truncation) month set."""
    wide = series.pivot_table(index="month", columns="aquifer_idx",
                              values="aux_cm")
    meta = series.drop_duplicates("aquifer_idx").set_index("aquifer_idx")
    rows = []
    for centre in CENTRES:
        for truncs, months, jd in valid_month_sets(centre):
            missing = [m for m in months if m not in wide.index]
            if missing:
                raise SystemExit(f"[{centre}] months missing from the "
                                 f"auxiliary series: {missing}")
            t_years = (jd - jd[0]) / 365.25
            sub = wide.loc[months]  # (n_t, n_aq), valid-month order
            for aq_idx in wide.columns:
                y = sub[aq_idx].values
                ok = np.isfinite(y)
                if ok.sum() >= 5:
                    t_c = t_years[ok] - t_years[ok].mean()
                    y_c = y[ok] - y[ok].mean()
                    slope = float((y_c * t_c).sum() / (t_c * t_c).sum())
                else:
                    slope = np.nan
                for tr in truncs:
                    rows.append({"aquifer_idx": int(aq_idx),
                                 "Study_area": meta.loc[aq_idx,
                                                        "Study_area"],
                                 "centre": centre.upper(),
                                 "truncation": tr,
                                 "aux_trend_cmyr": slope,
                                 "n_months": int(ok.sum())})
            print(f"[{centre}] truncation {truncs}: {len(months)} valid "
                  "months matched", flush=True)
    return pd.DataFrame(rows)


def window_trends(source: str) -> None:
    """Auxiliary trends per (paper-window, aquifer, centre, truncation).

    For the GWS-space scoring of the literature claims:
    the per-window per-recipe TWS trend tables under
    datasets/output/tier_a_windows/ were fitted on each recipe's valid
    GRACE months inside the window, so the auxiliary must be fitted on the
    SAME months for the per-window subtraction to be exact. Windows are
    discovered from the table file names
    (jasechko_per_recipe_trends__<start>_<end>.csv); the monthly auxiliary
    series must already exist (run the default mode first).
    """
    monthly_csv = ("gws_aux_gldas_monthly.csv" if source == "gldas"
                   else "gws_aux_era5land_monthly.csv")
    series = pd.read_csv(OUT_DIR / monthly_csv)
    wide = series.pivot_table(index="month", columns="aquifer_idx",
                              values="aux_cm")
    meta = series.drop_duplicates("aquifer_idx").set_index("aquifer_idx")
    win_dir = REPO / "datasets/output/tier_a_windows"
    wins = sorted(p.name.split("__")[1][:-4]
                  for p in win_dir.glob("jasechko_per_recipe_trends__*.csv"))
    print(f"{len(wins)} paper windows discovered", flush=True)

    epoch = dt.date(1970, 1, 1)
    rows = []
    for centre in CENTRES:
        npz = np.load(NPZ_DIR / f"monthly_{centre}.npz", allow_pickle=True)
        masks = np.isnan(npz["series"][0])
        keys = np.asarray(npz["recipe_keys"])
        ti = list(npz["recipe_cols"]).index("truncation")
        jd_all = npz["dates_jd"].astype(float)
        months_all = np.array(
            [(epoch + dt.timedelta(days=float(d))).strftime("%Y-%m")
             for d in jd_all])
        patterns, inv = np.unique(masks, axis=0, return_inverse=True)
        for p_i, pat in enumerate(patterns):
            truncs = sorted({int(v) for v in keys[inv == p_i, ti]})
            for win in wins:
                w_lo, w_hi = win.split("_")
                ok = (~pat) & (months_all >= w_lo) & (months_all <= w_hi)
                if ok.sum() < 5:
                    continue
                jd = jd_all[ok]
                t_years = (jd - jd[0]) / 365.25
                t_c = t_years - t_years.mean()
                sub = wide.loc[months_all[ok]]      # (n_t, n_aq)
                y_c = sub.values - sub.values.mean(axis=0, keepdims=True)
                slope = (t_c @ y_c) / (t_c @ t_c)
                for tr in truncs:
                    for aq_idx, s in zip(wide.columns, slope):
                        rows.append({"data_window": win,
                                     "aquifer_idx": int(aq_idx),
                                     "Study_area": meta.loc[aq_idx,
                                                            "Study_area"],
                                     "centre": centre.upper(),
                                     "truncation": tr,
                                     "aux_trend_cmyr": float(s),
                                     "n_months": int(ok.sum())})
        print(f"[{centre}] window auxiliary trends done", flush=True)
    out = pd.DataFrame(rows)
    target = OUT_DIR / f"gws_aux_window_trends_{source}.csv"
    out.to_csv(target, index=False, float_format="%.6g")
    print(f"wrote {target} ({len(out):,} rows)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=("gldas", "era5land"),
                    default="gldas")
    ap.add_argument("--windows", action="store_true",
                    help="emit per-paper-window auxiliary trends instead")
    args = ap.parse_args()
    src = args.source
    if args.windows:
        window_trends(src)
        return
    if src == "gldas":
        series = monthly_aux_series_gldas()
        monthly_csv, trends_csv = ("gws_aux_gldas_monthly.csv",
                                   "gws_aux_trends.csv")
    else:
        series = monthly_aux_series_era5land()
        monthly_csv, trends_csv = ("gws_aux_era5land_monthly.csv",
                                   "gws_aux_trends_era5land.csv")
    n_dry = series["aux_cm"].isna().sum()
    if n_dry:
        bad = series[series["aux_cm"].isna()]["Study_area"].unique()
        print(f"WARNING: {n_dry} aquifer-months all-NaN ({bad})")
    series.to_csv(OUT_DIR / monthly_csv, index=False)
    trends = month_matched_trends(series)
    trends.to_csv(OUT_DIR / trends_csv, index=False)
    print(f"\nwrote {OUT_DIR / monthly_csv}")
    print(f"wrote {OUT_DIR / trends_csv}")
    t1 = trends[trends["aquifer_idx"] < 73]
    print(f"\nTier-1 auxiliary trend summary ({src}, cm/yr):")
    print(t1.groupby("centre")["aux_trend_cmyr"]
            .describe()[["count", "min", "25%", "50%", "75%", "max"]])


if __name__ == "__main__":
    main()
