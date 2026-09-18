"""Rolling-window null + mission-window trends on the synthetic planets (D033).

The real-Earth analysis (code/grace_pipeline/
tier_a_rolling_windows.py + the `rolling-windows` chunk of
code_grace_uncertainty_clean.Rmd) finds that the recipe-median trend reverses
sign across 10-year windows for 68 of 73 Tier-1 aquifers, and that the
sign-paradox share falls from 77% on the full record to 26% on the GRACE-FO
window (SM10 Table: paradox robustness). On the synthetic planets the truth
is LINEAR (trend + seasonal cycle): every window flip and every window-to-
window paradox change is, by construction, noise + recipe arithmetic. This
script therefore provides the NULL distribution for both real results.

Mirrors the real code exactly:
  - 10-year calendar windows, April of y0 (2002..2016) -> March of y0+10,
    per-recipe OLS on valid months, window scored only when >= 84 of its 120
    nominal months are present (MIN_MONTHS, same as real);
  - per (aquifer, window): frac_pos / med_trend_cmyr / spans, pooled over the
    1,920 recipes (3 centres x 640);
  - mission windows exactly as tier_a_per_aquifer.py: w_grace <= 2017-06-30,
    w_fo >= 2018-06-01, w_full = everything. Per-recipe OLS on the aquifer-
    mean monthly series equals the pipeline's mean-of-cell-slopes for a fixed
    time design (equivalence proven in tier_a_window_trends.py; series parity
    vs the cubes verified at ~1e-7 in 04b).

Time coordinate gotcha (D022): xarray decodes time to datetime64[ns];
convert through datetime64[D] before casting to float.

Run (any planet subset):
  python code/07_rolling_windows_null.py

Outputs (datasets/output/rolling_null/):
  rolling_window_signs_<planet>.csv    one row per aquifer x window
  window_trends_<planet>_<win>.csv     Study_area x recipe trend (long)
"""
from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parent
OUT_DIR = ROOT / "datasets" / "output" / "rolling_null"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PLANETS = {
    "pilot_uniform_depletion_v2": "aquifer_series",
    "graded_floor_v2": "aquifer_series_graded_floor_v2",
    "graded_floor_v3": "aquifer_series_graded_floor_v3",
}
CENTRES = ("csr", "jpl", "gfz")
WIN_YEARS = 10
START_YEARS = range(2002, 2017)
MIN_MONTHS = 84

MISSION_WINDOWS = {
    "w_full": (None, None),
    "w_grace": (None, dt.date(2017, 6, 30)),
    "w_fo": (dt.date(2018, 6, 1), None),
}


def ols_slopes(series: np.ndarray, t_yr: np.ndarray,
               min_months: int = 2) -> np.ndarray:
    """Vectorised per-(aquifer, recipe) OLS slope over axis 2, NaN-aware.

    Same masked-OLS arithmetic as tier_a_rolling_windows.window_slopes.
    """
    finite = np.isfinite(series)
    n_ok = finite.sum(axis=2)
    y_f = np.where(finite, series, 0.0)
    t_b = t_yr[None, None, :]
    t_f = np.where(finite, t_b, 0.0)
    mean_y = y_f.sum(axis=2) / np.maximum(n_ok, 1)
    mean_t = t_f.sum(axis=2) / np.maximum(n_ok, 1)
    y_c = np.where(finite, series - mean_y[:, :, None], 0.0)
    t_c = np.where(finite, t_b - mean_t[:, :, None], 0.0)
    num = (y_c * t_c).sum(axis=2)
    den = (t_c * t_c).sum(axis=2)
    return np.where((n_ok >= min_months) & (den > 0),
                    num / np.maximum(den, 1e-30), np.nan)


def load_centre(series_dir: Path, centre: str):
    ds = xr.open_dataset(series_dir / f"aquifer_series_{centre}.nc")
    series = ds["ewh_anom"].values.astype(np.float64)   # (aq, rec, t)
    days = ds["time"].values.astype("datetime64[D]")
    dates = days.astype(dt.date)
    t_days = days.astype(float)
    names = [str(n) for n in ds["aquifer_name"].values]
    ds.close()
    return series, dates, t_days, names


def main() -> None:
    for planet, subdir in PLANETS.items():
        series_dir = ROOT / "datasets" / "output" / subdir
        per_centre = {}
        names = None
        for c in CENTRES:
            per_centre[c] = load_centre(series_dir, c)
            if names is None:
                names = per_centre[c][3]
            else:
                assert names == per_centre[c][3], f"aquifer order differs ({c})"
        n_aq = len(names)
        print(f"[{planet}] {n_aq} aquifers, "
              f"{per_centre['csr'][0].shape[1]} recipes/centre", flush=True)

        # Rolling 10-year windows ---------------------------------------------
        rows = []
        for y0 in START_YEARS:
            w_lo, w_hi = dt.date(y0, 4, 1), dt.date(y0 + WIN_YEARS, 4, 1)
            pooled = []
            for c in CENTRES:
                series, dates, t_days, _ = per_centre[c]
                in_w = np.array([(d >= w_lo) and (d < w_hi) for d in dates])
                t_yr = (t_days[in_w] - t_days[in_w][0]) / 365.25
                pooled.append(ols_slopes(series[:, :, in_w], t_yr,
                                         min_months=MIN_MONTHS))
            pooled = np.concatenate(pooled, axis=1)          # (aq, 1920)
            n_fin = np.isfinite(pooled).sum(axis=1)
            with np.errstate(invalid="ignore"):
                n_pos = np.nansum(pooled > 0, axis=1)
                n_neg = np.nansum(pooled < 0, axis=1)
                med = np.nanmedian(np.where(np.isfinite(pooled), pooled,
                                            np.nan), axis=1)
            for a in range(n_aq):
                if n_fin[a] == 0:
                    continue
                rows.append({"aquifer_idx": a, "Study_area": names[a],
                             "win_start": y0, "n_recipes": int(n_fin[a]),
                             "frac_pos": n_pos[a] / n_fin[a],
                             "med_trend_cmyr": med[a],
                             "spans": bool(n_pos[a] > 0 and n_neg[a] > 0)})
        df = pd.DataFrame(rows)
        out = OUT_DIR / f"rolling_window_signs_{planet}.csv"
        df.to_csv(out, index=False, float_format="%.6g")
        n_win = df.groupby("Study_area")["win_start"].nunique()
        print(f"  wrote {out.name}: {len(df):,} rows, windows/aquifer "
              f"{int(n_win.min())}-{int(n_win.max())}", flush=True)

        # Mission windows ------------------------------------------------------
        for wname, (lo, hi) in MISSION_WINDOWS.items():
            recs = []
            for c in CENTRES:
                series, dates, t_days, _ = per_centre[c]
                in_w = np.array([(lo is None or d >= lo) and
                                 (hi is None or d <= hi) for d in dates])
                t_yr = (t_days[in_w] - t_days[in_w][0]) / 365.25
                slopes = ols_slopes(series[:, :, in_w], t_yr)  # (aq, rec)
                n_rec = slopes.shape[1]
                recs.append(pd.DataFrame({
                    "centre": c.upper(),
                    "Study_area": np.repeat(names, n_rec),
                    "recipe_idx_per_centre": np.tile(np.arange(n_rec), n_aq),
                    "trend_cm_yr": slopes.ravel()}))
            wdf = pd.concat(recs, ignore_index=True)
            outw = OUT_DIR / f"window_trends_{planet}_{wname}.csv"
            wdf.to_csv(outw, index=False, float_format="%.8g")
            n_months = int(in_w.sum())
            print(f"  wrote {outw.name}: {wdf['trend_cm_yr'].notna().sum():,} "
                  f"trends ({wname}, {n_months} months, last centre)",
                  flush=True)


if __name__ == "__main__":
    main()
