"""Global trend-estimator invariance of the sign paradox.

The main text reports that ordinary least squares, Theil-Sen and seasonal
OLS agree on the sign of every recipe trend. That check was run at one
Indo-Gangetic grid cell. This script repeats it for the whole Tier-1 cohort:
73 aquifers x 1,920 recipes, refitting each aquifer-mean monthly series with
all three estimators and comparing signs.

Estimators
  ols      slope of a straight line in time
  theilsen median of all pairwise slopes (robust to outliers)
  seasonal OLS with annual and semi-annual harmonics, so that the trend is
           not absorbed by an aliased seasonal cycle

Input (banked, no cube pass needed):
  datasets/output/tier_a_windows/monthly_{csr,gfz,jpl}.npz
    series      (n_aquifers, n_recipes_per_centre, n_months) aquifer-mean EWH
    dates_jd    (n_months,) Julian days
    recipe_keys (n_recipes_per_centre, 6) axis levels
  datasets/jasechko_2024/cohort_50K_n73.csv   Tier-1 membership

Output:
  datasets/output/robustness/trend_estimator_global.csv
    one row per (aquifer, recipe): the three slopes and their signs
  plus a printed summary: recipes whose sign differs between estimators, and
  whether any aquifer changes sign-paradox status under a different estimator.
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
from scipy import stats

# Run from the repository root. The monthly NPZ archives (~160 MB per
# centre) are not shipped; build them from the recipe cubes with
# code/grace_pipeline/tier_a_window_trends.py (its step 1 writes
# monthly_<centre>.npz into datasets/output/tier_a_windows/).
WIN = pathlib.Path("datasets/output/tier_a_windows")
COHORT = pathlib.Path("datasets/jasechko_2024/cohort_50K_n73.csv")
OUT_DIRS = (pathlib.Path("datasets/output/robustness"),)
CENTRES = ("csr", "gfz", "jpl")


def seasonal_ols(t: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Trend of y on t with annual + semi-annual harmonics, columnwise."""
    x = np.column_stack([
        np.ones_like(t), t,
        np.sin(2 * np.pi * t), np.cos(2 * np.pi * t),
        np.sin(4 * np.pi * t), np.cos(4 * np.pi * t)])
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    return coef[1]


def main() -> int:
    tier1 = set(pd.read_csv(COHORT)["Study_area"])
    rows = []
    for centre in CENTRES:
        z = np.load(WIN / f"monthly_{centre}.npz", allow_pickle=True)
        series = z["series"]                       # (aq, recipe, month)
        study = np.asarray(z["study_area"], dtype=object)
        keys = np.asarray(z["recipe_keys"], dtype=object)
        cols = [str(c) for c in z["recipe_cols"]]
        t = z["dates_jd"].astype(float) / 365.25
        t = t - t.mean()

        keep = np.array([s in tier1 for s in study])
        print(f"[{centre}] {keep.sum()} Tier-1 aquifers x {series.shape[1]} recipes",
              flush=True)

        for ai in np.flatnonzero(keep):
            y = series[ai].astype(float)           # (recipe, month)
            ok = np.isfinite(y).all(axis=0)
            tt, yy = t[ok], y[:, ok]
            # OLS and seasonal OLS solve every recipe at once.
            b_ols = np.polyfit(tt, yy.T, 1)[0]
            b_sea = seasonal_ols(tt, yy.T)
            b_ts = np.array([stats.theilslopes(r, tt)[0] for r in yy])
            for ri in range(yy.shape[0]):
                rows.append({
                    "Study_area": study[ai], "centre": centre.upper(),
                    **{c: keys[ri, j] for j, c in enumerate(cols)},
                    "ols": b_ols[ri], "theilsen": b_ts[ri],
                    "seasonal": b_sea[ri]})
        del z, series

    df = pd.DataFrame(rows)
    for d in OUT_DIRS:
        d.mkdir(parents=True, exist_ok=True)
        df.to_csv(d / "trend_estimator_global.csv.gz", index=False,
                  compression="gzip")

    # --- sign agreement across estimators ------------------------------------
    s_ols = np.sign(df["ols"]); s_ts = np.sign(df["theilsen"])
    s_sea = np.sign(df["seasonal"])
    agree = (s_ols == s_ts) & (s_ols == s_sea)
    print(f"\nrecipe-level series           : {len(df):,}")
    print(f"all three estimators agree    : {agree.sum():,} "
          f"({100 * agree.mean():.2f}%)")
    print(f"  OLS vs Theil-Sen disagree   : {(s_ols != s_ts).sum():,}")
    print(f"  OLS vs seasonal disagree    : {(s_ols != s_sea).sum():,}")

    # --- does the paradox classification change? -----------------------------
    print("\nsign-paradoxical aquifers (of 73) by estimator:")
    for est in ("ols", "theilsen", "seasonal"):
        g = df.groupby("Study_area")[est]
        spans = ((g.min() < 0) & (g.max() > 0)).sum()
        print(f"  {est:9s}: {spans}")
    piv = df.groupby("Study_area").agg(
        ols_sp=("ols", lambda v: (v.min() < 0) & (v.max() > 0)),
        ts_sp=("theilsen", lambda v: (v.min() < 0) & (v.max() > 0)),
        sea_sp=("seasonal", lambda v: (v.min() < 0) & (v.max() > 0)))
    flip = piv[(piv.ols_sp != piv.ts_sp) | (piv.ols_sp != piv.sea_sp)]
    print(f"aquifers whose paradox status depends on the estimator: {len(flip)}")
    if len(flip):
        print(flip)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
