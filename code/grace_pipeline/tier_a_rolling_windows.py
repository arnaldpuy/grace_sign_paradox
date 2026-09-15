"""Rolling-window sign stability of the recipe ensemble.

The Murray-Darling worked example shows a verdict that reverses between a
drought window and the full record. This script promotes that anecdote to a
cohort-level analysis, and separates two phenomena: trends unstable
because of PRE-PROCESSING (the recipe
ensemble disagrees within a window) and trends unstable because storage is
genuinely NON-MONOTONIC (the recipe-median sign itself flips across
windows). The two have different remedies and different policy implications.

Design: 10-year calendar windows slid annually (window w starts in April of
year y, 2002..2016, and ends in March of y+10), fitted on each recipe's
valid months within the window; a window is scored only when at least 84 of
its 120 nominal months are present (70%; the sparsest scored window, the
2009 start spanning the 2017-2018 mission bridge, retains 88 months). For every (aquifer, window) the 1,920
per-recipe OLS signs are pooled across the three centres:

  frac_pos   fraction of recipes with a positive trend in that window
  med_trend  cross-recipe median trend (cm/yr)
  spans      whether the window ensemble contains both signs

Output: datasets/output/tier_a/rolling_window_signs.csv
        (one row per aquifer x window, all 276 aquifers; Tier-1 selection
        happens R-side)

Run:  .venv python code/grace_pipeline/tier_a_rolling_windows.py
"""
from __future__ import annotations

import datetime as dt
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as cfg

REPO = cfg.REPO_ROOT
NPZ_DIR = cfg.REPO_ROOT / "datasets/output/tier_a_windows"
OUT = REPO / "datasets/output/tier_a/rolling_window_signs.csv"
CENTRES = ("csr", "jpl", "gfz")
WIN_YEARS = 10
START_YEARS = range(2002, 2017)
MIN_MONTHS = 84


def window_slopes(centre: str) -> dict[tuple[int, int], np.ndarray]:
    """{(start_year, n_aq-index-base): slopes} -> {start_year: (aq, rec)}"""
    npz = np.load(NPZ_DIR / f"monthly_{centre}.npz", allow_pickle=True)
    series = npz["series"]                      # (n_aq, n_rec, n_t)
    n_aq, n_rec, _ = series.shape
    jd = npz["dates_jd"].astype(float)
    epoch = dt.date(1970, 1, 1)
    dates = np.array([epoch + dt.timedelta(days=float(d)) for d in jd])
    out: dict[int, np.ndarray] = {}
    for y0 in START_YEARS:
        w_lo = dt.date(y0, 4, 1)
        w_hi = dt.date(y0 + WIN_YEARS, 4, 1)
        in_w = (dates >= w_lo) & (dates < w_hi)
        sub = series[:, :, in_w]                # (aq, rec, n_w) with NaNs
        t = (jd[in_w] - jd[in_w][0]) / 365.25
        finite = np.isfinite(sub)
        n_ok = finite.sum(axis=2)
        y_f = np.where(finite, sub, 0.0)
        t_b = t[None, None, :]
        t_f = np.where(finite, t_b, 0.0)
        mean_y = y_f.sum(axis=2) / np.maximum(n_ok, 1)
        mean_t = t_f.sum(axis=2) / np.maximum(n_ok, 1)
        y_c = np.where(finite, sub - mean_y[:, :, None], 0.0)
        t_c = np.where(finite, t_b - mean_t[:, :, None], 0.0)
        num = (y_c * t_c).sum(axis=2)
        den = (t_c * t_c).sum(axis=2)
        slope = np.where((n_ok >= MIN_MONTHS) & (den > 0),
                         num / np.maximum(den, 1e-30), np.nan)
        out[y0] = slope                          # (n_aq, n_rec)
    meta = {"aquifer_idx": npz["aquifer_idx"],
            "study_area": npz["study_area"]}
    return out, meta


def main() -> None:
    per_centre = {}
    meta = None
    for c in CENTRES:
        per_centre[c], meta = window_slopes(c)
        print(f"[{c}] {len(per_centre[c])} windows fitted", flush=True)

    rows = []
    for y0 in START_YEARS:
        pooled = np.concatenate([per_centre[c][y0] for c in CENTRES],
                                axis=1)          # (n_aq, 1920)
        n_fin = np.isfinite(pooled).sum(axis=1)
        with np.errstate(invalid="ignore"):
            n_pos = np.nansum(pooled > 0, axis=1)
            n_neg = np.nansum(pooled < 0, axis=1)
            med = np.nanmedian(np.where(np.isfinite(pooled), pooled,
                                        np.nan), axis=1)
        for a in range(pooled.shape[0]):
            if n_fin[a] == 0:
                continue
            rows.append({"aquifer_idx": int(meta["aquifer_idx"][a]),
                         "Study_area": meta["study_area"][a],
                         "win_start": y0,
                         "n_recipes": int(n_fin[a]),
                         "frac_pos": n_pos[a] / n_fin[a],
                         "med_trend_cmyr": med[a],
                         "spans": bool(n_pos[a] > 0 and n_neg[a] > 0)})
    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False, float_format="%.6g")
    n_win = df.groupby("Study_area")["win_start"].nunique()
    print(f"wrote {OUT}: {len(df):,} rows, "
          f"{df['Study_area'].nunique()} aquifers, "
          f"windows per aquifer {int(n_win.min())}-{int(n_win.max())}",
          flush=True)


if __name__ == "__main__":
    main()
