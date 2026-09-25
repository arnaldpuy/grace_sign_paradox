"""Common-month sensitivity of the headline sign paradox (Supplementary
Materials, Table S1).

The three centres do not see the same months (CSR/GFZ 250 valid months, JPL
251 at degree 60 and 242 at degree 96), so the per-recipe trends of the
headline table are fitted on slightly different samplings. This script refits
every (aquifer, recipe) trend of the 73-aquifer cohort on the months valid in
EVERY recipe of every centre, from the banked per-aquifer monthly series
(tier_a_window_trends.py --emit-monthly), so the notebook can recompute the
sign-paradox count and the R = 1 boundary on identical sampling.

Output: datasets/output/robustness/common_months_per_recipe_trends.csv.gz
  centre, aquifer_idx (73-cohort index), Study_area, truncation, filter,
  gia_model, c20_treatment, c30_treatment, geocenter, trend_cm_yr, n_months
Validation: refitting on each recipe's own valid months reproduces the
canonical 73-aquifer table (max |diff| printed).

Usage
-----
    cd <repository>/code
    python -m grace_pipeline.tier_a_common_months

Needs datasets/output/tier_a_windows/monthly_<centre>.npz (about 170 MB per
centre, not distributed; one pass of tier_a_window_trends.py --emit-monthly
over the three recipe cubes rebuilds them).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from . import config as cfg

REPO = cfg.REPO_ROOT
DEFAULT_NPZ_DIR = REPO / "datasets" / "output" / "tier_a_windows"
DEFAULT_MASKS_CSV = REPO / "datasets" / "jasechko_2024" / "cell_masks_0p5deg.csv"
DEFAULT_CANON_PR = cfg.existing_table(
    REPO / "datasets" / "output" / "tier_a_full_73" / "jasechko_per_recipe_trends.csv")
DEFAULT_OUT = (REPO / "datasets" / "output" / "robustness"
               / "common_months_per_recipe_trends.csv.gz")
CENTRES = ("csr", "jpl", "gfz")
KEY = ["centre", "Study_area", "truncation", "filter", "gia_model",
       "c20_treatment", "c30_treatment", "geocenter"]


def ols(y: np.ndarray, t: np.ndarray) -> np.ndarray:
    """NaN-aware OLS slope along the last axis (t broadcast over y)."""
    ok = np.isfinite(y)
    n = ok.sum(-1)
    yf = np.where(ok, y, 0.0)
    tf = np.where(ok, t, 0.0)
    my = yf.sum(-1) / n
    mt = tf.sum(-1) / n
    yc = np.where(ok, y - my[..., None], 0.0)
    tc = np.where(ok, t - mt[..., None], 0.0)
    return (yc * tc).sum(-1) / (tc * tc).sum(-1)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--npz_dir", type=Path, default=DEFAULT_NPZ_DIR)
    ap.add_argument("--masks_csv", type=Path, default=DEFAULT_MASKS_CSV)
    ap.add_argument("--canon_pr", type=Path, default=DEFAULT_CANON_PR)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)

    cohort = pd.read_csv(args.masks_csv).groupby("Study_area")["aquifer_idx"].first()
    npz = {c: np.load(args.npz_dir / f"monthly_{c}.npz", allow_pickle=True)
           for c in CENTRES}

    # Months valid in every recipe of every centre (the validity pattern is
    # identical across aquifers: gaps are whole-month mission gaps).
    common = None
    for z in npz.values():
        fin = np.isfinite(z["series"][0])
        common = fin.all(0) if common is None else (common & fin.all(0))
    print("common months:", int(common.sum()), "of", len(common), flush=True)

    canon = pd.read_csv(args.canon_pr)
    canon["centre"] = canon["centre"].str.lower()
    rows, val = [], []
    for c, z in npz.items():
        sa = [str(s) for s in z["study_area"]]
        keys = pd.DataFrame(z["recipe_keys"], columns=list(z["recipe_cols"]))
        t_all = (z["dates_jd"] - z["dates_jd"][0]) / 365.25
        idx = [i for i, s in enumerate(sa) if s in cohort.index]
        series = z["series"][idx].astype(np.float64)        # (73, n_recipes, n_t)
        tr_common = ols(series[:, :, common], t_all[common])
        tr_own = ols(series, t_all)
        for a, i in enumerate(idx):
            df = keys.copy()
            df.insert(0, "Study_area", sa[i])
            df.insert(0, "aquifer_idx", int(cohort[sa[i]]))
            df.insert(0, "centre", c.upper())
            df["trend_cm_yr"] = tr_common[a]
            df["n_months"] = int(common.sum())
            rows.append(df)
            v = keys.copy()
            v["centre"] = c
            v["Study_area"] = sa[i]
            v["trend_own"] = tr_own[a]
            val.append(v)

    out = pd.concat(rows, ignore_index=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False, compression="gzip")
    print("wrote", args.out, len(out), "rows", flush=True)

    # Validation: own-month refit against the canonical table.
    v = pd.concat(val, ignore_index=True)
    v["truncation"] = v["truncation"].astype(int)
    canon["truncation"] = canon["truncation"].astype(int)
    m = canon.merge(v, on=KEY, how="inner")
    print("validation rows", len(m), "max|diff| own-months vs canonical:",
          float((m.trend_cm_yr - m.trend_own).abs().max()), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
