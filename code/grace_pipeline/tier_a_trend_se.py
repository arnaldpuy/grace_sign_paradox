"""Per-recipe stochastic trend standard errors.

The detectability boundary is defined within the recipe ensemble: it marks
where PRE-PROCESSING spread overtakes the median trend. GRACE trends also
carry stochastic error (measurement noise, residual dealiasing) that the
ensemble does not sample. This script supplies the conventional counterpart:
for every (centre, aquifer, recipe) it fits the aquifer-mean monthly series

    y(t) = a + b t + c1 cos(2 pi t) + s1 sin(2 pi t)
                 + c2 cos(4 pi t) + s2 sin(4 pi t)

and reports the trend b with its standard error, inflated for lag-1
autocorrelation of the residuals (Santer et al. 2000): the effective sample
size n_eff = n (1 - rho1) / (1 + rho1) replaces n in the residual variance,
so se_ar1 = se_ols * sqrt((n - p) / (n_eff - p)), p = 6 regressors. rho1 is
computed over consecutive-month residual pairs only (gaps > 45 days - the
2017-2018 mission bridge and battery dropouts - are excluded from pairing)
and clamped at 0 from below (negative rho1 would deflate the error).

A plain slope-plus-intercept fit (b_plain) is emitted alongside as a parity
gate: it must reproduce the committed per-recipe trend tables to float
tolerance, proving the series and enumeration match the published numbers.

Input : analysis/datasets/output/tier_a_windows/monthly_{c}.npz
        (series[n_aq=276, n_recipes=640, n_t=272] float32, dates_jd,
        aquifer_idx, study_area, recipe_keys, recipe_cols; one per centre)
Output: datasets/output/tier_a/per_recipe_trend_se.csv.gz
        (529,920 rows: centre, aquifer, axes, trend_plain_cmyr,
        trend_seas_cmyr, se_raw, rho1, n_eff, se_ar1)

Run:  .venv python code/grace_pipeline/tier_a_trend_se.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import config as cfg

REPO = cfg.REPO_ROOT
NPZ_DIR = cfg.REPO_ROOT / "datasets/output/tier_a_windows"
OUT = REPO / "datasets/output/tier_a/per_recipe_trend_se.csv.gz"
PARITY = cfg.existing_table(
    REPO / "datasets/output/tier_a_full/jasechko_per_recipe_trends.csv")
CENTRES = ("csr", "jpl", "gfz")
P = 6  # regressors in the seasonal model


def fit_centre(centre: str) -> pd.DataFrame:
    npz = np.load(NPZ_DIR / f"monthly_{centre}.npz", allow_pickle=True)
    series = npz["series"]                      # (n_aq, n_rec, n_t)
    n_aq, n_rec, n_t_u = series.shape
    jd_u = npz["dates_jd"].astype(float)
    # The time axis is the cross-centre union; months a recipe lacks (centre
    # gaps, truncation-specific file availability) are NaN for every aquifer.
    # Verify aquifer-independence, group recipes by their NaN pattern, and
    # fit each group on its valid months only (matching _ols_slopes'
    # NaN-aware fit in tier_a_per_aquifer.py).
    masks_rec = np.isnan(series[0])             # (n_rec, n_t)
    if not (np.isnan(series) == masks_rec[None, :, :]).all():
        raise SystemExit(f"[{centre}] NaN pattern varies across aquifers")
    patterns, pat_of_rec = np.unique(masks_rec, axis=0, return_inverse=True)
    print(f"[{centre}] {len(patterns)} NaN pattern(s); valid months: "
          f"{[int(n_t_u - p.sum()) for p in patterns]}", flush=True)

    b_plain = np.empty((n_aq, n_rec)); b_seas = np.empty((n_aq, n_rec))
    se_raw = np.empty((n_aq, n_rec)); rho1_a = np.empty((n_aq, n_rec))
    neff_a = np.empty((n_aq, n_rec)); se_ar1 = np.empty((n_aq, n_rec))
    for p_i, pat in enumerate(patterns):
        rec = np.where(pat_of_rec == p_i)[0]
        ok = ~pat
        n_t = int(ok.sum())
        jd = jd_u[ok]
        t = (jd - jd[0]) / 365.25
        x = np.column_stack([np.ones(n_t), t,
                             np.cos(2 * np.pi * t), np.sin(2 * np.pi * t),
                             np.cos(4 * np.pi * t), np.sin(4 * np.pi * t)])
        y = series[:, rec][:, :, ok].reshape(n_aq * len(rec), n_t)
        y = y.T.astype(np.float64)              # (n_t, N)
        xtx_inv = np.linalg.inv(x.T @ x)
        beta = xtx_inv @ (x.T @ y)
        resid = y - x @ beta
        dof = n_t - P
        sigma2 = (resid ** 2).sum(axis=0) / dof
        se = np.sqrt(sigma2 * xtx_inv[1, 1])
        pair_ok = np.diff(jd) <= 45.0
        r0, r1 = resid[:-1][pair_ok], resid[1:][pair_ok]
        rho = (r0 * r1).sum(axis=0) / np.maximum(
            (resid ** 2).mean(axis=0) * pair_ok.sum(), 1e-30)
        rho = np.clip(rho, 0.0, 0.99)
        neff = n_t * (1.0 - rho) / (1.0 + rho)
        sea = se * np.sqrt(dof / np.maximum(neff - P, 3.0))
        t_c = t - t.mean()
        bp = (t_c @ (y - y.mean(axis=0))) / (t_c @ t_c)
        sh = (n_aq, len(rec))
        b_plain[:, rec] = bp.reshape(sh); b_seas[:, rec] = beta[1].reshape(sh)
        se_raw[:, rec] = se.reshape(sh); rho1_a[:, rec] = rho.reshape(sh)
        neff_a[:, rec] = neff.reshape(sh); se_ar1[:, rec] = sea.reshape(sh)

    keys = pd.DataFrame(np.asarray(npz["recipe_keys"]),
                        columns=list(npz["recipe_cols"]))
    df = pd.DataFrame({
        "centre": centre.upper(),
        "aquifer_idx": np.repeat(npz["aquifer_idx"], n_rec),
        "Study_area": np.repeat(npz["study_area"], n_rec),
        "trend_plain_cmyr": b_plain.ravel(),
        "trend_seas_cmyr": b_seas.ravel(),
        "se_raw": se_raw.ravel(), "rho1": rho1_a.ravel(),
        "n_eff": np.round(neff_a.ravel(), 1), "se_ar1": se_ar1.ravel()})
    for c in keys.columns:
        df[c] = np.tile(keys[c].values, n_aq)
    print(f"[{centre}] {len(df):,} fits, median rho1 "
          f"{np.median(rho1_a):.2f}, median se_ar1 "
          f"{np.median(se_ar1):.3f} cm/yr", flush=True)
    return df


def main() -> None:
    out = pd.concat([fit_centre(c) for c in CENTRES], ignore_index=True)

    committed = pd.read_csv(PARITY)
    axes = ["truncation", "filter", "gia_model", "c20_treatment",
            "c30_treatment", "geocenter"]
    committed["truncation"] = committed["truncation"].astype(str)
    out["truncation"] = out["truncation"].astype(str)
    m = committed.merge(out[["centre", "aquifer_idx", *axes,
                             "trend_plain_cmyr"]],
                        on=["centre", "aquifer_idx", *axes], how="inner")
    if len(m) != len(committed):
        raise SystemExit(f"parity merge incomplete: {len(m):,} of "
                         f"{len(committed):,} committed rows matched")
    gap = float(np.nanmax(np.abs(m["trend_cm_yr"] - m["trend_plain_cmyr"])))
    print(f"parity gate: {len(m):,} rows, max |plain-slope - committed| "
          f"= {gap:.2e} cm/yr", flush=True)
    if gap > 1e-4:
        raise SystemExit("parity gate FAILED - enumeration mismatch")

    out.to_csv(OUT, index=False, float_format="%.6g")
    print(f"wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
