"""Attribute the rolling-window sign flips on the synthetic planets (D061).

The planted groundwater trends are linear, but the total-storage truth also
carries the GLDAS-Noah background with its interannual variability, so a
recovered window flip is not spurious by construction. This script scores
every 10-year window against the known truth and, on the noise-replicate
pair, separates the flips into three sources:

  hydrology        the truth's own window trend changes sign;
  shared           the recovered recipe-median sign is wrong against a truth
                   whose sign is stable, in BOTH members of the noise-replicate
                   pair;
  replicate-only   wrong against a stable truth in ONE member only.

The split is descriptive. An error shared by two noise realisations is not
thereby shown to be free of noise, and an error in one member only is not
thereby shown to be free of processing; the two classes are reported as what
they are, not as an identified partition into causes.

Inputs (per planet): truth_tws_<planet>.nc (GRACE_CE_CUBE_DIR), aquifer_series_<centre>.nc
(GRACE_CE_CUBE_DIR) or the banked rolling_window_signs_<planet>.csv, cell_masks_0p5deg.csv.
Windows, months threshold and OLS arithmetic mirror 07_rolling_windows_null.py.

Outputs (datasets/output/rolling_null/):
  rolling_window_signs_<planet>.csv          (new: graded_balanced_v3, _noise2)
  rolling_truth_window_trends.csv            planet x aquifer x window truth trend
  rolling_attribution_planets.csv            per-planet flip attribution
  rolling_attribution_replicates.csv         replicate-pair attribution
"""
from __future__ import annotations
import datetime as dt
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr

ROOT = Path(__file__).resolve().parents[1]
import os
CUBE_DIR = Path(os.environ.get("GRACE_CE_CUBE_DIR", ROOT / "datasets" / "synthetic_l2"))  # planet folders with truth_tws_<planet>.nc and aquifer_series/
OUT = ROOT / "datasets" / "output" / "rolling_null"
REPO_ROOT = ROOT.parent
MASKS = REPO_ROOT / "datasets" / "jasechko_2024" / "cell_masks_0p5deg.csv"
PLANETS = ["pilot_uniform_depletion_v2", "graded_floor_v2", "graded_floor_v3",
           "graded_balanced_v3", "graded_balanced_v3_noise2"]
CENTRES = ("csr", "jpl", "gfz"); WIN_YEARS = 10; START_YEARS = range(2002, 2017); MIN_MONTHS = 84

def ols_slopes(series, t_yr, min_months=2):
    finite = np.isfinite(series); n_ok = finite.sum(axis=-1)
    y_f = np.where(finite, series, 0.0); t_f = np.where(finite, t_yr, 0.0)
    mean_y = y_f.sum(-1) / np.maximum(n_ok, 1); mean_t = t_f.sum(-1) / np.maximum(n_ok, 1)
    y_c = np.where(finite, series - mean_y[..., None], 0.0); t_c = np.where(finite, t_yr - mean_t[..., None], 0.0)
    num = (y_c * t_c).sum(-1); den = (t_c * t_c).sum(-1)
    return np.where((n_ok >= min_months) & (den > 0), num / np.maximum(den, 1e-30), np.nan)

def windows():
    return [(y0, dt.date(y0, 4, 1), dt.date(y0 + WIN_YEARS, 4, 1)) for y0 in START_YEARS]

def recovered_signs(planet):
    """Recipe-median sign per aquifer x window (banked CSV or recomputed from the series)."""
    p = OUT / f"rolling_window_signs_{planet}.csv"
    if p.exists():
        return pd.read_csv(p)
    per_centre = {}
    for c in CENTRES:
        ds = xr.open_dataset(CUBE_DIR / planet / "aquifer_series" / f"aquifer_series_{c}.nc")
        per_centre[c] = (ds["ewh_anom"].values.astype(np.float64), ds["time"].values.astype("datetime64[D]"),
                         [str(n) for n in ds["aquifer_name"].values]); ds.close()
    names = per_centre["csr"][2]; rows = []
    for y0, lo, hi in windows():
        pooled = []
        for c in CENTRES:
            series, days, _ = per_centre[c]
            dates = days.astype(dt.date); in_w = np.array([(d >= lo) and (d < hi) for d in dates])
            t_yr = (days[in_w].astype(float) - days[in_w].astype(float)[0]) / 365.25
            pooled.append(ols_slopes(series[:, :, in_w], t_yr, MIN_MONTHS))
        pooled = np.concatenate(pooled, axis=1)
        for a, nm in enumerate(names):
            v = pooled[a][np.isfinite(pooled[a])]
            rows.append(dict(aquifer_idx=a, Study_area=nm, win_start=y0, n_recipes=len(v),
                             frac_pos=np.mean(v > 0) if len(v) else np.nan,
                             med_trend_cmyr=np.median(v) if len(v) else np.nan,
                             spans=bool(len(v) and (v > 0).any() and (v < 0).any())))
    df = pd.DataFrame(rows); df.to_csv(p, index=False); print(f"[{planet}] banked {p.name}"); return df

def truth_window_trends(planet, masks):
    ds = xr.open_dataset(CUBE_DIR / planet / f"truth_tws_{planet}.nc")
    lat = ds["lat"].values; lon = ds["lon"].values
    t = pd.to_datetime(ds["time"].values); tws = ds["tws_anom"]
    rows = []
    for nm, g in masks.groupby("Study_area", sort=False):
        li = g["lat_idx"].values; lo = g["lon_idx"].values
        assert np.allclose(lat[li], g["lat"].values) and np.allclose(lon[lo], g["lon"].values), nm
        s = np.nanmean(tws.values[:, li, lo], axis=1) if False else None
        sub = tws.isel(lat=xr.DataArray(li, dims="c"), lon=xr.DataArray(lo, dims="c")).values  # (t, c)
        s = np.nanmean(sub, axis=1)
        for y0, w_lo, w_hi in windows():
            in_w = (t >= pd.Timestamp(w_lo)) & (t < pd.Timestamp(w_hi))
            if in_w.sum() < MIN_MONTHS: continue
            days = (t[in_w] - t[in_w][0]).days.values.astype(float) / 365.25
            rows.append(dict(planet=planet, Study_area=nm, win_start=y0, n_months=int(in_w.sum()),
                             truth_trend_cmyr=float(ols_slopes(s[in_w], days, MIN_MONTHS))))
    ds.close(); return pd.DataFrame(rows)

def flips(sign_by_win):
    """True when the sign changes at least once across the ordered windows (zeros ignored)."""
    s = np.array([x for x in sign_by_win if x != 0]); return len(s) > 1 and (s != s[0]).any()

masks = pd.read_csv(MASKS)
truth_all, planet_rows, per_planet = [], [], {}
for planet in PLANETS:
    rec = recovered_signs(planet); tru = truth_window_trends(planet, masks); truth_all.append(tru)
    d = rec.merge(tru, on=["Study_area", "win_start"], how="inner")
    d["rec_sign"] = np.sign(d["med_trend_cmyr"]); d["truth_sign"] = np.sign(d["truth_trend_cmyr"])
    d["wrong"] = d["rec_sign"] != d["truth_sign"]
    per_planet[planet] = d
    aq = d.sort_values("win_start").groupby("Study_area").agg(
        rec_flip=("rec_sign", flips), truth_flip=("truth_sign", flips),
        n_wrong=("wrong", "sum"), n_win=("wrong", "size"), n_wrong_spanning=("wrong", lambda x: int((x & d.loc[x.index, "spans"]).sum())))
    planet_rows.append(dict(planet=planet, n_aq=len(aq), rec_flip=int(aq.rec_flip.sum()), truth_flip=int(aq.truth_flip.sum()),
        rec_flip_truth_flip=int((aq.rec_flip & aq.truth_flip).sum()), rec_flip_truth_stable=int((aq.rec_flip & ~aq.truth_flip).sum()),
        aq_any_wrong_window=int((aq.n_wrong > 0).sum()), windows=int(aq.n_win.sum()), wrong_windows=int(aq.n_wrong.sum()),
        wrong_windows_spanning=int(aq.n_wrong_spanning.sum()),
        median_frac_wrong=float((aq.n_wrong / aq.n_win).median())))
pd.concat(truth_all).to_csv(OUT / "rolling_truth_window_trends.csv", index=False)
pl = pd.DataFrame(planet_rows); pl.to_csv(OUT / "rolling_attribution_planets.csv", index=False)
pd.set_option("display.width", 220); print("\n=== per-planet attribution"); print(pl.to_string(index=False))

# Replicate pair ---------------------------------------------------------------
a = per_planet["graded_balanced_v3"]; b = per_planet["graded_balanced_v3_noise2"]
m = a.merge(b, on=["Study_area", "win_start"], suffixes=("_a", "_b"))
assert np.allclose(m.truth_trend_cmyr_a, m.truth_trend_cmyr_b)   # same truth
m["wrong_both"] = m.wrong_a & m.wrong_b; m["wrong_one"] = m.wrong_a ^ m.wrong_b
aq = m.sort_values("win_start").groupby("Study_area").agg(
    truth_flip=("truth_sign_a", flips), flip_a=("rec_sign_a", flips), flip_b=("rec_sign_b", flips),
    wrong_both=("wrong_both", "sum"), wrong_one=("wrong_one", "sum"), n_win=("wrong_a", "size"),
    sign_agree=("rec_sign_a", lambda s: float((s.values == m.loc[s.index, "rec_sign_b"].values).mean())))
rep = dict(n_aq=len(aq), truth_flip=int(aq.truth_flip.sum()), flip_a=int(aq.flip_a.sum()), flip_b=int(aq.flip_b.sum()),
           flip_both=int((aq.flip_a & aq.flip_b).sum()), flip_one=int((aq.flip_a ^ aq.flip_b).sum()),
           flip_a_truth_stable=int((aq.flip_a & ~aq.truth_flip).sum()), flip_b_truth_stable=int((aq.flip_b & ~aq.truth_flip).sum()),
           flip_both_truth_stable=int((aq.flip_a & aq.flip_b & ~aq.truth_flip).sum()),
           windows=int(aq.n_win.sum()), wrong_both_windows=int(aq.wrong_both.sum()), wrong_one_windows=int(aq.wrong_one.sum()),
           replicate_sign_agreement=float(m.rec_sign_a.eq(m.rec_sign_b).mean()))
pd.DataFrame([rep]).to_csv(OUT / "rolling_attribution_replicates.csv", index=False)
print("\n=== replicate pair (graded_balanced_v3 vs noise2)"); print(pd.Series(rep).to_string())
# truth stable + recovered flips: how many aquifers' truth window trend never changes sign?
print("\n=== truth |window trend| quantiles per planet (cm/yr)")
for p, d in per_planet.items(): print(p, np.round(np.quantile(np.abs(d.truth_trend_cmyr), [0.1, 0.5, 0.9]), 3))
