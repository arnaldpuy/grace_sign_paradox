"""13_make_balanced_spec.py -- spec for the sign-balanced graded planet (D048).

The 0.75 cm/yr recovery floor is fitted on graded_floor_v3, where 69 of 73
aquifers are negative: an always-depleting rule scores 95% there, so coherent
same-sign leakage could flatter recovery at low magnitudes. This planet keeps
graded_floor_v3's per-aquifer |tau| assignment EXACTLY (magnitudes already
permuted independently of area, D032) and re-randomises only the SIGNS,
balanced 37/36 and stratified across |tau| quartiles so that sign is
independent of magnitude by construction. Refitting the floor here isolates
the sign-imbalance effect: the two planets differ in nothing else.

Background stays GLDAS-2.1 Noah (the standard chain; the crossed planet of
D045 covers the auxiliary question separately).

  python3 code/forward/13_make_balanced_spec.py        # writes + verifies
"""

from __future__ import annotations

import json
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
REPO_ROOT = ROOT.parent
SEED = 20260805
PLANET = "graded_balanced_v1"

v3 = json.loads((ROOT / "datasets/truth/truth_spec_graded_floor_v3.json")
                .read_text())
mags = {int(k): abs(v) for k, v in v3["aquifer_trends_cm_yr"].items()
        if k != "__default__"}
assert len(mags) == 73

# Stratified sign assignment: quartiles of |tau|, half positive within each.
# Stratification guarantees sign is independent of magnitude; a deterministic
# rejection loop over sub-seeds additionally forces near-zero correlation with
# aquifer area, so sign carries no accidental structure at all.
masks = pd.read_csv(REPO_ROOT / "datasets/jasechko_2024/cell_masks_0p5deg.csv")
n_cells = masks.groupby("aquifer_idx").size()

df = pd.DataFrame({"aq": list(mags), "mag": [mags[a] for a in mags]})
df["stratum"] = pd.qcut(df["mag"].rank(method="first"), 4, labels=False)
df["n_cells"] = df["aq"].map(n_cells)

for sub in range(1000):
    rng = np.random.default_rng(SEED + sub)
    df["sign"] = -1
    for s, grp in df.groupby("stratum"):
        n_pos = len(grp) // 2 + (1 if (s % 2 == 0 and len(grp) % 2) else 0)
        pos = rng.choice(grp.index.to_numpy(), size=n_pos, replace=False)
        df.loc[pos, "sign"] = 1
    rho_mag = df[["sign", "mag"]].corr(method="spearman").iloc[0, 1]
    rho_area = df[["sign", "n_cells"]].corr(method="spearman").iloc[0, 1]
    if abs(rho_mag) < 0.10 and abs(rho_area) < 0.10:
        break
else:
    raise RuntimeError("no sub-seed satisfied the independence gates")
df["tau"] = df["sign"] * df["mag"]

n_pos, n_neg = (df["sign"] > 0).sum(), (df["sign"] < 0).sum()
per_stratum = df.groupby("stratum")["sign"].apply(lambda s: (s > 0).sum())
print(f"[{PLANET}] sub-seed {SEED}+{sub}; signs: {n_pos} pos / {n_neg} neg; "
      f"per-quartile pos: {per_stratum.tolist()}")
print(f"  spearman(sign, |tau|) = {rho_mag:+.3f}; "
      f"spearman(sign, n_cells) = {rho_area:+.3f}")
assert abs(n_pos - n_neg) <= 1

spec = dict(v3)
spec["planet_id"] = PLANET
spec["seed"] = SEED
spec["comment"] = (
    "Sign-balanced graded planet (D048). Identical to graded_floor_v3 in "
    "every respect -- same per-aquifer |tau| (already permuted independently "
    "of area), same seasonal terms, edge taper and accepted noise parameters, "
    "same GLDAS-2.1 Noah background -- EXCEPT that the trend signs are "
    "re-randomised, balanced 37/36 and stratified across |tau| quartiles so "
    "sign is independent of magnitude and area by construction. Purpose: the "
    "recovery floor of D035 is fitted where 69 of 73 aquifers are negative; "
    "refitting it here tests whether coherent same-sign leakage flattered "
    "recovery at low magnitudes. Built by code/forward/13_make_balanced_spec.py.")
spec["aquifer_trends_cm_yr"] = {"__default__": 0.0,
                                **{str(a): round(t, 4) for a, t in
                                   zip(df["aq"], df["tau"])}}
out = ROOT / f"datasets/truth/truth_spec_{PLANET}.json"
out.write_text(json.dumps(spec, indent=1))
print(f"  wrote {out.name}; tau range "
      f"{df['tau'].min():+.2f} .. {df['tau'].max():+.2f}")
