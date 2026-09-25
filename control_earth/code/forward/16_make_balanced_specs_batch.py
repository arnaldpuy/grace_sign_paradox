"""16_make_balanced_specs_batch.py -- specs for balanced layouts 2, 3 and the
sign-reversed counterpart of layout 1 (AMENDMENT_01_sign_balance.md section 3;
D051).

Layouts 2 and 3 repeat the D048 construction of graded_balanced_v1 with new
seeds: same per-aquifer |tau| as graded_floor_v3, signs balanced and
stratified across |tau| quartiles, deterministic sub-seed rejection forcing
|spearman(sign, |tau|)| and |spearman(sign, area)| below 0.10. The spec seed
also drives the forward operator's noise, so each layout carries a fresh
noise realisation.

The reversed planet flips EVERY sign of graded_balanced_v1 and keeps its seed
(20260805), so the forward operator reproduces layout 1's noise exactly:
identical spatial arrangement, magnitudes and noise, opposite directions.

  python3 code/forward/16_make_balanced_specs_batch.py   # writes + verifies
"""

from __future__ import annotations

import json
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
REPO_ROOT = ROOT.parent

v3 = json.loads((ROOT / "datasets/truth/truth_spec_graded_floor_v3.json")
                .read_text())
mags = {int(k): abs(v) for k, v in v3["aquifer_trends_cm_yr"].items()
        if k != "__default__"}
assert len(mags) == 73

masks = pd.read_csv(REPO_ROOT / "datasets/jasechko_2024/cell_masks_0p5deg.csv")
n_cells = masks.groupby("aquifer_idx").size()


def signs_of(spec: dict) -> dict:
    return {k: np.sign(v) for k, v in spec["aquifer_trends_cm_yr"].items()
            if k != "__default__"}


v1_signs = signs_of(json.loads(
    (ROOT / "datasets/truth/truth_spec_graded_balanced_v1.json").read_text()))
written_signs = {"graded_balanced_v1": v1_signs}


def make_layout(planet: str, seed: int) -> None:
    """The D048 construction verbatim, parameterised by planet id and seed."""
    df = pd.DataFrame({"aq": list(mags), "mag": [mags[a] for a in mags]})
    df["stratum"] = pd.qcut(df["mag"].rank(method="first"), 4, labels=False)
    df["n_cells"] = df["aq"].map(n_cells)

    for sub in range(1000):
        rng = np.random.default_rng(seed + sub)
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
        raise RuntimeError(f"[{planet}] no sub-seed satisfied the gates")
    df["tau"] = df["sign"] * df["mag"]

    n_pos, n_neg = (df["sign"] > 0).sum(), (df["sign"] < 0).sum()
    assert abs(n_pos - n_neg) <= 1
    print(f"[{planet}] sub-seed {seed}+{sub}; signs: {n_pos} pos / {n_neg} neg; "
          f"spearman(sign,|tau|) = {rho_mag:+.3f}; "
          f"spearman(sign,area) = {rho_area:+.3f}")

    spec = dict(v3)
    spec["planet_id"] = planet
    spec["seed"] = seed
    spec["comment"] = (
        f"Balanced sign layout ({planet}, amendment 01 / D051). Same "
        "construction as graded_balanced_v1 (D048) with a new seed: identical "
        "per-aquifer |tau| to graded_floor_v3, signs balanced and stratified "
        "across |tau| quartiles, independence of magnitude and area enforced "
        "by sub-seed rejection. Fresh noise realisation (the spec seed drives "
        "the forward operator). Built by "
        "code/forward/16_make_balanced_specs_batch.py.")
    spec["aquifer_trends_cm_yr"] = {"__default__": 0.0,
                                    **{str(a): round(t, 4) for a, t in
                                       zip(df["aq"], df["tau"])}}
    # A layout must be a NEW sign assignment: reject any duplicate (or exact
    # flip) of layout 1 or of a layout already written. Layout 1 accepted
    # sub-seed 20260805+1, so a base seed of 20260806 silently reproduces its
    # RNG stream -- this happened on the first attempt and is why the bases
    # below are separated by more than the rejection loop's range.
    new_signs = signs_of(spec)
    for other, osigns in written_signs.items():
        n_same = sum(new_signs[k] == osigns[k] for k in new_signs)
        assert 0 < n_same < 73, f"{planet} duplicates (or exactly flips) {other}"
    written_signs[planet] = new_signs

    out = ROOT / f"datasets/truth/truth_spec_{planet}.json"
    out.write_text(json.dumps(spec, indent=1))
    print(f"  wrote {out.name}; tau range "
          f"{df['tau'].min():+.2f} .. {df['tau'].max():+.2f}")


def make_reversed() -> None:
    """Flip every sign of graded_balanced_v1; keep its seed (same noise)."""
    v1 = json.loads((ROOT / "datasets/truth/truth_spec_graded_balanced_v1.json")
                    .read_text())
    planet = "graded_balanced_v1_reversed"
    spec = dict(v1)
    spec["planet_id"] = planet
    spec["comment"] = (
        "Sign-reversed counterpart of graded_balanced_v1 (amendment 01 / "
        "D051): every aquifer trend negated, seed UNCHANGED (20260805) so the "
        "forward operator reproduces layout 1's noise exactly. Identical "
        "spatial arrangement, magnitudes and noise; any recovery difference "
        "is a pure advantage of one sign direction over the other. Built by "
        "code/forward/16_make_balanced_specs_batch.py.")
    spec["aquifer_trends_cm_yr"] = {
        k: -v for k, v in v1["aquifer_trends_cm_yr"].items()}

    # The reversal must differ from layout 1 ONLY in tau signs, id and comment.
    assert spec["seed"] == v1["seed"] == 20260805
    same = {k: v for k, v in spec.items()
            if k not in ("planet_id", "comment", "aquifer_trends_cm_yr")}
    assert same == {k: v for k, v in v1.items()
                    if k not in ("planet_id", "comment", "aquifer_trends_cm_yr")}
    taus = spec["aquifer_trends_cm_yr"]
    assert all(taus[k] == -v for k, v in v1["aquifer_trends_cm_yr"].items())
    assert taus["__default__"] == 0.0

    out = ROOT / f"datasets/truth/truth_spec_{planet}.json"
    out.write_text(json.dumps(spec, indent=1))
    vals = [v for k, v in taus.items() if k != "__default__"]
    print(f"[{planet}] seed kept at {spec['seed']}; "
          f"{sum(v > 0 for v in vals)} pos / {sum(v < 0 for v in vals)} neg; "
          f"tau range {min(vals):+.2f} .. {max(vals):+.2f}")
    print(f"  wrote {out.name}")


def make_noise_replicate(src: str, planet: str, seed: int) -> None:
    """Same sign layout and magnitudes as `src`, new forward-operator noise
    (amendment 01 section 5: the stopping rule tripped on 2026-08-07, D055,
    so one same-layout different-noise replicate attributes the layout-to-
    layout spread to noise vs sign geometry)."""
    base = json.loads((ROOT / f"datasets/truth/truth_spec_{src}.json").read_text())
    spec = dict(base)
    spec["planet_id"] = planet
    spec["seed"] = seed
    spec["comment"] = (
        f"Noise replicate of {src} (amendment 01 section 5; D055/D056): "
        "IDENTICAL aquifer trends (same signs, same magnitudes), new seed so "
        "the forward operator draws a fresh noise realisation. If this planet "
        "scores like its source, the layout-to-layout spread is sign "
        "geometry; if it recovers toward layout 1, the spread is noise. "
        "Built by code/forward/16_make_balanced_specs_batch.py.")
    assert spec["aquifer_trends_cm_yr"] == base["aquifer_trends_cm_yr"]
    assert spec["seed"] != base["seed"]
    out = ROOT / f"datasets/truth/truth_spec_{planet}.json"
    out.write_text(json.dumps(spec, indent=1))
    print(f"[{planet}] trends identical to {src}; noise seed {base['seed']} "
          f"-> {seed}")
    print(f"  wrote {out.name}")


# Bases separated by >> 1000 (the sub-seed loop range) from layout 1's
# 20260805 and from each other, so no two layouts can share an RNG stream.
make_layout("graded_balanced_v2", seed=20270805)
make_layout("graded_balanced_v3", seed=20280805)
make_reversed()
make_noise_replicate("graded_balanced_v3", "graded_balanced_v3_noise2",
                     seed=20290805)
