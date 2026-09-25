"""Per-recipe exterior-ring trends (Supplementary Materials, leakage test).

tier_a_leakage.py banks per-aquifer recipe-ensemble METRICS on three masks.
The notebook then inferred the ring trend from the buffered and original
medians with the cell-count identity, which holds for one recipe's means but
not for medians across recipes. This script writes the quantity that identity
was standing in for: for every (centre, aquifer, recipe), the mean per-cell
OLS slope over the ORIGINAL cells and over the one-cell RING (buffered minus
original), so the ring-minus-interior contrast can be summarised across
recipes directly.

Masks, connectivity, cube streaming and the NaN-aware OLS are imported from
tier_a_leakage.py unchanged (the original-mask trends reproduce
jasechko_per_recipe_trends.csv, checked below).

Output:
    datasets/output/robustness/leakage_ring_per_recipe.csv.gz
      centre, aquifer_idx, Study_area, recipe_idx_per_centre, n_cells_original,
      n_cells_ring, trend_original, trend_ring, ring_contrast (= ring - original)

Usage
-----
    cd <repository>/code
    python -m grace_pipeline.tier_a_leakage_ring [--centres csr jpl gfz]

Needs the three per-centre recipe cubes (track 3 of the README); about 35 min.
"""
from __future__ import annotations
import argparse, sys, time as _time
from pathlib import Path
import numpy as np, pandas as pd, netCDF4 as nc4
from . import config as cfg
from .tier_a_leakage import (build_masks, slopes_for_centre, DEFAULT_CUBE_DIR,
                             DEFAULT_MASKS_CSV, DEFAULT_CANON_PR)

DEFAULT_OUT = (cfg.REPO_ROOT / "datasets" / "output" / "robustness"
               / "leakage_ring_per_recipe.csv.gz")

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--cube_dir", type=Path, default=DEFAULT_CUBE_DIR)
    ap.add_argument("--masks_csv", type=Path, default=DEFAULT_MASKS_CSV)
    ap.add_argument("--canon_pr", type=Path, default=DEFAULT_CANON_PR)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--centres", nargs="+", default=["csr", "jpl", "gfz"])
    args = ap.parse_args(argv)
    masks_df = pd.read_csv(args.masks_csv)
    with nc4.Dataset(args.cube_dir / f"grace_ewh_global_{args.centres[0]}_appendix.nc", "r") as nc:
        nlat, nlon = len(nc.dimensions["lat"]), len(nc.dimensions["lon"])
    aq_masks, urows, ucols, col_of = build_masks(masks_df, nlat, nlon)
    print(f"[ring] {len(aq_masks)} aquifers; {len(urows)} union cells", flush=True)
    rows = []
    for centre in args.centres:
        cube = args.cube_dir / f"grace_ewh_global_{centre}_appendix.nc"
        print(f"[{centre}] streaming {cube.name} ...", flush=True); ct0 = _time.time()
        slope, _ = slopes_for_centre(cube, urows, ucols)              # (n_recipes, nC)
        print(f"[{centre}] slopes computed ({(_time.time()-ct0)/60:.1f} min)", flush=True)
        for a in aq_masks:
            ro, co = a["cells"]["original"]; rb, cb = a["cells"]["buffered"]
            orig = set(zip(ro.tolist(), co.tolist()))
            ring = [(int(r), int(c)) for r, c in zip(rb.tolist(), cb.tolist()) if (r, c) not in orig]
            cols_o = np.array([col_of[(int(r), int(c))] for r, c in orig])
            cols_r = np.array([col_of[rc] for rc in ring])
            with np.errstate(invalid="ignore"):
                t_o = np.nanmean(slope[:, cols_o], axis=1)
                t_r = np.nanmean(slope[:, cols_r], axis=1) if len(cols_r) else np.full(slope.shape[0], np.nan)
            for ridx in range(slope.shape[0]):
                rows.append((centre, a["aquifer_idx"], a["Study_area"], ridx, len(cols_o), len(cols_r),
                             t_o[ridx], t_r[ridx], t_r[ridx] - t_o[ridx]))
    out = pd.DataFrame(rows, columns=["centre", "aquifer_idx", "Study_area", "recipe_idx_per_centre",
                                      "n_cells_original", "n_cells_ring", "trend_original", "trend_ring", "ring_contrast"])
    args.out.parent.mkdir(parents=True, exist_ok=True); out.to_csv(args.out, index=False, compression="gzip")
    print(f"\nWrote {args.out} ({len(out)} rows)", flush=True)
    if args.canon_pr.exists():
        canon = pd.read_csv(args.canon_pr); canon["centre"] = canon["centre"].str.lower()
        m = canon.merge(out, on=["centre", "aquifer_idx", "recipe_idx_per_centre"], how="inner")
        d = (m["trend_cm_yr"] - m["trend_original"]).abs()
        print(f"[validation] original vs canonical: n={len(m)}, max|diff|={np.nanmax(d):.3e} cm/yr", flush=True)
    return 0

if __name__ == "__main__":
    sys.exit(main())
