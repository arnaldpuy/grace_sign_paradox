"""Phase 4 -- recovery: the unchanged pipeline on the synthetic planet.

Runs `pipeline_global.compute_centre_global` + `tier_a_per_aquifer` per
centre, sequentially (disk budget D006: a synthetic appendix cube is
~111 GB and free space ~150 GB, so each centre's cube is reduced to the
per-aquifer per-recipe CSV and then DELETED before the next centre runs).

Path redirection (D005) -- the ONLY interaction with the pipeline:
  cfg.CENTRE_DIRS       -> synthetic GSM tree (per centre)
  cfg.TN14_FILE         -> synthetic TN-14 (truth + SLR-like noise, D016)
  cfg.GRAVIS_GEOCENTER_FILE -> synthetic GravIS (truth deg-1 + noise, D016)
GIA model files, DDK kernels and Love numbers stay REAL: recipes must
correct the synthetic Earth with the same tools they use on the real one.

CRITICAL import-order note: `io_aux.load_tn14` binds `cfg.TN14_FILE` as a
DEFAULT ARGUMENT, i.e. at import time. The rebinding below therefore
happens after `import config` but BEFORE any other grace_pipeline import.

Run (background, ~5 h per centre + reduction):
    python \
        code/forward/04_recover.py --centres CSR JPL GFZ
"""
from __future__ import annotations

import argparse
import pathlib
import shutil
import subprocess
import sys
import time as _time

PIPELINE_CODE = pathlib.Path(__file__).resolve().parents[3] / "code"
sys.path.insert(0, str(PIPELINE_CODE))

import os

CONTROL_ROOT = pathlib.Path(__file__).resolve().parents[2]
REPO_ROOT = CONTROL_ROOT.parent
# Planet selected via env var (D030) because the cfg path rebinding below must
# happen at import time, before argparse could run.
_PLANET = os.environ.get("CE_PLANET", "pilot")
L2_DIR = CONTROL_ROOT / "datasets" / "synthetic_l2" / _PLANET
AUX_DIR = L2_DIR / "auxiliary"
# Synthetic cubes (~113 GB each) are written to the external cube volume, not
# the internal system disk (D021). Time Machine snapshots only the internal
# disk, so a cube on an external volume is never captured by a local snapshot and deleting
# it reclaims space immediately -- the internal disk stays clear and the
# snapshot-retention trap does not recur. External backups are untouched.
CUBE_DIR = pathlib.Path(os.environ.get("GRACE_CE_CUBE_DIR",
                                       CONTROL_ROOT / "datasets" / "synthetic_l2"))
TIER_A_OUT = CONTROL_ROOT / "datasets" / "output" / (
    "tier_a_synth" if _PLANET == "pilot" else f"tier_a_synth_{_PLANET}")
MASKS_CSV = pathlib.Path(REPO_ROOT / "datasets/jasechko_2024/cell_masks_0p5deg.csv")
PYBIN = sys.executable

MIN_FREE_GB = 125.0

# ---- D005 path redirection (BEFORE the other pipeline imports) -------------
from grace_pipeline import config as cfg  # noqa: E402

cfg.CENTRE_DIRS = {c: L2_DIR / c.lower() for c in cfg.CENTRES}
cfg.TN14_FILE = AUX_DIR / "TN-14_C30_C20_GSFC_SLR.txt"
cfg.GRAVIS_GEOCENTER_FILE = AUX_DIR / "GRAVIS-2B_GFZOP_GEOCENTER_0004.dat"

from grace_pipeline.pipeline_global import compute_centre_global  # noqa: E402


def free_gb(path: pathlib.Path) -> float:
    st = shutil.disk_usage(path)
    return st.free / 1e9


def run_centre(centre: str) -> None:
    t0 = _time.time()
    if not CUBE_DIR.parent.exists():
        raise RuntimeError(
            f"external cube volume {CUBE_DIR.parent} not mounted (D021); "
            f"set GRACE_CE_CUBE_DIR or mount the volume before recovering")
    if free_gb(CUBE_DIR.parent) < MIN_FREE_GB:
        raise RuntimeError(
            f"only {free_gb(CUBE_DIR.parent):.0f} GB free on "
            f"{CUBE_DIR.parent}; need ~{MIN_FREE_GB:.0f} GB for a synthetic "
            f"cube (D006/D021)")
    # sanity: synthetic tree present for this centre
    n = len(list((L2_DIR / centre.lower()).glob("GSM-2_*")))
    if n == 0:
        raise FileNotFoundError(f"no synthetic GSM files for {centre}")
    print(f"[recover {centre}] {n} synthetic GSM files; "
          f"{free_gb(CONTROL_ROOT):.0f} GB free; starting cube", flush=True)

    CUBE_DIR.mkdir(parents=True, exist_ok=True)
    cube = compute_centre_global(centre, window="appendix", out_dir=CUBE_DIR)
    print(f"[recover {centre}] cube done "
          f"({(_time.time()-t0)/3600:.1f} h): {cube}", flush=True)

    out_dir = TIER_A_OUT / centre.lower()
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [PYBIN, "-m", "grace_pipeline.tier_a_per_aquifer",
           "--cube_dir", str(CUBE_DIR),
           "--masks_csv", str(MASKS_CSV),
           "--window", "appendix",
           "--out_dir", str(out_dir),
           "--centres", centre.lower()]
    print(f"[recover {centre}] reducing: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, cwd=PIPELINE_CODE, check=True)

    pr = out_dir / "jasechko_per_recipe_trends.csv"
    n_rows = sum(1 for _ in pr.open()) - 1
    expected = 73 * 640
    if n_rows != expected:
        raise RuntimeError(
            f"{centre}: per-recipe CSV has {n_rows} rows, expected "
            f"{expected}; NOT deleting the cube -- inspect first")

    # D018: extract per-aquifer per-recipe monthly SERIES (~47 MB) before
    # deleting the 111 GB cube -- every later window analysis (rolling null,
    # w_grace/w_fo decomposition, seasonal fits) becomes cheap R work. The
    # parity gate against the canonical tier_a trends must pass before the
    # cube may be deleted.
    series_dir = CONTROL_ROOT / "datasets" / "output" / (
        "aquifer_series" if _PLANET == "pilot"
        else f"aquifer_series_{_PLANET}")
    cmd = [PYBIN, str(CONTROL_ROOT / "code" / "forward" /
                      "04b_extract_series.py"),
           "--cube", str(cube), "--out", str(series_dir),
           "--parity_csv", str(pr)]
    print(f"[recover {centre}] extracting series (D018)", flush=True)
    subprocess.run(cmd, check=True)

    # Hard guard: synthetic cubes share the REAL cubes' filenames (only the
    # directory differs). Refuse to delete anything not inside our own cube
    # output directory -- this protects the real-Earth cubes, which
    # live elsewhere, wherever CUBE_DIR points.
    if CUBE_DIR.resolve() not in cube.resolve().parents:
        raise RuntimeError(
            f"refusing to delete {cube}: not inside {CUBE_DIR}")
    cube.unlink()
    print(f"[recover {centre}] verified {n_rows} rows + series parity; "
          f"cube deleted; total {(_time.time()-t0)/3600:.1f} h", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--centres", nargs="+", default=["CSR", "JPL", "GFZ"])
    args = ap.parse_args()
    for c in args.centres:
        run_centre(c)
    print("[recover] ALL CENTRES DONE", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
