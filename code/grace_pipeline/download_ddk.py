"""Fetch the Bonn DDK decorrelation kernels the pipeline applies.

The kernels (Kusche 2009; Kusche et al. 2009, doi:10.1007/s00190-009-0308-3)
are binary BDFULLV0 files distributed by Roelof Rietbroek's GRACE-filter
repository (MIT licence). The pipeline applies five of the eight levels
(DDK2, DDK3, DDK5, DDK7, DDK8; config.FILTERS). This script downloads them
from a pinned commit of that repository into config.DDK_KERNEL_DIR and
verifies every file against the MD5 checksum of the copy used for the
published results.

Usage
-----
    cd <repository>/code
    python -m grace_pipeline.download_ddk          # fetch what is missing
    python -m grace_pipeline.download_ddk --check  # verify only, no download
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.request

from . import config as cfg
from .ddk import DDK_FILENAME

REPO_URL = "https://github.com/strawpants/GRACE-filter"
COMMIT = "45b991fdbb4407378b362336247e5ca6ac480111"   # data/DDK last changed 2016-08-23
RAW = f"https://raw.githubusercontent.com/strawpants/GRACE-filter/{COMMIT}/data/DDK/"
SIZE_BYTES = 9_793_592

# MD5 of the files used for the published results (identical to the GitHub
# blobs at the pinned commit).
MD5 = {
    "DDK2": "1c7b392ee97ed0c9de122a9157936b6e",   # Wbd_2-120.a_1d13p_4
    "DDK3": "61572deebecb77ae9f8af0efd047839c",   # Wbd_2-120.a_1d12p_4
    "DDK5": "6b91f82649f30467489892ade2a4de32",   # Wbd_2-120.a_1d11p_4
    "DDK7": "37e4066710f17bab5060e58e86da9244",   # Wbd_2-120.a_1d10p_4
    "DDK8": "24b6b90b6368390a7edcffc700d34ec7",   # Wbd_2-120.a_5d9p_4
}


def md5_of(path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true", help="verify only")
    args = ap.parse_args(argv)
    cfg.DDK_KERNEL_DIR.mkdir(parents=True, exist_ok=True)
    bad = 0
    for level, expected in MD5.items():
        name = DDK_FILENAME[level]
        dest = cfg.DDK_KERNEL_DIR / name
        if not dest.exists() and not args.check:
            print(f"[{level}] downloading {name} ...", flush=True)
            urllib.request.urlretrieve(RAW + name, dest)
        if not dest.exists():
            print(f"[{level}] MISSING {dest}")
            bad += 1
            continue
        got = md5_of(dest)
        ok = got == expected and dest.stat().st_size == SIZE_BYTES
        print(f"[{level}] {name}: md5 {got} {'OK' if ok else 'MISMATCH'}")
        bad += 0 if ok else 1
    print(f"{len(MD5) - bad}/{len(MD5)} kernels verified in {cfg.DDK_KERNEL_DIR}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
