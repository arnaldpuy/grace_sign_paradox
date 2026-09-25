"""NetCDF writer for the per-centre EWH variant cube, with provenance."""
from __future__ import annotations

import datetime as dt
import subprocess
from pathlib import Path
from typing import Iterable

import gravity_toolkit as gt
import xarray as xr

from . import config as cfg
from .io_shm import Month


def _git_commit() -> str:
    try:
        out = subprocess.check_output(
            ["git", "-C", str(cfg.REPO_ROOT),
             "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
        )
        return out.decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return "unknown"


def write_centre(
    ds: xr.Dataset,
    centre: str,
    months: list[Month],
    ran_filters: Iterable[str],
    out_dir: Path = cfg.OUT_DIR,
    window: str = "poc",
    region: str = "igp",
) -> Path:
    """Write the per-centre Dataset to NetCDF with provenance metadata.

    The output filename composes two orthogonal suffixes:
      * ``region``: 'igp' → bare name (backward-compatible — every existing
        IGP script + the parity harness assume the bare name); 'cv' → '_cv'.
      * ``window``: 'poc' → no suffix; 'appendix' → '_appendix'.
    e.g. ``grace_ewh_variants_csr.nc`` (igp/poc),
    ``grace_ewh_variants_csr_cv.nc`` (cv/poc),
    ``grace_ewh_variants_csr_cv_appendix.nc`` (cv/appendix).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    region_suffix = cfg.REGION_SUFFIX[region]
    window_suffix = "_appendix" if window == "appendix" else ""
    fname = (f"grace_ewh_variants_{centre.lower()}"
             f"{region_suffix}{window_suffix}.nc")
    path = out_dir / fname

    ds = ds.copy()
    ds.attrs.update({
        "window": window,
        "region": region,
        "gravity_toolkit_version": gt.__version__,
        "code_git_commit": _git_commit(),
        "processing_date": dt.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "input_gsm_files": ",".join(m.file_path.name for m in months),
        "tn14_file": cfg.TN14_FILE.name,
        "gia_files": ",".join(p.name for p in cfg.GIA_FILES.values()),
        "love_numbers_source": "PREM (gravity_toolkit/data/love_numbers)",
        "filters_attempted": ",".join(cfg.FILTERS),
        "filters_realised": ",".join(ran_filters),
        "filters_nan_in_output": ",".join(
            f for f in cfg.FILTERS if f not in set(ran_filters)
        ),
        "geocenter_treatments": ",".join(cfg.GEOCENTER_TREATMENTS),
        "geocenter_gravis_file": cfg.GRAVIS_GEOCENTER_FILE.name,
    })

    encoding = {
        v: {"zlib": True, "complevel": 4, "_FillValue": float("nan")}
        for v in ds.data_vars
    }
    ds.to_netcdf(path, encoding=encoding, format="NETCDF4")
    return path
