"""Per-recipe stripe-residual index.

Some recipe combinations leave north-south stripes a practitioner would
reject on inspection (e.g. a 300-km Gaussian with no destriping); this
script computes the sign-paradox counts on the subset of recipes that
pass a basic stripe-residual criterion.

Stripes are a property of the damping chain -- the (centre, truncation,
filter) combination -- not of the GIA/low-degree/geocentre axes, which add
smooth large-scale fields. This script therefore computes one stripe index
per (centre, truncation, filter) from the global anomaly fields and lets
every recipe inherit its combination's index.

Index: for three well-separated GRACE-era months (April 2005/2010/2015),
take the global 0.5 deg anomaly field, FFT each latitude row within
+-60 deg, and measure the fraction of zonal spectral power at wavenumbers
m >= 16 (the meridional-stripe band; GRACE stripes are sub-monthly-
correlated errors elongated north-south, i.e. high zonal wavenumber).
The index is the band-mean fraction averaged over the three months.
Other axes are held at their first level (gia = none, C20/C30 original,
geocentre none) -- these alter the field by smooth low-degree patterns
that contribute negligibly to m >= 16 power.

Output: datasets/output/tier_a/recipe_stripe_index.csv
        (centre, truncation, filter, stripe_idx)

Run:  .venv python code/grace_pipeline/tier_a_stripe_index.py  (~2 min)
"""
from __future__ import annotations

from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

from . import config as cfg

REPO = cfg.REPO_ROOT
NPZ_DIR = cfg.REPO_ROOT / "datasets/output/tier_a_windows"
CUBE_DIR = cfg.OUT_DIR / "global"
OUT = REPO / "datasets/output/tier_a/recipe_stripe_index.csv"
CENTRES = ("csr", "jpl", "gfz")
TARGET_MONTHS = ("2005-04", "2010-04", "2015-04")
M_MIN = 16          # stripe band: zonal wavenumbers >= m_min
LAT_BAND = 60.0


def stripe_fraction(field: np.ndarray, lat: np.ndarray) -> float:
    """Fraction of zonal spectral power at m >= M_MIN within +-LAT_BAND."""
    band = np.abs(lat) <= LAT_BAND
    rows = field[band]
    rows = rows - rows.mean(axis=1, keepdims=True)
    spec = np.abs(np.fft.rfft(rows, axis=1)) ** 2
    tot = spec[:, 1:].sum(axis=1)
    hi = spec[:, M_MIN:].sum(axis=1)
    ok = tot > 0
    return float(np.mean(hi[ok] / tot[ok]))


def main() -> None:
    rows = []
    for centre in CENTRES:
        npz = np.load(NPZ_DIR / f"monthly_{centre}.npz", allow_pickle=True)
        keys = np.asarray(npz["recipe_keys"])
        cols = list(npz["recipe_cols"])
        ti, fi = cols.index("truncation"), cols.index("filter")
        with netCDF4.Dataset(
                CUBE_DIR / f"grace_ewh_global_{centre}_appendix.nc") as nc:
            tvar = nc.variables["time"]
            dates = netCDF4.num2date(tvar[:], tvar.units,
                                     getattr(tvar, "calendar", "standard"))
            months = [f"{d.year:04d}-{d.month:02d}" for d in dates]
            t_idx = [months.index(m) for m in TARGET_MONTHS]
            lat = nc.variables["lat"][:].filled(np.nan)
            ewh = nc.variables["ewh_anom"]
            n_tr = ewh.shape[1]
            n_f = ewh.shape[2]
            # recipe enumeration: truncation outer, then filter; stride of
            # one filter step = n_gia * n_c20 * n_c30 * n_gc
            stride = int(np.prod(ewh.shape[3:7]))
            n_per_tr = n_f * stride
            for tr in range(n_tr):
                for f in range(n_f):
                    key = keys[tr * n_per_tr + f * stride]
                    fracs = []
                    for t in t_idx:
                        field = np.asarray(
                            ewh[t, tr, f, 0, 0, 0, 0, :, :], dtype=float)
                        if np.isnan(field).all():
                            continue
                        fracs.append(stripe_fraction(
                            np.nan_to_num(field), lat))
                    rows.append({"centre": centre.upper(),
                                 "truncation": int(key[ti]),
                                 "filter": key[fi],
                                 "stripe_idx": float(np.mean(fracs))})
        print(f"[{centre}] {n_tr * n_f} damping combinations indexed",
              flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False, float_format="%.6g")
    print(f"wrote {OUT}")
    print("\nstripe index by filter (mean over centres x truncations):")
    print(df.groupby("filter")["stripe_idx"].mean()
            .sort_values().to_string(float_format="%.3f"))


if __name__ == "__main__":
    main()
