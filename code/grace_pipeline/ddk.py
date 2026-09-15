"""Reader and applier for Bonn DDK (Kusche 2009) decorrelation filter kernels.

Kernels are distributed as binary BDFULLV0 files (Bonn / Frommle convention,
vnum 2.1) at https://github.com/strawpants/GRACE-filter under MIT licence.
The format is reverse-engineered from the file header here; a higher-level
reference (LGPL-licensed) is Roelof Rietbroek's `frommle.io.BINV.readBIN`.
We re-implemented the parse so this project stays under its own licence.

Mathematical structure: a DDK kernel is a per-order block-diagonal linear
operator on Stokes coefficients. For each spherical-harmonic order m it
mixes degrees l, l' within that m. Cosine and sine blocks are independent
for m ≥ 1; m = 0 has only a cosine block.
"""
from __future__ import annotations

import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


# ---------------------------------------------------------------------------
# Mapping from Bonn filename convention to DDK level.
# Regularization parameter alpha encoded as `<n>d<p>p_4`  ⇒  alpha = n × 10^p.
# ---------------------------------------------------------------------------
DDK_FILENAME = {
    "DDK1": "Wbd_2-120.a_1d14p_4",
    "DDK2": "Wbd_2-120.a_1d13p_4",
    "DDK3": "Wbd_2-120.a_1d12p_4",
    "DDK4": "Wbd_2-120.a_5d11p_4",
    "DDK5": "Wbd_2-120.a_1d11p_4",
    "DDK6": "Wbd_2-120.a_5d10p_4",
    "DDK7": "Wbd_2-120.a_1d10p_4",
    "DDK8": "Wbd_2-120.a_5d9p_4",
}


@dataclass(frozen=True)
class DDKKernel:
    """Parsed DDK kernel: per-order (m, trig) → block matrix, plus L bounds."""
    lmax: int
    lmin: int
    # blocks[(m, trig)] is the (lmax-m+1) × (lmax-m+1) weight matrix for that
    # order and trig (0 = cosine, 1 = sine).
    blocks: dict[tuple[int, int], np.ndarray]
    source_file: str


def read_ddk_kernel(path: Path) -> DDKKernel:
    """Parse a Bonn-format BDFULLV0 binary DDK kernel file.

    Tested against vnum 2.1 / type BDFULLV0 / nvec 0 (the only flavour the
    project actually consumes). Other variants of the BIN format are not
    supported here; raise rather than silently misread.
    """
    with open(path, "rb") as fh:
        data = fh.read()

    # Endianness magic (B, I → 18754 little-endian)
    (endian_check,) = struct.unpack_from("<H", data, 0)
    if endian_check != 18754:
        raise ValueError(
            f"{path}: endian magic {endian_check} != 18754; big-endian or "
            f"corrupt — not handled."
        )
    version = "BI" + data[2:8].decode("utf-8")
    vnum = float(version[4:7])
    type_ = data[8:16].decode("utf-8")
    if type_ != "BDFULLV0":
        raise ValueError(f"{path}: type {type_!r} ≠ 'BDFULLV0' — not handled.")
    if not (2.0 <= vnum <= 2.1):
        raise ValueError(f"{path}: vnum {vnum} outside tested range [2.0, 2.1].")

    nints, ndbls, nval1, _nval2 = struct.unpack_from("<IIII", data, 96)
    pval1, _pval2 = struct.unpack_from("<II", data, 112)  # vnum < 2.4 ⇒ 4-byte

    # vnum ≤ 2.1 + BDFULLV0  ⇒  nvec=0, nread=0, nval2=nval1, pval2=1
    cursor = 120
    (nblocks,) = struct.unpack_from("<I", data, cursor)
    cursor += 4

    # Integer metadata: 24-byte names + (vnum ≤ 2.4) 4-byte uints
    int_names = np.frombuffer(data, dtype="|S24", count=nints, offset=cursor)
    cursor += 24 * nints
    int_vals = np.frombuffer(data, dtype="<i4", count=nints, offset=cursor)
    cursor += 4 * nints

    # Double metadata: 24-byte names + 8-byte doubles
    dbl_names = np.frombuffer(data, dtype="|S24", count=ndbls, offset=cursor)
    cursor += 24 * ndbls
    dbl_vals = np.frombuffer(data, dtype="<f8", count=ndbls, offset=cursor)
    cursor += 8 * ndbls

    meta: dict[str, float | int] = {}
    for n, v in zip(int_names, int_vals):
        meta[n.decode("utf-8").strip()] = int(v)
    for n, v in zip(dbl_names, dbl_vals):
        meta[n.decode("utf-8").strip()] = float(v)

    if "Lmax" not in meta or "Lmin" not in meta:
        raise ValueError(f"{path}: kernel header missing Lmax/Lmin in meta {meta}")
    lmax = int(meta["Lmax"])
    lmin = int(meta["Lmin"])

    # side1_d (24-byte coefficient labels), then blockind (4-byte uints)
    cursor += 24 * nval1
    blockind = np.frombuffer(data, dtype="<i4", count=nblocks, offset=cursor).copy()
    cursor += 4 * nblocks
    # vnum > 2.2 has side2_d; we are at vnum 2.1 so side2_d == side1_d implicitly.

    # Packed double matrix: pval1 doubles (one per matrix element across all blocks)
    pack = np.frombuffer(data, dtype="<f8", count=pval1, offset=cursor)
    if pack.size != pval1:
        raise ValueError(
            f"{path}: read {pack.size} doubles, expected pval1={pval1}"
        )

    # Unpack into per-(m, trig) blocks. Block 0 → m=0/cos; blocks 1,2 → m=1
    # cos/sin; blocks 3,4 → m=2 cos/sin; etc. (matches Frommle/Bonn convention.)
    blocks: dict[tuple[int, int], np.ndarray] = {}
    last_ind = 0
    pos = 0
    for iblck in range(nblocks):
        if iblck == 0:
            m, trig = 0, 0
        else:
            m = (iblck + 1) // 2
            trig = (iblck - 1) % 2  # 0 cos, 1 sin
        sz = int(blockind[iblck]) - last_ind
        last_ind = int(blockind[iblck])

        # Identity-padded full block of size (lmax - m + 1). The packed data
        # populates the (lmax - max(lmin, m) + 1) × (...) sub-block aligned at
        # [shift:, shift:]; the leading rows/cols are identity (degrees below
        # max(lmin, m) are passed through unchanged).
        full = np.eye(lmax - m + 1, dtype=np.float64)
        nminblk = max(lmin, m)
        shift = nminblk - m
        # Frommle loads as C order (so np.dot(C_in, block) gives correct output).
        full[shift:, shift:] = pack[pos : pos + sz * sz].reshape((sz, sz), order="C")
        blocks[(m, trig)] = full
        pos += sz * sz

    if pos != pval1:
        raise ValueError(
            f"{path}: consumed {pos} doubles, expected {pval1} (block sum mismatch)"
        )

    return DDKKernel(lmax=lmax, lmin=lmin, blocks=blocks, source_file=path.name)


def apply_ddk(harm: Any, kernel: DDKKernel) -> Any:
    """Apply DDK kernel to a gravity_toolkit.harmonics object in place on a copy.

    Each spherical-harmonic order m has a (lmax_field-m+1) × (lmax_field-m+1)
    sub-block taken from the kernel; clm[m:, m] and slm[m:, m] are multiplied
    (vector @ matrix) by it. m = 0 has cosine only.
    """
    h = harm.copy()
    lmax_field = int(h.lmax)
    if lmax_field > kernel.lmax:
        raise ValueError(
            f"DDK kernel lmax={kernel.lmax} < field lmax={lmax_field}"
        )

    for m in range(lmax_field + 1):
        n = lmax_field - m + 1  # length of the per-m vector in the field

        block_c = kernel.blocks[(m, 0)][:n, :n]
        c_vec = h.clm[m : lmax_field + 1, m]
        h.clm[m : lmax_field + 1, m] = c_vec @ block_c

        if m >= 1:
            block_s = kernel.blocks[(m, 1)][:n, :n]
            s_vec = h.slm[m : lmax_field + 1, m]
            h.slm[m : lmax_field + 1, m] = s_vec @ block_s
        # m=0 has no sine; slm[:, 0] is already 0 in any GSM.

    return h
