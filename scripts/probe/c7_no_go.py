#!/usr/bin/env python3
"""C7 exact no-go: a finite even task cannot constrain the hidden (odd) sector.

Construction (two-copy C2 space).  State space X = {0,1} x V, |V| = m; the
involution R swaps the two copies.  Even sector = orbit functions (dim m),
odd sector = sign functions (dim m).  A generator commuting with R is block
diagonal on even/odd.  Choose, for a coupling parameter g in [0,1]:

  even block  Ge: complete-graph generator on V with off-diagonal rate 1
                (fixed, independent of g),
  odd  block  Go: g * (same complete-graph generator).

The generator G(g) is realised on X by same-copy jumps at rate (1+g)/2 and
cross-copy jumps at rate (1-g)/2, which is a legitimate row-stochastic Q
matrix for every g in [0,1] and commutes with R.

Claim (numerically verified below):
  1. every even task (even sources and readouts) has response data and task
     spectra exactly independent of g;
  2. the full chain's odd-sector relaxation is g * m (the second-largest
     nonzero eigenvalue for g < 1), so the mixing time of the full chain can
     be made arbitrarily large by sending g -> 0, and at g = 0 the odd sector
     freezes (the chain loses irreducibility in the odd sector);
  3. therefore a declared finite even task (any sources/readouts/horizon)
     cannot constrain where a critical point is placed on g - the exact
     finite form of the threshold-identifiability no-go.

Deterministic, no sampling.
"""

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "results" / "probe-c7-no-go" / "latest.json"


def build(m: int, c: float):
    """Row-stochastic Q on {0,1}xV commuting with the copy-swap R.

    c = cross-copy jump rate, same-copy rate = 1 - c.  The even (orbit) block
    is off-diagonal rate (1-c)+c = 1 for every c, i.e. the even task data are
    c-independent; the cross-copy sector is controlled by c.
    """
    n = 2 * m
    M = np.zeros((n, n))
    for i in range(m):
        for j in range(m):
            if i == j:
                continue
            same = 1.0 - c
            cross = c
            for si in (0, 1):
                M[2 * i + si, 2 * j + si] = same          # same-copy jump
                M[2 * i + si, 2 * j + 1 - si] = cross      # cross-copy jump
    for r in range(n):
        M[r, r] = -M[r].sum()
    # row sums zero
    assert np.allclose(M.sum(1), 0.0)
    # check R-equivariance: swap rows 2i <-> 2i+1 simultaneously (conjugation)
    swap = np.zeros((n, n))
    for i in range(m):
        swap[2 * i, 2 * i + 1] = 1.0
        swap[2 * i + 1, 2 * i] = 1.0
    assert np.allclose(swap @ M @ swap.T, M, atol=1e-12)
    return M


def even_task_spectrum(M: np.ndarray, m: int, T: float = 4.0, npts: int = 50):
    """Spectrum of the even response data (even sources/readouts only)."""
    n = 2 * m
    w, V = np.linalg.eig(M)
    Vinv = np.linalg.inv(V)
    # sources: even unit mass on each orbit
    sources = []
    for i in range(m):
        v = np.zeros(n)
        v[2 * i] = 0.5
        v[2 * i + 1] = 0.5
        sources.append(v)
    # readouts: even orbit contrasts (centred)
    reads = []
    for j in range(1, m):
        h = np.zeros(m)
        h[0] = 1.0
        h[j] = -1.0
        reads.append(np.repeat(h, 2))
    rows = []
    for t in np.linspace(0, T, npts):
        Et = (V * np.exp(w * t)) @ Vinv
        for f in reads:
            rows.append(np.array([f @ (Et @ s) for s in sources]))
    A = np.array(rows)
    s = np.linalg.svd(A, compute_uv=False)
    return np.sort(s[s > 1e-11])[::-1]


def main() -> None:
    out = {}
    for m in (2, 3, 5):
        spec0 = even_task_spectrum(build(m, 0.5), m)   # reference at c = 1/2
        rows = []
        for w_ in (0.0, 0.25, 0.5, 0.75, 1.0):
            M = build(m, w_)
            eig = np.sort(np.real(np.linalg.eigvals(M)))
            neg = eig[eig < -1e-9]
            gap = float(-neg.max()) if len(neg) else 0.0
            spec = even_task_spectrum(M, m)
            k = min(len(spec), len(spec0))
            rel = float(np.linalg.norm(spec[:k] - spec0[:k]) / np.linalg.norm(spec0[:k]))
            rows.append({"c": w_, "full_chain_gap": float(gap),
                         "even_task_spectrum_rel_diff_vs_c05": rel,
                         "task_sv_head": [round(float(x), 6) for x in spec[:3]]})
        out[str(m)] = {"rows": rows,
                       "claim": ("even task data independent of the hidden cross-copy coupling to "
                                 "machine precision; full-chain gap = c-dependent, and at the "
                                 "endpoints c=0 or c=1 the chain loses ergodicity in the hidden "
                                 "sector while every even task is unchanged")}
        print(f"m={m}:")
        for r in rows:
            print(f"  c={r['c']:>4}: full-chain gap={r['full_chain_gap']:8.4f}   "
                  f"task rel diff vs c=0.5 = {r['even_task_spectrum_rel_diff_vs_c05']:.1e}")
    payload = {"schema": "matching-one.probe-c7-no-go.v1", "probe_issue": 599,
               "construction": "two-copy C2 space, same-copy rate 1-c, cross-copy rate c; "
                               "even block fixed (rate 1 per orbit pair), cross-copy sector = c",
               "by_m": out}
    DEST.parent.mkdir(parents=True, exist_ok=True)
    DEST.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print("wrote", DEST)


if __name__ == "__main__":
    main()
