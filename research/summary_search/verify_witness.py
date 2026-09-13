#!/usr/bin/env python3
"""Search-independent hard-coded exact-rational check of the n=7 r=1 witness.

Hardcoded incidence lists (no search). Recomputes S(z), r=1 neighborhood
equality, planarity, and the frozen E2_c2 probabilities from scratch using
the same stdlib primitives as bounded_summary_search.py. This is not a
fully independent implementation of the enumerator.
"""
from __future__ import annotations

import sys
from pathlib import Path
from fractions import Fraction
from math import comb

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bounded_summary_search import (
    Net,
    S_coeffs,
    connected_carrier,
    experiments,
    is_plane_two_terminal,
    neighborhood_key,
    safe_table,
)

# Graph A: two 4-paths sharing only corridor ports, plus a triangle bypass.
# Vertices 0..6; 5 = L-port, 6 = R-port.
# L-5-2-1-6-R and L-5-4-0-6-R, plus chord 2-3-1.
A = Net(
    n=7,
    adj_L=1 << 5,
    adj_R=1 << 6,
    adj=(
        (1 << 4) | (1 << 6),  # 0: 4, 6
        (1 << 2) | (1 << 3) | (1 << 6),  # 1: 2, 3, 6
        (1 << 1) | (1 << 3) | (1 << 5),  # 2: 1, 3, 5
        (1 << 1) | (1 << 2),  # 3: 1, 2
        (1 << 0) | (1 << 5),  # 4: 0, 5
        (1 << 2) | (1 << 4),  # 5: 2, 4
        (1 << 0) | (1 << 1),  # 6: 0, 1
    ),
    lr_edge=False,
    name="witnessA",
)

# Graph B: two 4-paths that also share the R-adjacent core vertex 0.
# L-5-2-0-6-R and L-5-4-0-6-R, plus a 5-path 5-2-1-3-6.
B = Net(
    n=7,
    adj_L=1 << 5,
    adj_R=1 << 6,
    adj=(
        (1 << 2) | (1 << 4) | (1 << 6),  # 0: 2, 4, 6
        (1 << 2) | (1 << 3),  # 1: 2, 3
        (1 << 0) | (1 << 1) | (1 << 5),  # 2: 0, 1, 5
        (1 << 1) | (1 << 6),  # 3: 1, 6
        (1 << 0) | (1 << 5),  # 4: 0, 5
        (1 << 2) | (1 << 4),  # 5: 2, 4
        (1 << 0) | (1 << 3),  # 6: 0, 3
    ),
    lr_edge=False,
    name="witnessB",
)


def connecting_ksets(table, n, k):
    out = []
    for mask, safe in enumerate(table):
        if mask.bit_count() == k and not safe:
            verts = tuple(i for i in range(n) if mask & (1 << i))
            out.append(verts)
    return out


def main() -> int:
    for g, tag in ((A, "A"), (B, "B")):
        assert g.n == 7
        assert connected_carrier(g)
        assert is_plane_two_terminal(g), tag
        assert g.adj_L.bit_count() == 1 and g.adj_R.bit_count() == 1

    tA, tB = safe_table(A), safe_table(B)
    sA, sB = S_coeffs(tA, 7), S_coeffs(tB, 7)
    assert sA == sB == (1, 7, 21, 35, 33, 15, 2, 0), (sA, sB)
    assert sA[1] == 7 and sA[2] == comb(7, 2)  # H2=b2=0

    n1A, n1B = neighborhood_key(A, 1), neighborhood_key(B, 1)
    n2A, n2B = neighborhood_key(A, 2), neighborhood_key(B, 2)
    assert n1A == n1B, (n1A, n1B)
    assert n2A != n2B

    eA, eB = experiments(tA, 7), experiments(tB, 7)
    for k in ("E0_c1", "E0_c2", "E0_mix", "E1_c1", "E1_c2", "E1_mix", "E2_c1", "E2_mix"):
        assert eA[k] == eB[k], (k, eA[k], eB[k])
    assert eA["E2_c2"] == Fraction(937, 1050)
    assert eB["E2_c2"] == Fraction(313, 350)
    gap = eA["E2_c2"] - eB["E2_c2"]
    assert gap == Fraction(-1, 525) or gap == Fraction(1, 525)
    assert abs(gap) == Fraction(1, 525)

    cutsA = connecting_ksets(tA, 7, 4)
    cutsB = connecting_ksets(tB, 7, 4)
    assert len(cutsA) == 2 and len(cutsB) == 2
    # A mincuts share only the ports {5,6}; B mincuts share {0,5,6}.
    sa, sb = set(cutsA[0]), set(cutsA[1])
    assert sa & sb == {5, 6}, sa & sb
    ta, tb = set(cutsB[0]), set(cutsB[1])
    assert ta & tb == {0, 5, 6}, ta & tb

    print("VERIFY_OK")
    print("S", sA)
    print("E2_c2", eA["E2_c2"], eB["E2_c2"], "gap", abs(gap))
    print("4-mincuts A", cutsA)
    print("4-mincuts B", cutsB)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
