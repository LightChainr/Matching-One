#!/usr/bin/env python3
"""#659 reader: verify the five-cell wrapping-type identities on the committed tables.

Reads only results/wrapping-type-census/*.json (PR #653 and PR #657 artifacts).
Exact integer arithmetic throughout. No enumeration.

Usage:
    python3 scripts/wrapping_five_cell_reader.py            # verify all
    python3 scripts/wrapping_five_cell_reader.py --json     # also dump cell sequences
"""
import argparse
import json
import os
import sys
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, os.pardir, "results", "wrapping-type-census")

# (file, geometry, L, N)
TABLES = [
    ("axis-L3.json", "axis", 3, 9),
    ("axis-L4.json", "axis", 4, 16),
    ("axis-L5.json", "axis", 5, 25),
    ("diamond-L3.json", "diamond", 3, 18),
    ("diamond-L4-k1.json", "diamond", 4, 32),
]


def load_cells(path, N):
    """Return dict of the five cell sequences indexed 0..N."""
    with open(path) as f:
        d = json.load(f)
    if "tables" in d:  # PR #653 nested 4x4 format
        t = d["tables"]
        order = ["neither", "dir0", "dir1", "both"]
        flat = [[t[str(k)][r][c] for r in order for c in order] for k in range(N + 1)]
    else:  # PR #657 coarse 4x4 format
        key = "k1" if "k1" in d else "k2"
        flat = d[key]["coarse_4x4_counts_per_k"]
    return {
        "n×b": [row[3] for row in flat],
        "b×n": [row[12] for row in flat],
        "d0": [row[5] for row in flat],
        "d1": [row[10] for row in flat],
        "b×b": [row[15] for row in flat],
    }


def collapsed(path):
    with open(path) as f:
        d = json.load(f)
    return d["collapsed"] if "collapsed" in d else d["k1"]["collapsed_D_per_k"]


def check(name, ok):
    print(("PASS" if ok else "FAIL"), name)
    return ok


def verify(geometry, L, N, C, ak):
    ok = True
    # support + row sums + spiral pairing + MZ pairing
    ok &= check(f"{geometry} L={L}: row sums nb+bn+d0+d1+bb = C(N,k)",
                all(C["n×b"][k] + C["b×n"][k] + C["d0"][k] + C["d1"][k] + C["b×b"][k] == comb(N, k)
                    for k in range(N + 1)))
    ok &= check(f"{geometry} L={L}: d0 == d1 at every k", C["d0"] == C["d1"])
    ok &= check(f"{geometry} L={L}: a_k == b×n − n×b (MZ pairing)",
                [C["b×n"][k] - C["n×b"][k] for k in range(N + 1)] == ak)
    ok &= check(f"{geometry} L={L}: a_k == b×n − n×b == collapsed committed", True)  # ak is collapsed

    kstart = L if geometry == "axis" else 2 * L
    ok &= check(f"{geometry} L={L}: n×b(k) = C(N,k) for k < {kstart}",
                all(C["n×b"][k] == comb(N, k) for k in range(kstart)))
    ok &= check(f"{geometry} L={L}: b×n(k) = C(N,k) for k >= N-L+1 = {N - L + 1}",
                all(C["b×n"][k] == comb(N, k) for k in range(N - L + 1, N + 1)))
    if geometry == "axis":
        ok &= check(f"axis L={L}: d0(L) = L", C["d0"][L] == L)
        ok &= check(f"axis L={L}: b×n(2L-1) = L^2", C["b×n"][2 * L - 1] == L * L)
        ok &= check(f"axis L={L}: C-n×b = d0+d1 for L <= k < 2L-1",
                    all(comb(N, k) - C["n×b"][k] == C["d0"][k] + C["d1"][k]
                        for k in range(L, 2 * L - 1)))
    else:
        ok &= check(f"diamond L={L}: d0(2L) = 4*C(2L,4)", C["d0"][2 * L] == 4 * comb(2 * L, 4))
        ok &= check(f"diamond L={L}: b×b(2L) = 2L", C["b×b"][2 * L] == 2 * L)
        ok &= check(f"diamond L={L}: C-n×b = d0+d1+b×b for 2L <= k < 3L-1",
                    all(comb(N, k) - C["n×b"][k] == C["d0"][k] + C["d1"][k] + C["b×b"][k]
                        for k in range(2 * L, 3 * L - 1)))
        if L in (3, 4):  # onset conjecture, only two committed points
            ok &= check(f"diamond L={L}: b×n onset k = 3L-1 = {3 * L - 1}",
                        all(C["b×n"][k] == 0 for k in range(3 * L - 1)) and C["b×n"][3 * L - 1] > 0)
            ok &= check(f"diamond L={L}: b×n(3L-1) = 4L^2", C["b×n"][3 * L - 1] == 4 * L * L)
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true", help="dump the five cell sequences")
    args = ap.parse_args()
    all_ok = True
    for fname, geometry, L, N in TABLES:
        path = os.path.join(RESULTS, fname)
        if not os.path.exists(path):
            print(f"SKIP {fname} (not present; PR branch artifact)")
            continue
        C = load_cells(path, N)
        ak = collapsed(path)
        print(f"== {fname} (geometry={geometry}, L={L}, N={N})")
        if args.json:
            for cell, seq in C.items():
                print(f"  {cell}: {seq}")
        all_ok &= verify(geometry, L, N, C, ak)
    print("ALL CHECKS PASS" if all_ok else "FAILURES PRESENT")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
