#!/usr/bin/env python3
"""Path B: brute-force enumeration of all 2^(nm) occupancy patterns.

Completely independent of the transfer matrix in p11_square_dp.py: each
configuration is materialised as a grid, clusters are found with BFS on the
4-neighbour graph, and spanning (a cluster touching rows 1 and m) is decided
directly.  Counts A_{n,m}(k) are tallied by number of occupied sites.

Used to validate the transfer-matrix state canonicalization on small widths
(acceptance requirement of Issue 11).  Feasible for nm <= 16 or so in pure
Python (2^16 = 65536 configurations).
"""

from __future__ import annotations

import argparse
import json
import sys


def brute_force_A(n, m):
    nm = n * m
    A = [0] * (nm + 1)
    for mask in range(1 << nm):
        k = bin(mask).count("1")
        occ = [(mask >> (i * n + j)) & 1 for i in range(m) for j in range(n)]
        # BFS over occupied sites; row index i in 0..m-1 corresponds to
        # paper rows 1..m; spanning touches i=0 and i=m-1.
        seen = [False] * nm
        spans = False
        for start in range(nm):
            if not occ[start] or seen[start]:
                continue
            stack = [start]
            seen[start] = True
            touch_top = touch_bottom = False
            while stack:
                v = stack.pop()
                i, j = divmod(v, n)
                if i == 0:
                    touch_top = True
                if i == m - 1:
                    touch_bottom = True
                for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    ii, jj = i + di, j + dj
                    if 0 <= ii < m and 0 <= jj < n:
                        w = ii * n + jj
                        if occ[w] and not seen[w]:
                            seen[w] = True
                            stack.append(w)
            if touch_top and touch_bottom:
                spans = True
                break
        if spans:
            A[k] += 1
    return A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--m", type=int, default=None)
    args = ap.parse_args()
    m = args.m or args.n
    A = brute_force_A(args.n, m)
    assert A[m] == args.n, f"A[{m}] = {A[m]} != n={args.n}"
    parity = sum(c * (-1) ** k for k, c in enumerate(A) if c)
    print(json.dumps({
        "path": "brute_force",
        "n": args.n, "m": m,
        "A": A,
        "F_minus1": str(parity),
    }))


if __name__ == "__main__":
    main()
