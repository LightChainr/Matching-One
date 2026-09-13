#!/usr/bin/env python3
"""Validate the Akhunzhanov-Eserkepov-Tarasevich (2022, arXiv:2204.01517)
exact torus wrapping polynomials against an independent exact enumeration, and
against the repository's committed continuum baseline.

The ancillary file ``torus.txt`` (https://arxiv.org/src/2204.01517v1/anc/torus.txt)
lists, for each L = 3..12, the integers c_0..c_{L^2} where c_k is the number of
L x L torus site configurations with exactly k occupied sites that wind once
around the horizontal (period-1) cycle.  The wrapping probability polynomial is

    P_L(p) = sum_{k=0}^{L^2} c_k * p^k * (1-p)^{L^2-k}.

At p = 1/2 this is (sum c_k) / 2^{L^2}, which must equal the independent exact
enumeration in ``scripts/exact_wrapping_enum.cpp``.  As L -> infinity, P_L(1/2)
must converge to the repository's committed continuum baseline
``pi_i({1,0})`` = 0.16941543532134688938... (Pinson/Arguin, tau = i).
"""

from __future__ import annotations
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

import mpmath as mp


def parse_torus_txt(path: str) -> dict[int, list[int]]:
    coeffs: dict[int, list[int]] = {}
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("L="):
                current = int(line[2:])
                coeffs[current] = []
            elif current is not None:
                coeffs[current].append(int(line))
    return coeffs


def wrap_probability(coeffs: list[int], p: Fraction) -> Fraction:
    L2 = len(coeffs) - 1
    total = Fraction(0)
    one_minus_p = 1 - p
    for k, c in enumerate(coeffs):
        if c == 0:
            continue
        total += c * p**k * one_minus_p**(L2 - k)
    return total


def main() -> int:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/torus.txt")
    coeffs = parse_torus_txt(str(src))

    # independent exact enumeration counts (scripts/exact_wrapping_enum.cpp)
    enum = {2: 7, 3: 175, 4: 19571, 5: 8853291}

    print("=== Part 1: Akhunzhanov torus wrapping polynomials vs independent enumeration ===")
    print(f"{'L':>3} {'coeff_blocks':>12} {'sum c_k':>12} {'enum count':>12} {'match':>6} {'P_L(1/2) exact':>22}")
    for L in sorted(coeffs):
        c = coeffs[L]
        total = sum(c)
        e = enum.get(L)
        match = (e is not None and total == e)
        pl_half = Fraction(total, 2**(L*L)) if e is not None else None
        print(f"{L:>3} {len(c):>12} {total:>12} {str(e):>12} "
              f"{('YES' if match else '-'):>6} "
              f"{(str(pl_half) if pl_half is not None else '-'):>22}")

    # evaluate at several p with high precision; compare to continuum baseline
    baseline = mp.mpf("0.16941543532134688938260796919875445000145337645375")
    print()
    print("=== P_L(p) for selected p, and convergence to continuum baseline ===")
    print(f"continuum baseline pi_i({{1,0}}) = {mp.nstr(baseline, 22)}")
    ps = [Fraction(1, 2), Fraction(11, 20), Fraction(592746, 1000000)]
    for L in sorted(coeffs):
        c = coeffs[L]
        row = [f"L={L:>2}"]
        for p in ps:
            val = wrap_probability(c, p)
            if p.denominator <= 100000:
                row.append(f"p={float(p):.3f}:{float(val):.10f}")
            else:
                row.append(f"p~0.593:{float(val):.10f}")
        print("  " + "  ".join(row))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
