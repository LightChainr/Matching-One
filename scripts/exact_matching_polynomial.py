#!/usr/bin/env python3
"""Enumerate exact finite matching polynomials on tiny periodic quotients.

For a geometry with N sites, let

    D(C) = 1{black NN wraps} - 1{white NN+NNN wraps}.

For configurations with exactly k black sites, accumulate the integer
Bernstein coefficient

    a_k = sum_{|C|=k} D(C).

Then the Mertens-Ziff matching function is exactly

    M(p) = sum_{k=0}^N a_k p^k (1-p)^(N-k).

The script also converts this to an integer power-basis polynomial.  This is a
brute-force 2^N reference calculation, intended for exact small-size algebra,
factor/GCD experiments, and regression tests. It is not a frontier algorithm.

Examples:

    python scripts/exact_matching_polynomial.py --geometry axis --L 4
    python scripts/exact_matching_polynomial.py --geometry diamond --L 3 --factor

`--factor` requires SymPy. Numerical roots use mpmath and do not affect the
exact integer coefficients.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import mpmath as mp

from matched_torus_reference import (
    Geometry,
    axis_geometry,
    cluster_stats,
    diamond_geometry,
)


def bernstein_counts(geometry: Geometry) -> list[int]:
    if geometry.n > 26:
        raise ValueError(
            f"N={geometry.n} requires 2^{geometry.n} configurations; "
            "this reference implementation refuses N>26"
        )

    counts = [0] * (geometry.n + 1)
    for mask in range(1 << geometry.n):
        k = mask.bit_count()
        black = [bool((mask >> i) & 1) for i in range(geometry.n)]
        _nb, black_wrap = cluster_stats(black, geometry.primal_edges)

        white = [not value for value in black]
        _nw, white_wrap = cluster_stats(white, geometry.matching_edges)

        counts[k] += int(black_wrap) - int(white_wrap)
    return counts


def bernstein_to_power(counts: list[int]) -> list[int]:
    """Convert sum a_k p^k(1-p)^(N-k) to integer power coefficients."""

    degree_bound = len(counts) - 1
    coefficients = [0] * (degree_bound + 1)
    for k, value in enumerate(counts):
        if value == 0:
            continue
        for degree in range(k, degree_bound + 1):
            coefficients[degree] += (
                value
                * (-1) ** (degree - k)
                * math.comb(degree_bound - k, degree - k)
            )

    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def polynomial_string(coefficients: list[int]) -> str:
    terms: list[str] = []
    for degree in range(len(coefficients) - 1, -1, -1):
        value = coefficients[degree]
        if value == 0:
            continue

        sign = "-" if value < 0 else "+"
        magnitude = abs(value)
        if degree == 0:
            body = str(magnitude)
        elif degree == 1:
            body = "p" if magnitude == 1 else f"{magnitude}*p"
        else:
            body = f"p^{degree}" if magnitude == 1 else f"{magnitude}*p^{degree}"

        if not terms:
            terms.append(body if value > 0 else f"-{body}")
        else:
            terms.append(f" {sign} {body}")
    return "".join(terms) if terms else "0"


def evaluate_power(coefficients: list[int], p: mp.mpf) -> mp.mpf:
    result = mp.mpf(0)
    for value in reversed(coefficients):
        result = result * p + value
    return result


def physical_roots(coefficients: list[int], dps: int) -> list[mp.mpf]:
    mp.mp.dps = dps
    # mp.polyroots wants coefficients highest degree first. Small exact systems
    # are generally well behaved, but the method remains a numerical diagnostic.
    try:
        roots = mp.polyroots(
            [mp.mpf(value) for value in reversed(coefficients)],
            maxsteps=1000,
            error=False,
        )
    except (ValueError, ZeroDivisionError):
        return []

    real = []
    tolerance = mp.power(10, -(dps // 2))
    for root in roots:
        if abs(mp.im(root)) <= tolerance:
            value = mp.re(root)
            if 0 < value < 1:
                real.append(value)
    return sorted(real)


def factor_string(coefficients: list[int]) -> str | None:
    try:
        import sympy as sp
    except ImportError:
        return None

    p = sp.symbols("p")
    expression = sum(
        sp.Integer(value) * p**degree
        for degree, value in enumerate(coefficients)
    )
    return str(sp.factor(expression))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", choices=("axis", "diamond"), required=True)
    parser.add_argument("--L", type=int, required=True)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--factor", action="store_true")
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args()

    if args.dps < 40:
        raise SystemExit("--dps must be at least 40")

    geometry = axis_geometry(args.L) if args.geometry == "axis" else diamond_geometry(args.L)
    try:
        bernstein = bernstein_counts(geometry)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    power = bernstein_to_power(bernstein)
    roots = physical_roots(power, args.dps)
    factor = factor_string(power) if args.factor else None

    print(f"geometry: {geometry.name}")
    print(f"L: {geometry.L}")
    print(f"N: {geometry.n}")
    print(f"physical_period: {geometry.physical_period}")
    print(f"bernstein_counts: {bernstein}")
    print(f"power_coefficients_ascending: {power}")
    print(f"polynomial: {polynomial_string(power)}")
    if roots:
        print("physical_roots:")
        for root in roots:
            print(f"  {mp.nstr(root, 40)}")
    else:
        print("physical_roots: []  # numerical root finder returned none")
    if args.factor:
        print(f"factorization: {factor if factor is not None else 'SymPy not installed'}")

    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "geometry": geometry.name,
            "L": geometry.L,
            "N": geometry.n,
            "physical_period": geometry.physical_period,
            "bernstein_counts": bernstein,
            "power_coefficients_ascending": power,
            "polynomial": polynomial_string(power),
            "physical_roots": [mp.nstr(root, 50) for root in roots],
            "factorization": factor,
        }
        args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
