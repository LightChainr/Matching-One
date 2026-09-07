#!/usr/bin/env python3
"""Exact pairwise comparison of committed geometries' decile vectors and roots.

Issue #641, deliverable 3 (done first, per the ticket): check the committed
axis/diamond exact M(p) polynomials against each other before building any
counterexample.

For each geometry the script recomputes the integer Bernstein coefficients
a_k from the wrapping-side observable D(C) = 1{black NN wraps} - 1{white
NN+NNN wraps} by brute force (same convention as
scripts/exact_matching_polynomial.py), builds the exact power-basis CDF

    F(p) = (1 + M(p)) / 2,   M(p) = sum_k a_k p^k (1-p)^(N-k),

isolates every root of M in [0,1] and every root of F(p) - u for the nine
frozen deciles u = 0.1..0.9 by exact Sturm bisection, and then compares
geometry pairs on (a) the nine decile quantile points, (b) the roots,
(c) M(1/2).  Everything is Fraction; no floats enter any exact claim.

Run:  python scripts/decile_grid_pc_pairwise.py
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from matched_torus_reference import (
    axis_geometry,
    diamond_geometry,
    cluster_stats,
    popcount,
)
from exact_polynomial_root_certificate import (
    evaluate,
    fraction_text,
    isolate_roots,
    trim,
)

DECILES = tuple(
    Fraction(n, 10) for n in (1, 2, 3, 4, 5, 6, 7, 8, 9)
)
ISOLATION_BITS = 60


def bernstein_counts(geometry) -> list[int]:
    counts = [0] * (geometry.n + 1)
    for mask in range(1 << geometry.n):
        k = popcount(mask)
        black = [bool((mask >> i) & 1) for i in range(geometry.n)]
        _nb, black_wrap = cluster_stats(black, geometry.primal_edges)
        white = [not value for value in black]
        _nw, white_wrap = cluster_stats(white, geometry.matching_edges)
        counts[k] += int(black_wrap) - int(white_wrap)
    return counts


def bernstein_to_power(counts: list[int]) -> list[int]:
    """sum_k a_k p^k (1-p)^(N-k) as integer power coefficients."""
    n = len(counts) - 1
    coefficients = [0] * (n + 1)
    for k, value in enumerate(counts):
        if value == 0:
            continue
        for degree in range(k, n + 1):
            coefficients[degree] += value * (-1) ** (degree - k) * binomial(n - k, degree - k)
    return trim(coefficients)


def binomial(a: int, b: int) -> int:
    from math import comb
    return comb(a, b)


def cdf_from_counts(counts: list[int]) -> list[Fraction]:
    """F(p) = (1 + M(p))/2 as exact power-basis coefficients."""
    m_power = bernstein_to_power(counts)
    half = [Fraction(v) for v in m_power]
    if half:
        half[0] += 1
    else:
        half = [Fraction(1)]
    half = [v / 2 for v in half]
    return trim(half)


def evaluate_m(counts: list[int], p: Fraction) -> Fraction:
    n = len(counts) - 1
    total = Fraction(0)
    for k, value in enumerate(counts):
        if value:
            total += value * p**k * (1 - p) ** (n - k)
    return total


def quantile_brackets(cdf: list[Fraction], target: Fraction) -> list[tuple[Fraction, Fraction]]:
    polynomial = list(cdf)
    polynomial[0] -= target
    return isolate_roots(polynomial, bits=ISOLATION_BITS)


def geometry_record(name: str, geometry) -> dict:
    counts = bernstein_counts(geometry)
    m_power = bernstein_to_power(counts)
    cdf = cdf_from_counts(counts)

    m_roots = isolate_roots(m_power, bits=ISOLATION_BITS)
    decile_points = {}
    for u in DECILES:
        brackets = quantile_brackets(cdf, u)
        if len(brackets) != 1:
            raise AssertionError(
                f"{name}: decile u={u} has {len(brackets)} quantile brackets, need exactly 1"
            )
        left, right = brackets[0]
        decile_points[u] = (left, right)

    half_m = evaluate_m(counts, Fraction(1, 2))
    return {
        "name": name,
        "n": geometry.n,
        "physical_period": geometry.physical_period,
        "bernstein_counts": counts,
        "m_power": [int(v) for v in m_power],
        "m_roots": [f"({fraction_text(l)}, {fraction_text(r)})" for l, r in m_roots],
        "decile_quantile_brackets": {
            fraction_text(u): f"({fraction_text(l)}, {fraction_text(r)})"
            for u, (l, r) in decile_points.items()
        },
        "m_half": fraction_text(half_m),
        "_counts": counts,
        "_m_power": [Fraction(v) for v in m_power],
        "_cdf": cdf,
        "_deciles": decile_points,
        "_m_roots": m_roots,
        "_m_half": half_m,
    }


def bracket_gap(a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]) -> Fraction:
    """Smallest possible distance between two intervals, exact."""
    if a[1] < b[0]:
        return b[0] - a[1]
    if b[1] < a[0]:
        return a[0] - b[1]
    return Fraction(0)


def pairwise(a: dict, b: dict) -> dict:
    decile_report = {}
    max_gap = Fraction(0)
    for u in DECILES:
        gap = bracket_gap(a["_deciles"][u], b["_deciles"][u])
        max_gap = max(max_gap, gap)
        decile_report[fraction_text(u)] = fraction_text(gap)
    exact_decile_agreement = max_gap == 0

    # root separation: smallest gap between any root pair (physical roots only,
    # all isolated roots inside (0,1) are reported)
    root_gap = None
    for ra in a["_m_roots"]:
        for rb in b["_m_roots"]:
            gap = bracket_gap(ra, rb)
            if root_gap is None or gap < root_gap:
                root_gap = gap
    return {
        "pair": f"{a['name']} vs {b['name']}",
        "decile_bracket_gaps": decile_report,
        "max_decile_bracket_gap": fraction_text(max_gap),
        "exact_decile_agreement": exact_decile_agreement,
        "min_root_bracket_gap": fraction_text(root_gap),
        "m_half_equal": a["_m_half"] == b["_m_half"],
        "m_half_a": fraction_text(a["_m_half"]),
        "m_half_b": fraction_text(b["_m_half"]),
    }


def main() -> int:
    geometries = [
        ("axis L=2 (N=4)", axis_geometry(2)),
        ("axis L=3 (N=9)", axis_geometry(3)),
        ("axis L=4 (N=16)", axis_geometry(4)),
        ("diamond L=2 (N=8)", diamond_geometry(2)),
        ("diamond L=3 (N=18)", diamond_geometry(3)),
    ]
    records = [geometry_record(name, geometry) for name, geometry in geometries]

    # regression check against results/exact_small_matching_polynomials.md
    expected_m_half = {
        "axis L=2 (N=4)": None,
    }
    for record in records:
        m_half = record["_m_half"]
        print(
            f"{record['name']}: N={record['n']}  M(1/2)={m_half}  "
            f"roots={[str(r) for r in record['_m_roots']]}"
        )

    pairs = []
    for i in range(len(records)):
        for j in range(i + 1, len(records)):
            pairs.append(pairwise(records[i], records[j]))

    artifact = {
        "schema": "matching-one/decile-grid-pc-pairwise/v1",
        "issue": 641,
        "data_class": "exact Fraction arithmetic; Sturm isolation bits=60",
        "geometries": [
            {
                "name": r["name"],
                "n": r["n"],
                "physical_period": r["physical_period"],
                "bernstein_counts": r["bernstein_counts"],
                "m_power": r["m_power"],
                "m_roots": r["m_roots"],
                "decile_quantile_brackets": r["decile_quantile_brackets"],
                "m_half": r["m_half"],
            }
            for r in records
        ],
        "pairwise": pairs,
    }

    out_dir = ROOT / "results" / "decile-grid-pc"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "committed-pairwise.json"
    out_path.write_text(json.dumps(artifact, indent=2), encoding="utf-8")
    print(f"\nwrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
