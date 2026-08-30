#!/usr/bin/env python3
"""Exact two-lift oracle for the Q tangent of the matching projector.

For a periodic FK graph, let ``W_r(Q,v)`` be the restricted state sum whose
open subgraph has ambient homology rank ``r``.  This script compares

    H_Q = W_2 - W_0,
    C_Q = W_2 - Q W_0.

The two lifts have the same Q=1 endpoint but different tangents.  The
difference is algebraic and independent of the finite graph:

    C_Q - H_Q = -(Q-1) W_0.

We nevertheless enumerate the smallest square-bond tori so that the identity,
the normalization term, the fixed-v tangent, and the critical-manifold
``v=sqrt(Q)`` tangent are all recorded in the repository's exact homology
convention.  No fitted or continuum input enters this oracle.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path
from typing import Iterable, Mapping, MutableMapping, Sequence

from square_bond_kappa3 import BondPair, square_bond_pairs
from torus_homology import HomologyUnionFind


Monomial = tuple[int, int]  # powers of (Q, v)
Polynomial = dict[Monomial, int]


def _clean(poly: Mapping[Monomial, int]) -> Polynomial:
    return {key: value for key, value in poly.items() if value}


def add(*polys: Mapping[Monomial, int]) -> Polynomial:
    out: MutableMapping[Monomial, int] = defaultdict(int)
    for poly in polys:
        for monomial, coefficient in poly.items():
            out[monomial] += coefficient
    return _clean(out)


def scale(poly: Mapping[Monomial, int], factor: int) -> Polynomial:
    return _clean({monomial: factor * coefficient for monomial, coefficient in poly.items()})


def multiply_q(poly: Mapping[Monomial, int]) -> Polynomial:
    return {(q_power + 1, v_power): coefficient for (q_power, v_power), coefficient in poly.items()}


def value_at_one(poly: Mapping[Monomial, int]) -> int:
    return sum(poly.values())


def path_derivative_at_one(poly: Mapping[Monomial, int], path: str) -> Fraction:
    """Return d/dQ at Q=1 along fixed v=1 or v=sqrt(Q)."""

    if path == "fixed_v_1":
        return Fraction(sum(coefficient * q_power for (q_power, _), coefficient in poly.items()))
    if path == "critical_square_bond_v_sqrt_Q":
        # Put Q=t^2 and v=t, differentiate in t, then divide by dQ/dt=2t.
        return Fraction(
            sum(
                coefficient * (2 * q_power + v_power)
                for (q_power, v_power), coefficient in poly.items()
            ),
            2,
        )
    raise ValueError(f"unknown path {path!r}")


def normalized_derivative_at_one(
    numerator: Mapping[Monomial, int],
    denominator: Mapping[Monomial, int],
    path: str,
) -> Fraction:
    n0 = value_at_one(numerator)
    z0 = value_at_one(denominator)
    dn = path_derivative_at_one(numerator, path)
    dz = path_derivative_at_one(denominator, path)
    return (dn * z0 - Fraction(n0) * dz) / (z0 * z0)


def _configuration_rank_and_clusters(
    length: int, mask: int, pairs: Sequence[BondPair]
) -> tuple[int, int, int]:
    union_find = HomologyUnionFind(length * length, (length, length))
    occupied = 0
    for index, pair in enumerate(pairs):
        if (mask >> index) & 1:
            occupied += 1
            union_find.add_edge(*pair.primal)
    roots = set()
    max_rank = 0
    for vertex in range(length * length):
        root, _, _ = union_find.find(vertex)
        roots.add(root)
    for root in roots:
        max_rank = max(max_rank, len(union_find.basis[root]))
    return max_rank, len(roots), occupied


def restricted_state_sums(length: int) -> tuple[Polynomial, Polynomial, Polynomial]:
    bond_count = 2 * length * length
    if bond_count > 20:
        raise ValueError("exact two-lift oracle is capped at 20 bonds")
    sums: list[MutableMapping[Monomial, int]] = [defaultdict(int) for _ in range(3)]
    pairs = square_bond_pairs(length)
    for mask in range(1 << bond_count):
        rank, clusters, occupied = _configuration_rank_and_clusters(
            length, mask, pairs
        )
        sums[rank][(clusters, occupied)] += 1
    return tuple(_clean(poly) for poly in sums)  # type: ignore[return-value]


def _fraction(value: Fraction) -> dict[str, int | str]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "fraction": str(value),
    }


def _sparse(poly: Mapping[Monomial, int]) -> list[dict[str, int]]:
    return [
        {"Q_power": q_power, "v_power": v_power, "coefficient": coefficient}
        for (q_power, v_power), coefficient in sorted(poly.items())
    ]


def analyze_length(length: int) -> dict[str, object]:
    w0, w1, w2 = restricted_state_sums(length)
    total = add(w0, w1, w2)
    homology_lift = add(w2, scale(w0, -1))
    critical_polynomial_lift = add(w2, scale(multiply_q(w0), -1))
    lift_difference = add(critical_polynomial_lift, scale(homology_lift, -1))
    expected_difference = add(w0, scale(multiply_q(w0), -1))
    if lift_difference != expected_difference:
        raise AssertionError("C_Q-H_Q != -(Q-1)W_0")

    endpoint_w0 = value_at_one(w0)
    endpoint_total = value_at_one(total)
    paths: dict[str, object] = {}
    for path in ("fixed_v_1", "critical_square_bond_v_sqrt_Q"):
        d_h = path_derivative_at_one(homology_lift, path)
        d_c = path_derivative_at_one(critical_polynomial_lift, path)
        nd_h = normalized_derivative_at_one(homology_lift, total, path)
        nd_c = normalized_derivative_at_one(critical_polynomial_lift, total, path)
        paths[path] = {
            "unnormalized_dH": _fraction(d_h),
            "unnormalized_dC": _fraction(d_c),
            "unnormalized_dC_minus_dH": _fraction(d_c - d_h),
            "expected_unnormalized_counterterm": _fraction(Fraction(-endpoint_w0)),
            "normalized_dh": _fraction(nd_h),
            "normalized_dc": _fraction(nd_c),
            "normalized_dc_minus_dh": _fraction(nd_c - nd_h),
            "expected_normalized_counterterm": _fraction(
                Fraction(-endpoint_w0, endpoint_total)
            ),
            "passed": (
                d_c - d_h == -endpoint_w0
                and nd_c - nd_h == Fraction(-endpoint_w0, endpoint_total)
            ),
        }

    return {
        "L": length,
        "vertices": length * length,
        "bonds": 2 * length * length,
        "configurations": 1 << (2 * length * length),
        "restricted_state_sums": {
            "W_0D": _sparse(w0),
            "W_1D": _sparse(w1),
            "W_2D": _sparse(w2),
        },
        "Q1_v1_sector_counts": {
            "W_0D": endpoint_w0,
            "W_1D": value_at_one(w1),
            "W_2D": value_at_one(w2),
            "Z": endpoint_total,
            "pi_0D": _fraction(Fraction(endpoint_w0, endpoint_total)),
        },
        "endpoint_H_equals_C": value_at_one(homology_lift) == value_at_one(critical_polynomial_lift),
        "exact_polynomial_difference_passed": lift_difference == expected_difference,
        "paths": paths,
        "passed": all(row["passed"] for row in paths.values()),  # type: ignore[union-attr]
    }


def analyze(lengths: Sequence[int]) -> dict[str, object]:
    rows = [analyze_length(length) for length in lengths]
    return {
        "schema": "matching-one/q-lift-covariance-oracle/v1",
        "issue": 333,
        "model": "square_bond_FK_on_square_torus",
        "lifts": {
            "L_hom": "H_Q=W_2D-W_0D",
            "L_CP": "C_Q=W_2D-Q*W_0D",
            "exact_relation": "C_Q-H_Q=-(Q-1)*W_0D",
        },
        "path_rule": (
            "For either declared path, the first tangent difference at Q=1 is "
            "-W_0D before normalization and -pi_0D after normalization. The "
            "individual tangents remain path dependent."
        ),
        "finite_tori": rows,
        "passed": all(row["passed"] for row in rows),
        "scientific_boundary": (
            "This proves a finite counterterm ambiguity between two natural "
            "generic-Q lifts. It does not choose a canonical connection or "
            "identify a logarithmic field."
        ),
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=int, nargs="+", default=[2, 3])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(list(argv) if argv is not None else None)
    result = analyze(args.L)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
