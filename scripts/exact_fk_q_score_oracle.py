#!/usr/bin/env python3
"""Exact critical-manifold Q-score oracle for tiny square FK tori."""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import sys
from typing import Callable, Dict, Iterable, Sequence

from square_bond_kappa3 import BondPair, square_bond_pairs
from torus_homology import HomologyUnionFind


OBSERVABLES = (
    "open_wrap",
    "open_cross",
    "closed_dual_wrap",
    "wrap_difference",
    "open_homology_rank",
)


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _roots(union_find: HomologyUnionFind, vertices: int) -> list[int]:
    return sorted({union_find.find(vertex)[0] for vertex in range(vertices)})


def _rank(union_find: HomologyUnionFind, roots: Iterable[int]) -> int:
    return max((len(union_find.basis[root]) for root in roots), default=0)


def configuration_row(length: int, mask: int, pairs: Sequence[BondPair]) -> Dict[str, int]:
    vertices = length * length
    primal = HomologyUnionFind(vertices, (length, length))
    dual = HomologyUnionFind(vertices, (length, length))
    bonds = 0
    for index, pair in enumerate(pairs):
        if (mask >> index) & 1:
            bonds += 1
            primal.add_edge(*pair.primal)
        else:
            dual.add_edge(*pair.dual)
    primal_roots = _roots(primal, vertices)
    dual_roots = _roots(dual, vertices)
    primal_rank = _rank(primal, primal_roots)
    dual_rank = _rank(dual, dual_roots)
    clusters = len(primal_roots)
    return {
        "b": bonds,
        "k": clusters,
        "J": 2 * clusters + bonds,
        "open_wrap": int(primal_rank > 0),
        "open_cross": int(primal_rank == 2),
        "closed_dual_wrap": int(dual_rank > 0),
        "wrap_difference": int(primal_rank > 0) - int(dual_rank > 0),
        "open_homology_rank": primal_rank,
    }


def mean(rows: Sequence[Dict[str, int]], function: Callable[[Dict[str, int]], Fraction]) -> Fraction:
    return sum((function(row) for row in rows), Fraction(0)) / len(rows)


def ratio_derivatives(rows: Sequence[Dict[str, int]], observable: str, order: int = 3) -> list[Fraction]:
    """Differentiate N(t)/Z(t), with weights exp(t*T), at t=0."""

    z = [mean(rows, lambda row, n=n: Fraction(row["J"], 2) ** n) for n in range(order + 1)]
    numerator = [
        mean(
            rows,
            lambda row, n=n: Fraction(row[observable]) * Fraction(row["J"], 2) ** n,
        )
        for n in range(order + 1)
    ]
    derivatives = [numerator[0]]
    for n in range(1, order + 1):
        derivatives.append(
            numerator[n]
            - sum(comb(n, k) * derivatives[k] * z[n - k] for k in range(n))
        )
    return derivatives


def score_derivatives(rows: Sequence[Dict[str, int]], observable: str) -> list[Fraction]:
    mu_t = mean(rows, lambda row: Fraction(row["J"], 2))
    kappa2 = mean(rows, lambda row: (Fraction(row["J"], 2) - mu_t) ** 2)
    kappa3 = mean(rows, lambda row: (Fraction(row["J"], 2) - mu_t) ** 3)

    def score(row: Dict[str, int], order: int) -> Fraction:
        x = Fraction(row["J"], 2) - mu_t
        if order == 1:
            return x
        if order == 2:
            return x * x - kappa2
        if order == 3:
            return x**3 - 3 * kappa2 * x - kappa3
        raise ValueError(order)

    return [
        mean(rows, lambda row, n=n: Fraction(row[observable]) * score(row, n))
        for n in (1, 2, 3)
    ]


def mixed_eta_t_derivative(rows: Sequence[Dict[str, int]], observable: str) -> Fraction:
    mu_b = mean(rows, lambda row: Fraction(row["b"]))
    mu_t = mean(rows, lambda row: Fraction(row["J"], 2))
    covariance = mean(
        rows,
        lambda row: (Fraction(row["b"]) - mu_b) * (Fraction(row["J"], 2) - mu_t),
    )
    return mean(
        rows,
        lambda row: Fraction(row[observable])
        * (
            (Fraction(row["b"]) - mu_b)
            * (Fraction(row["J"], 2) - mu_t)
            - covariance
        ),
    )


def histogram(rows: Sequence[Dict[str, int]]) -> list[dict]:
    counts: Dict[int, int] = defaultdict(int)
    sums: Dict[str, Dict[int, int]] = {
        observable: defaultdict(int) for observable in OBSERVABLES
    }
    for row in rows:
        j = row["J"]
        counts[j] += 1
        for observable in OBSERVABLES:
            sums[observable][j] += row[observable]
    return [
        {
            "J": j,
            "count": counts[j],
            "observable_sums": {observable: sums[observable][j] for observable in OBSERVABLES},
        }
        for j in sorted(counts)
    ]


def render(length: int = 2) -> dict:
    pairs = square_bond_pairs(length)
    if len(pairs) > 20:
        raise ValueError("exact Q-score oracle is capped at 20 bonds")
    rows = [configuration_row(length, mask, pairs) for mask in range(1 << len(pairs))]
    mu_b = mean(rows, lambda row: Fraction(row["b"]))
    mu_k = mean(rows, lambda row: Fraction(row["k"]))
    mu_t = mean(rows, lambda row: Fraction(row["J"], 2))
    derivatives = {}
    for observable in OBSERVABLES:
        ratio = ratio_derivatives(rows, observable)
        scores = score_derivatives(rows, observable)
        derivatives[observable] = {
            "expectation_at_Q1": fraction_text(ratio[0]),
            "t_derivatives_direct_ratio": [fraction_text(value) for value in ratio[1:]],
            "t_derivatives_score_H1_H2_H3": [fraction_text(value) for value in scores],
            "direct_equals_score": ratio[1:] == scores,
            "mixed_eta_t_derivative": fraction_text(mixed_eta_t_derivative(rows, observable)),
        }
    return {
        "schema": "matching-one.exact-fk-critical-manifold-q-score.v1",
        "issue": 258,
        "geometry": {
            "lattice": "square_bond",
            "L": length,
            "vertices": length * length,
            "bonds": len(pairs),
            "configurations": len(rows),
            "periods": [[length, 0], [0, length]],
        },
        "critical_coordinates": {
            "eta": "log(v/sqrt(Q))",
            "t": "log(Q)",
            "weight": "exp(eta*b+t*T)",
            "T": "k+b/2=J/2",
            "J": "2*k+b",
            "Q1_p": "1/2",
        },
        "exact_means": {
            "b": fraction_text(mu_b),
            "k": fraction_text(mu_k),
            "T": fraction_text(mu_t),
            "T_equals_k_plus_b_over_2": mu_t == mu_k + mu_b / 2,
        },
        "sqrtQ_histogram": histogram(rows),
        "observable_derivatives": derivatives,
        "exact_checks": {
            "histogram_counts_sum_to_configurations": sum(row["count"] for row in histogram(rows)) == len(rows),
            "all_score_orders_match_ratio_differentiation": all(
                row["direct_equals_score"] for row in derivatives.values()
            ),
            "open_closed_wrap_equal_at_Q1": (
                derivatives["open_wrap"]["expectation_at_Q1"]
                == derivatives["closed_dual_wrap"]["expectation_at_Q1"]
            ),
        },
        "interpretation": {
            "exact_result": "One Q=1 ensemble reconstructs the complete finite-volume critical-manifold Q tangent through J=2k+b.",
            "field_warning": "These are derivatives of the FK measure for fixed observables; Q-dependent projectors, normalizations, and insertion definitions remain separate terms.",
            "next_gate": "Add the VJS energy/two-cluster field derivative before calling any T-score direction a logarithmic partner.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=int, default=2)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    payload = render(args.L)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
