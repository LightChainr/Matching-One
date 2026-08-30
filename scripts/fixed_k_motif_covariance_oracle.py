#!/usr/bin/env python3
"""Exact fixed-K covariance oracle for signed motif-count combinations."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
from math import comb
from pathlib import Path
from typing import Iterable, Mapping, Sequence


Embedding = frozenset[int]
Family = tuple[Embedding, ...]


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def inclusion_probability(n: int, k: int, size: int) -> Fraction:
    if n < 0 or not 0 <= k <= n:
        raise ValueError("require N >= 0 and 0 <= K <= N")
    if not 0 <= size <= n:
        raise ValueError("union size must lie in [0,N]")
    return Fraction(comb(k, size), comb(n, size)) if k >= size else Fraction(0)


def normalize_family(n: int, embeddings: Iterable[Iterable[int]]) -> Family:
    result = []
    for raw in embeddings:
        vertices = tuple(raw)
        if not vertices or len(set(vertices)) != len(vertices):
            raise ValueError("each embedding must contain distinct vertices")
        if any(not isinstance(vertex, int) or not 0 <= vertex < n for vertex in vertices):
            raise ValueError("embedding vertex lies outside [0,N)")
        result.append(frozenset(vertices))
    return tuple(result)


def overlap_union_histogram(first: Family, second: Family) -> Counter[int]:
    return Counter(len(left | right) for left in first for right in second)


def expected_count(n: int, k: int, family: Family) -> Fraction:
    return sum((inclusion_probability(n, k, len(item)) for item in family), Fraction(0))


def expected_product(n: int, k: int, first: Family, second: Family) -> Fraction:
    histogram = overlap_union_histogram(first, second)
    return sum(
        (multiplicity * inclusion_probability(n, k, size) for size, multiplicity in histogram.items()),
        Fraction(0),
    )


def count_covariance(n: int, k: int, first: Family, second: Family) -> Fraction:
    return expected_product(n, k, first, second) - expected_count(n, k, first) * expected_count(n, k, second)


def signed_moments(
    n: int,
    k: int,
    families: Mapping[str, Family],
    contrasts: Mapping[str, Mapping[str, Fraction]],
) -> tuple[dict[str, Fraction], dict[str, dict[str, Fraction]]]:
    unknown = {
        family for weights in contrasts.values() for family in weights if family not in families
    }
    if unknown:
        raise ValueError("contrast references unknown motif families: " + ", ".join(sorted(unknown)))
    means = {
        name: sum(
            (Fraction(weight) * expected_count(n, k, families[family]) for family, weight in weights.items()),
            Fraction(0),
        )
        for name, weights in contrasts.items()
    }
    covariance: dict[str, dict[str, Fraction]] = {}
    for left_name, left_weights in contrasts.items():
        covariance[left_name] = {}
        for right_name, right_weights in contrasts.items():
            covariance[left_name][right_name] = sum(
                (
                    Fraction(left_weight)
                    * Fraction(right_weight)
                    * count_covariance(n, k, families[left], families[right])
                    for left, left_weight in left_weights.items()
                    for right, right_weight in right_weights.items()
                ),
                Fraction(0),
            )
    return means, covariance


def brute_force_signed_moments(
    n: int,
    k: int,
    families: Mapping[str, Family],
    contrasts: Mapping[str, Mapping[str, Fraction]],
) -> tuple[dict[str, Fraction], dict[str, dict[str, Fraction]]]:
    rows = []
    for occupied_tuple in combinations(range(n), k):
        occupied = frozenset(occupied_tuple)
        counts = {
            name: sum(item <= occupied for item in family) for name, family in families.items()
        }
        rows.append({
            name: sum(
                (Fraction(weight) * counts[family] for family, weight in weights.items()),
                Fraction(0),
            )
            for name, weights in contrasts.items()
        })
    denominator = len(rows)
    means = {
        name: sum((row[name] for row in rows), Fraction(0)) / denominator
        for name in contrasts
    }
    covariance = {
        left: {
            right: sum(
                ((row[left] - means[left]) * (row[right] - means[right]) for row in rows),
                Fraction(0),
            ) / denominator
            for right in contrasts
        }
        for left in contrasts
    }
    return means, covariance


def fixture() -> tuple[int, dict[str, Family], dict[str, dict[str, Fraction]]]:
    n = 6
    families = {
        "ring_edges": normalize_family(n, ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0))),
        "ring_diagonals": normalize_family(n, ((0, 2), (1, 3), (2, 4), (3, 5), (4, 0), (5, 1))),
        "triangles_a": normalize_family(n, ((0, 1, 2), (2, 3, 4), (4, 5, 0))),
        "triangles_b": normalize_family(n, ((1, 2, 3), (3, 4, 5), (5, 0, 1))),
    }
    contrasts = {
        "edge_difference": {"ring_edges": Fraction(1), "ring_diagonals": Fraction(-1)},
        "triangle_difference": {"triangles_a": Fraction(1), "triangles_b": Fraction(-1)},
    }
    return n, families, contrasts


def build_artifact() -> dict:
    n, families, contrasts = fixture()
    rows = []
    failures = 0
    for k in range(n + 1):
        exact = signed_moments(n, k, families, contrasts)
        brute = brute_force_signed_moments(n, k, families, contrasts)
        failures += exact != brute
        means, covariance = exact
        rows.append({
            "K": k,
            "subset_count": comb(n, k),
            "means": {name: fraction_text(value) for name, value in means.items()},
            "covariance": {
                left: {right: fraction_text(value) for right, value in values.items()}
                for left, values in covariance.items()
            },
        })
    return {
        "schema": "matching-one/fixed-k-motif-covariance-oracle/v1",
        "issue": 40,
        "data_class": "synthetic complete K-subset enumeration",
        "N": n,
        "family_sizes": {name: len(family) for name, family in families.items()},
        "contrasts": {
            name: {family: fraction_text(weight) for family, weight in weights.items()}
            for name, weights in contrasts.items()
        },
        "checked_K_values": n + 1,
        "checked_subsets": 1 << n,
        "oracle_vs_enumeration_failures": failures,
        "rows": rows,
        "theorem_checked": "E[1_{A subset X} 1_{B subset X}|K]=binom(K,|A union B|)/binom(N,|A union B|)",
        "boundary": (
            "Synthetic exact covariance only: no production target covariance, fitted control coefficient, "
            "variance-reduction estimate, wall-time measurement, or >=2x promotion claim."
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    rendered = json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
