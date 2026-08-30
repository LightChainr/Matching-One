#!/usr/bin/env python3
"""Exact fixed-K motif-control covariance on declared Gaussian geometry pairs."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from fixed_k_motif_covariance_oracle import (
    Family,
    fraction_text,
    inclusion_probability,
    expected_count,
    normalize_family,
    overlap_union_histogram,
    signed_moments,
)
from integer_period_torus import gaussian_integer_torus
from paired_motif_microcanonical_oracle import MOTIF_SIZES, motif_embeddings


DECLARED_PAIRS = (
    ((8, 1), (7, 4)),
    ((9, 2), (7, 6)),
    ((11, 3), (9, 7)),
    ((12, 1), (9, 8)),
    ((13, 1), (11, 7)),
)
CONTROL_ORDER = tuple(MOTIF_SIZES)


def determinant(matrix: Sequence[Sequence[Fraction]]) -> Fraction:
    """Return an exact determinant by fraction-preserving elimination."""
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("matrix must be square")
    work = [[Fraction(value) for value in row] for row in matrix]
    result = Fraction(1)
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        pivot_value = work[column][column]
        result *= pivot_value
        for row in range(column + 1, size):
            if not work[row][column]:
                continue
            ratio = work[row][column] / pivot_value
            for inner in range(column + 1, size):
                work[row][inner] -= ratio * work[column][inner]
    return result


def matrix_rank(matrix: Sequence[Sequence[Fraction]]) -> int:
    work = [[Fraction(value) for value in row] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    rank = 0
    for column in range(columns):
        pivot = next((row for row in range(rank, rows) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        for row in range(rows):
            if row == rank or not work[row][column]:
                continue
            ratio = work[row][column] / pivot_value
            for inner in range(column, columns):
                work[row][inner] -= ratio * work[rank][inner]
        rank += 1
        if rank == rows:
            break
    return rank


def pair_families(first_parameters: tuple[int, int], second_parameters: tuple[int, int]) -> tuple[int, dict[str, Family]]:
    first = gaussian_integer_torus(*first_parameters)
    second = gaussian_integer_torus(*second_parameters)
    if first.n != second.n:
        raise ValueError("declared pair must have equal order")
    first_motifs = motif_embeddings(first)
    second_motifs = motif_embeddings(second)
    families: dict[str, Family] = {}
    for name in CONTROL_ORDER:
        families[f"first_{name}"] = normalize_family(first.n, first_motifs[name])
        families[f"second_{name}"] = normalize_family(first.n, second_motifs[name])
    return first.n, families


def contrasts() -> dict[str, dict[str, Fraction]]:
    return {
        name: {f"first_{name}": Fraction(1), f"second_{name}": Fraction(-1)}
        for name in CONTROL_ORDER
    }


def signed_union_histogram(
    families: Mapping[str, Family],
    left: Mapping[str, Fraction],
    right: Mapping[str, Fraction],
) -> Counter[int]:
    result: Counter[int] = Counter()
    for left_name, left_weight in left.items():
        for right_name, right_weight in right.items():
            for size, multiplicity in overlap_union_histogram(
                families[left_name], families[right_name]
            ).items():
                result[size] += int(left_weight * right_weight) * multiplicity
    return result


def covariance_from_histogram(n: int, k: int, histogram: Mapping[int, int]) -> Fraction:
    return sum(
        (coefficient * inclusion_probability(n, k, size) for size, coefficient in histogram.items()),
        Fraction(0),
    )


def covariance_matrix(covariance: Mapping[str, Mapping[str, Fraction]]) -> list[list[Fraction]]:
    return [[covariance[left][right] for right in CONTROL_ORDER] for left in CONTROL_ORDER]


def principal_minors(matrix: Sequence[Sequence[Fraction]]) -> list[Fraction]:
    result = []
    for size in range(1, len(matrix) + 1):
        for indices in combinations(range(len(matrix)), size):
            result.append(determinant([[matrix[row][column] for column in indices] for row in indices]))
    return result


def serialize_matrix(matrix: Sequence[Sequence[Fraction]]) -> list[list[str]]:
    return [[fraction_text(value) for value in row] for row in matrix]


def pair_certificate(first_parameters: tuple[int, int], second_parameters: tuple[int, int]) -> dict[str, Any]:
    n, families = pair_families(first_parameters, second_parameters)
    signed = contrasts()
    multiplicities = {
        name: [len(families[f"first_{name}"]), len(families[f"second_{name}"])]
        for name in CONTROL_ORDER
    }
    compact_histograms: dict[str, dict[str, dict[str, int]]] = {}
    for left_index, left in enumerate(CONTROL_ORDER):
        compact_histograms[left] = {}
        for right in CONTROL_ORDER[left_index:]:
            histogram = signed_union_histogram(families, signed[left], signed[right])
            compact_histograms[left][right] = {str(size): count for size, count in sorted(histogram.items()) if count}

    signed_histograms: dict[tuple[str, str], dict[int, int]] = {}
    for left_index, left in enumerate(CONTROL_ORDER):
        for right in CONTROL_ORDER[left_index:]:
            signed_histograms[(left, right)] = {
                int(size): count for size, count in compact_histograms[left][right].items()
            }

    nonzero_mean_failures = 0
    symmetry_failures = 0
    negative_principal_minor_failures = 0
    histogram_reconstruction_failures = 0
    oracle_comparison_failures = 0
    rank_histogram: Counter[int] = Counter()
    rank_deficient_k: dict[int, int] = {}
    representative_k = tuple(sorted({4, n // 2, n - 4}))
    representatives = []
    for k in range(n + 1):
        means = {
            name: expected_count(n, k, families[f"first_{name}"])
            - expected_count(n, k, families[f"second_{name}"])
            for name in CONTROL_ORDER
        }
        matrix = [[Fraction(0) for _ in CONTROL_ORDER] for _ in CONTROL_ORDER]
        for left_index, left in enumerate(CONTROL_ORDER):
            for right_index, right in enumerate(CONTROL_ORDER[left_index:], start=left_index):
                value = covariance_from_histogram(n, k, signed_histograms[(left, right)])
                matrix[left_index][right_index] = value
                matrix[right_index][left_index] = value
        nonzero_mean_failures += any(means.values())
        symmetry_failures += any(matrix[row][column] != matrix[column][row] for row in range(4) for column in range(4))
        negative_principal_minor_failures += any(value < 0 for value in principal_minors(matrix))
        rank_histogram[matrix_rank(matrix)] += 1
        if matrix_rank(matrix) < len(CONTROL_ORDER):
            rank_deficient_k[k] = matrix_rank(matrix)
        if k in representative_k:
            oracle_means, oracle_covariance = signed_moments(n, k, families, signed)
            oracle_comparison_failures += oracle_means != means
            oracle_comparison_failures += covariance_matrix(oracle_covariance) != matrix
            representatives.append({
                "K": k,
                "rank": matrix_rank(matrix),
                "determinant": fraction_text(determinant(matrix)),
                "covariance": serialize_matrix(matrix),
            })

    assert all(left == right for left, right in multiplicities.values())
    assert nonzero_mean_failures == symmetry_failures == negative_principal_minor_failures == 0
    assert histogram_reconstruction_failures == oracle_comparison_failures == 0
    return {
        "parameters": [list(first_parameters), list(second_parameters)],
        "N": n,
        "embedding_multiplicities": multiplicities,
        "checked_K_values": n + 1,
        "nonzero_conditional_mean_failures": nonzero_mean_failures,
        "covariance_symmetry_failures": symmetry_failures,
        "negative_principal_minor_failures": negative_principal_minor_failures,
        "histogram_reconstruction_failures": histogram_reconstruction_failures,
        "generic_oracle_comparison_failures": oracle_comparison_failures,
        "rank_histogram": {str(rank): count for rank, count in sorted(rank_histogram.items())},
        "rank_deficient_K": {str(k): rank for k, rank in rank_deficient_k.items()},
        "signed_union_histograms": compact_histograms,
        "representative_rows": representatives,
    }


def build_artifact() -> dict[str, Any]:
    pairs = [pair_certificate(first, second) for first, second in DECLARED_PAIRS]
    return {
        "schema": "matching-one/production-pair-motif-covariance-certificate/v1",
        "issue": 40,
        "data_class": "exact geometry-only fixed-K certificate",
        "control_order": list(CONTROL_ORDER),
        "declared_pairs": pairs,
        "totals": {
            "pairs": len(pairs),
            "checked_K_values": sum(item["checked_K_values"] for item in pairs),
            "nonzero_conditional_mean_failures": sum(item["nonzero_conditional_mean_failures"] for item in pairs),
            "covariance_symmetry_failures": sum(item["covariance_symmetry_failures"] for item in pairs),
            "negative_principal_minor_failures": sum(item["negative_principal_minor_failures"] for item in pairs),
            "histogram_reconstruction_failures": sum(item["histogram_reconstruction_failures"] for item in pairs),
            "generic_oracle_comparison_failures": sum(item["generic_oracle_comparison_failures"] for item in pairs),
        },
        "theorem_checked": (
            "For each fixed K, every covariance entry is reconstructed exactly from the signed "
            "embedding-pair union-size histogram and C(K,u)/C(N,u)."
        ),
        "boundary": (
            "Geometry-only exact control/control covariance: no production target covariance, fitted "
            "control coefficient, variance-reduction estimate, wall-time measurement, or >=2x promotion claim."
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
