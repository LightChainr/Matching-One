#!/usr/bin/env python3
"""Exact Johnson-slice swap spectrum and polynomial-degree cutoff oracle."""

from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence


def popcount(mask: int) -> int:
    return bin(mask).count("1")


def slice_states(n: int, k: int) -> list[int]:
    if not 0 <= k <= n:
        raise ValueError("slice weight k must lie in [0,N]")
    return [mask for mask in range(1 << n) if popcount(mask) == k]


def swap_matrix(n: int, k: int) -> tuple[list[int], list[list[Fraction]]]:
    states = slice_states(n, k)
    if k in (0, n):
        return states, [[Fraction(1)]]
    index = {mask: row for row, mask in enumerate(states)}
    denominator = k * (n - k)
    matrix = [[Fraction(0) for _ in states] for _ in states]
    for row, mask in enumerate(states):
        occupied = [site for site in range(n) if (mask >> site) & 1]
        empty = [site for site in range(n) if not ((mask >> site) & 1)]
        for remove, add in itertools.product(occupied, empty):
            target = mask ^ (1 << remove) ^ (1 << add)
            matrix[row][index[target]] += Fraction(1, denominator)
    return states, matrix


def matrix_rank(matrix: Sequence[Sequence[Fraction]]) -> int:
    if not matrix:
        return 0
    rows = [list(row) for row in matrix]
    height = len(rows)
    width = len(rows[0])
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (row for row in range(pivot_row, height) if rows[row][column]), None
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [value / scale for value in rows[pivot_row]]
        for row in range(height):
            if row == pivot_row or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(rows[row], rows[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == height:
            break
    return pivot_row


def shifted_nullity(matrix: Sequence[Sequence[Fraction]], eigenvalue: Fraction) -> int:
    shifted = [
        [value - (eigenvalue if row == column else 0) for column, value in enumerate(line)]
        for row, line in enumerate(matrix)
    ]
    return len(matrix) - matrix_rank(shifted)


def subsets_up_to(n: int, degree: int) -> list[int]:
    return [mask for mask in range(1 << n) if popcount(mask) <= degree]


def monomial_evaluation(states: Sequence[int], subsets: Sequence[int]) -> list[list[Fraction]]:
    return [
        [Fraction(1 if state & subset == subset else 0) for subset in subsets]
        for state in states
    ]


def multiply(
    left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]
) -> list[list[Fraction]]:
    return [
        [
            sum(left[row][inner] * right[inner][column] for inner in range(len(right)))
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def concatenate_columns(
    left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]
) -> list[list[Fraction]]:
    return [list(a) + list(b) for a, b in zip(left, right)]


def expected_swap_eigenvalue(n: int, k: int, degree: int) -> Fraction:
    if k in (0, n):
        if degree != 0:
            raise ValueError("endpoint slice has only the constant degree")
        return Fraction(1)
    return Fraction(1) - Fraction(degree * (n - degree + 1), k * (n - k))


def expected_multiplicity(n: int, degree: int) -> int:
    return math.comb(n, degree) - (math.comb(n, degree - 1) if degree else 0)


def build_report(manifest: Mapping[str, object]) -> dict[str, object]:
    n = int(manifest["N"])
    k = int(manifest["k"])
    cutoff = int(manifest["degree_cutoff"])
    states, matrix = swap_matrix(n, k)
    if any(sum(row) != 1 for row in matrix):
        raise AssertionError("swap matrix is not row stochastic")
    spectrum = []
    for degree in range(min(k, n - k) + 1):
        eigenvalue = expected_swap_eigenvalue(n, k, degree)
        multiplicity = expected_multiplicity(n, degree)
        nullity = shifted_nullity(matrix, eigenvalue)
        if nullity != multiplicity:
            raise AssertionError("Johnson eigenvalue multiplicity mismatch")
        spectrum.append(
            {
                "degree": degree,
                "swap_eigenvalue": str(eigenvalue),
                "generator_eigenvalue": str(eigenvalue - 1),
                "expected_multiplicity": multiplicity,
                "exact_nullity": nullity,
            }
        )
    if sum(row["exact_nullity"] for row in spectrum) != len(states):
        raise AssertionError("Johnson eigenspaces do not span the slice")

    degree_rows = []
    for degree in range(cutoff + 1):
        evaluation = monomial_evaluation(states, subsets_up_to(n, degree))
        rank = matrix_rank(evaluation)
        expected = math.comb(n, degree)
        transformed = multiply(matrix, evaluation)
        invariant_rank = matrix_rank(concatenate_columns(evaluation, transformed))
        if rank != expected or invariant_rank != rank:
            raise AssertionError("slice polynomial degree space is not invariant")
        degree_rows.append(
            {
                "degree_at_most": degree,
                "evaluation_rank": rank,
                "expected_rank": expected,
                "rank_after_adjoining_P_images": invariant_rank,
            }
        )
    endpoint_rows = []
    for endpoint in (0, n):
        endpoint_states, endpoint_matrix = swap_matrix(n, endpoint)
        endpoint_rows.append(
            {
                "k": endpoint,
                "states": len(endpoint_states),
                "swap_matrix": [[str(value) for value in row] for row in endpoint_matrix],
                "generator_eigenvalue": "0",
            }
        )
    return {
        "schema": manifest["schema"],
        "status": "exact_johnson_spectrum_and_degree_cutoff_verified",
        "N": n,
        "k": k,
        "slice_states": len(states),
        "swap_denominator": k * (n - k),
        "spectrum": spectrum,
        "degree_spaces": degree_rows,
        "endpoint_slices": endpoint_rows,
        "boundary": manifest["boundary"],
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "analysis/observer_bandwidth_johnson_swap_manifest.json",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_report(json.loads(args.manifest.read_text(encoding="utf-8")))
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
