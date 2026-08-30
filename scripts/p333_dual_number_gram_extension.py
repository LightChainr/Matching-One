#!/usr/bin/env python3
"""Exact first-order Q=1 extension of the connectivity Gram module.

The calculation is over Q[epsilon]/(epsilon^2), with Q=1+epsilon.  It
separates the unit direction from the radical at Q=1 and computes the exact
first-order bilinear form induced on that radical.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import factorial
from pathlib import Path

from p262_confluent_potts_projectors import join_block_count, set_partitions


Matrix = list[list[Fraction]]


def transpose(matrix: Matrix) -> Matrix:
    return [list(row) for row in zip(*matrix)]


def multiply(left: Matrix, right: Matrix) -> Matrix:
    return [
        [sum(left[i][k] * right[k][j] for k in range(len(right)))
         for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def first_jet_gram(n: int) -> tuple[Matrix, Matrix]:
    """Return G0,G1 for G(Q)=G0+epsilon G1 at Q=1+epsilon."""
    partitions = set_partitions(n)
    size = len(partitions)
    g0 = [[Fraction(1) for _ in range(size)] for _ in range(size)]
    g1 = [
        [Fraction(join_block_count(left, right)) for right in partitions]
        for left in partitions
    ]
    return g0, g1


def radical_basis(size: int) -> Matrix:
    """Columns e_i-e_last span the Q=1 Gram radical (coordinate sum zero)."""
    return [
        [Fraction((row == col) - (row == size - 1)) for col in range(size - 1)]
        for row in range(size)
    ]


def restricted_first_form(n: int) -> Matrix:
    _, g1 = first_jet_gram(n)
    basis = radical_basis(len(g1))
    return multiply(transpose(basis), multiply(g1, basis))


def symmetric_inertia(matrix: Matrix) -> tuple[int, int, int]:
    """Exact inertia under rational congruence elimination."""
    rows = [row[:] for row in matrix]
    positive = negative = zero = 0
    while rows:
        size = len(rows)
        diagonal = next((i for i in range(size) if rows[i][i]), None)
        if diagonal is not None:
            if diagonal:
                rows[0], rows[diagonal] = rows[diagonal], rows[0]
                for row in rows:
                    row[0], row[diagonal] = row[diagonal], row[0]
            pivot = rows[0][0]
            positive += pivot > 0
            negative += pivot < 0
            rows = [
                [rows[i][j] - rows[i][0] * rows[0][j] / pivot
                 for j in range(1, size)]
                for i in range(1, size)
            ]
            continue

        off_diagonal = next(
            ((i, j) for i in range(size) for j in range(i + 1, size)
             if rows[i][j]),
            None,
        )
        if off_diagonal is None:
            zero += size
            break

        i, j = off_diagonal
        order = [i, j] + [k for k in range(size) if k not in (i, j)]
        rows = [[rows[a][b] for b in order] for a in order]
        block = [[rows[0][0], rows[0][1]], [rows[1][0], rows[1][1]]]
        determinant = block[0][0] * block[1][1] - block[0][1] * block[1][0]
        inverse = [
            [block[1][1] / determinant, -block[0][1] / determinant],
            [-block[1][0] / determinant, block[0][0] / determinant],
        ]
        positive += 1
        negative += 1
        rows = [
            [
                rows[a][b]
                - sum(rows[a][u] * inverse[u][v] * rows[v][b]
                      for u in range(2) for v in range(2))
                for b in range(2, size)
            ]
            for a in range(2, size)
        ]
    return positive, negative, zero


def stirling_second(n: int, k: int) -> int:
    table = [[0] * (k + 1) for _ in range(n + 1)]
    table[0][0] = 1
    for i in range(1, n + 1):
        for j in range(1, min(i, k) + 1):
            table[i][j] = table[i - 1][j - 1] + j * table[i - 1][j]
    return table[n][k]


def mobius_pivot_jet(k: int) -> int:
    """Derivative at Q=1 of (Q)_k, for k>=2."""
    if k < 2:
        raise ValueError("radical pivot requires k>=2")
    return (-1) ** (k - 2) * factorial(k - 2)


def expected_inertia(n: int) -> tuple[int, int, int]:
    positive = sum(stirling_second(n, k) for k in range(2, n + 1)
                   if mobius_pivot_jet(k) > 0)
    negative = sum(stirling_second(n, k) for k in range(2, n + 1)
                   if mobius_pivot_jet(k) < 0)
    return positive, negative, 0


def sharp_jordan_gate_oracle() -> dict:
    """A two-dimensional exact oracle showing the isotropic-bottom gate is sharp."""
    gram = [[Fraction(0), Fraction(1)], [Fraction(1), Fraction(0)]]
    operator = [[Fraction(2), Fraction(1)], [Fraction(0), Fraction(2)]]
    compatible = multiply(gram, operator) == multiply(transpose(operator), gram)
    bottom = [Fraction(1), Fraction(0)]
    bottom_norm = sum(bottom[i] * gram[i][j] * bottom[j]
                      for i in range(2) for j in range(2))
    return {
        "gram": [[str(value) for value in row] for row in gram],
        "operator": [[str(value) for value in row] for row in operator],
        "gram_self_adjoint": compatible,
        "jordan_eigenvalue": 2,
        "bottom_norm": str(bottom_norm),
    }


def analyze(max_points: int = 5) -> dict:
    rows = []
    for n in range(2, max_points + 1):
        partitions = set_partitions(n)
        observed = symmetric_inertia(restricted_first_form(n))
        expected = expected_inertia(n)
        rows.append({
            "marked_points": n,
            "bell_dimension": len(partitions),
            "endpoint_rank": 1,
            "radical_dimension": len(partitions) - 1,
            "first_radical_form_inertia": {
                "positive": observed[0],
                "negative": observed[1],
                "zero": observed[2],
            },
            "mobius_prediction": {
                "positive": expected[0],
                "negative": expected[1],
                "zero": expected[2],
            },
            "prediction_exact": observed == expected,
            "pivot_multiplicities": {
                str(k): {
                    "multiplicity": stirling_second(n, k),
                    "epsilon_coefficient": mobius_pivot_jet(k),
                }
                for k in range(2, n + 1)
            },
        })
    return {
        "schema_version": 1,
        "issue": 333,
        "ring": "Q[epsilon]/(epsilon^2), Q=1+epsilon",
        "exact_statement": (
            "The Q=1 connectivity Gram module has one unit pivot and, for "
            "each k-block partition with k>=2, one epsilon pivot with leading "
            "coefficient (-1)^(k-2)(k-2)!."
        ),
        "interpretation_boundary": (
            "This is a basis-invariant first-order pairing extension. It does "
            "not by itself identify a transfer-matrix Jordan block or LCFT field."
        ),
        "jordan_selector": {
            "necessary_condition": (
                "For every Gram-compatible regular operator, the endpoint "
                "radical action is self-adjoint for the first-jet form. The "
                "bottom of any nontrivial Jordan chain must therefore be "
                "isotropic in that form."
            ),
            "proof_identity": (
                "<v,v>_H=<v,(T-lambda)w>_H="
                "<(T-lambda)v,w>_H=0"
            ),
            "sharp_two_dimensional_oracle": sharp_jordan_gate_oracle(),
        },
        "checks": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-points", type=int, default=5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.dumps(analyze(args.max_points), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
