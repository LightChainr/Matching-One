#!/usr/bin/env python3
"""Exact low-leg Potts connectivity and pair-projector tomography at Q=1."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Iterable


Matrix = list[list[Fraction]]


def set_partitions(n: int) -> list[tuple[tuple[int, ...], ...]]:
    """Return canonical set partitions of range(n)."""
    out: list[tuple[tuple[int, ...], ...]] = []

    def visit(i: int, blocks: list[list[int]]) -> None:
        if i == n:
            out.append(tuple(tuple(block) for block in blocks))
            return
        for j in range(len(blocks)):
            blocks[j].append(i)
            visit(i + 1, blocks)
            blocks[j].pop()
        blocks.append([i])
        visit(i + 1, blocks)
        blocks.pop()

    visit(0, [])
    return out


def join_block_count(left: tuple[tuple[int, ...], ...], right: tuple[tuple[int, ...], ...]) -> int:
    n = sum(map(len, left))
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        a, b = find(a), find(b)
        if a != b:
            parent[b] = a

    for partition in (left, right):
        for block in partition:
            for x in block[1:]:
                union(block[0], x)
    return len({find(x) for x in range(n)})


def connectivity_gram(n: int, q: int) -> Matrix:
    partitions = set_partitions(n)
    return [
        [Fraction(q ** join_block_count(left, right)) for right in partitions]
        for left in partitions
    ]


def matrix_rank(matrix: Matrix) -> int:
    rows = [row[:] for row in matrix]
    if not rows:
        return 0
    rank = 0
    for col in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][col]
        rows[rank] = [value / scale for value in rows[rank]]
        for i in range(len(rows)):
            if i == rank or not rows[i][col]:
                continue
            scale = rows[i][col]
            rows[i] = [a - scale * b for a, b in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def determinant(matrix: Matrix) -> Fraction:
    rows = [row[:] for row in matrix]
    result = Fraction(1)
    for col in range(len(rows)):
        pivot = next((i for i in range(col, len(rows)) if rows[i][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            rows[col], rows[pivot] = rows[pivot], rows[col]
            result *= -1
        value = rows[col][col]
        result *= value
        for i in range(col + 1, len(rows)):
            scale = rows[i][col] / value
            for j in range(col + 1, len(rows)):
                rows[i][j] -= scale * rows[col][j]
    return result


def falling(q: int, k: int) -> int:
    value = 1
    for j in range(k):
        value *= q - j
    return value


def expected_gram_determinant(n: int, q: int) -> int:
    multiplicities: dict[int, int] = {}
    for partition in set_partitions(n):
        k = len(partition)
        multiplicities[k] = multiplicities.get(k, 0) + 1
    value = 1
    for k, count in multiplicities.items():
        value *= falling(q, k) ** count
    return value


def identity(n: int) -> Matrix:
    return [[Fraction(i == j) for j in range(n)] for i in range(n)]


def add(*matrices: Matrix) -> Matrix:
    return [[sum(matrix[i][j] for matrix in matrices) for j in range(len(matrices[0]))] for i in range(len(matrices[0]))]


def scale(matrix: Matrix, scalar: Fraction) -> Matrix:
    return [[scalar * value for value in row] for row in matrix]


def multiply(left: Matrix, right: Matrix) -> Matrix:
    return [
        [sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def trace(matrix: Matrix) -> Fraction:
    return sum(matrix[i][i] for i in range(len(matrix)))


def unordered_pair_projectors(q: int) -> dict[str, Matrix]:
    """Project C[2-subsets] into [], [1], and [2] at integer q >= 4."""
    if q < 4:
        raise ValueError("integer realization requires q >= 4")
    edges = list(combinations(range(q), 2))
    size = len(edges)
    one = [[Fraction(1) for _ in edges] for _ in edges]
    overlap = [[Fraction(len(set(a) & set(b))) for b in edges] for a in edges]
    p0 = scale(one, Fraction(2, q * (q - 1)))
    p1 = scale(add(overlap, scale(one, Fraction(-4, q))), Fraction(1, q - 2))
    p2 = add(identity(size), scale(p0, -1), scale(p1, -1))
    return {"singlet": p0, "standard": p1, "two_row_2": p2}


def projectors_are_orthogonal(projectors: Iterable[Matrix]) -> bool:
    values = list(projectors)
    zero = [[Fraction(0) for _ in values[0]] for _ in values[0]]
    for i, left in enumerate(values):
        if multiply(left, left) != left:
            return False
        for j, right in enumerate(values):
            if i != j and multiply(left, right) != zero:
                return False
    return True


def fraction_record(value: Fraction) -> dict[str, int | str]:
    return {"numerator": value.numerator, "denominator": value.denominator, "text": str(value)}


def analyze() -> dict:
    integer_checks = []
    for q in (4, 5, 6):
        projectors = unordered_pair_projectors(q)
        integer_checks.append({
            "Q": q,
            "edge_space_dimension": q * (q - 1) // 2,
            "orthogonal_idempotents": projectors_are_orthogonal(projectors.values()),
            "traces": {name: fraction_record(trace(matrix)) for name, matrix in projectors.items()},
        })
    gram_checks = []
    for n in (2, 4):
        gram_checks.append({
            "marked_points": n,
            "bell_dimension": len(set_partitions(n)),
            "rank_at_Q1": matrix_rank(connectivity_gram(n, 1)),
            "determinant_factorization": {
                "formula": "product_{k=1}^n (Q)_k^S(n,k)",
                "stirling_multiplicities": {
                    str(k): sum(len(p) == k for p in set_partitions(n)) for k in range(1, n + 1)
                },
            },
            "integer_oracle_Q5": {
                "direct_determinant": str(determinant(connectivity_gram(n, 5))),
                "factorized_determinant": str(expected_gram_determinant(n, 5)),
            },
        })
    return {
        "schema_version": 1,
        "issue": 262,
        "scope": "exact low-leg algebra; no continuum overlap inferred",
        "connectivity_gram": gram_checks,
        "unordered_pair_representation": {
            "meaning": "two distinct FK cluster colours, the minimal four-leg/two-cluster carrier",
            "generic_decomposition": "C[2-subsets]=[] + [1] + [2]",
            "diagram_basis": ["I_pair", "X_shared_colour", "J_all_ones"],
            "projectors": {
                "singlet": "2 J/[Q(Q-1)]",
                "standard_[1]": "(X-4J/Q)/(Q-2)",
                "two_row_[2]": "I-P_singlet-P_[1]",
            },
            "categorical_traces": {
                "singlet": "1",
                "standard_[1]": "Q-1",
                "two_row_[2]": "Q(Q-3)/2",
            },
            "Q_to_1_laurent": {
                "residue_(Q-1)P_singlet": "2 J",
                "residue_(Q-1)P_[2]": "-2 J",
                "P_[1]_at_Q1": "-X+4J",
                "regular_confluent_block_at_Q1": "lim(P_singlet+P_[2])=I+X-4J",
                "regular_block_derivative_at_Q1": "d_Q(P_singlet+P_[2])=X",
                "trace_collision": "tr(P_singlet)+tr(P_[2])=(Q-1)(Q-2)/2",
            },
            "integer_checks": integer_checks,
        },
        "logarithmic_gate": {
            "exact_statement": "projector poles alone do not imply an LCFT logarithm",
            "required_extra_condition": "x_singlet(1)=x_[2](1) with nonzero derivative difference",
            "conditional_log_coefficient": "2 R [x_[2]'(1)-x_singlet'(1)] multiplying log(r) r^(-2x*)",
            "tensor_residue": "R=2J in this normalization",
        },
        "separation_of_derivatives": {
            "measure_score_issue_258": "Cov(P O_bare,T_Q)",
            "projector_geometry_issue_262": "finite collision combination of (d_Q P)O_bare; individual d_QP has a double pole",
            "explicit_bare_operator": "<P d_Q O_bare>",
            "virasoro_gate_issue_252": "the 4:-6:3 Ward row distinguishes thermal Q4 from a spin-4 primary and is not part of Q differentiation",
        },
        "frozen_risk_predictions": [
            "Any singlet/[2] confluent tensor residue in the unordered-pair carrier is proportional to J with exact residue ratio -1.",
            "The [1] projector is coefficientwise regular at Q=1; its trace vanishes linearly, so a pole in a normalized [1] observable is not forced by projector algebra.",
            "A claimed VJS logarithmic coefficient must factor into the exact tensor residue and an independently nonzero scaling-dimension slope difference.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.dumps(analyze(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
