#!/usr/bin/env python3
"""Exact Hochschild-rigidity certificate for the relative-source algebra."""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path


SECTORS = (0, 1, 2)  # e_minus, e_zero, e_plus


def cochain_basis(degree: int) -> list[tuple[int, ...]]:
    """Basis (output, input_1, ..., input_degree) for Hom(A^n,A)."""

    return list(itertools.product(SECTORS, repeat=degree + 1))


def product(left: int, right: int) -> int | None:
    """Primitive-idempotent multiplication e_i e_j=delta_ij e_i."""

    return left if left == right else None


def hochschild_matrix(degree: int) -> list[list[Fraction]]:
    """Matrix of delta:C^degree -> C^(degree+1) over the rationals."""

    source = cochain_basis(degree)
    target = cochain_basis(degree + 1)
    source_index = {basis: index for index, basis in enumerate(source)}
    matrix = [[Fraction(0) for _ in source] for _ in target]

    for target_row, (output, *inputs) in enumerate(target):
        # a_1 f(a_2,...,a_{n+1})
        left_input = inputs[0]
        if product(left_input, output) == output:
            key = (output, *inputs[1:])
            matrix[target_row][source_index[key]] += 1

        # sum_i (-1)^i f(...,a_i a_{i+1},...)
        for position in range(degree):
            merged = product(inputs[position], inputs[position + 1])
            if merged is None:
                continue
            reduced_inputs = (
                inputs[:position] + [merged] + inputs[position + 2 :]
            )
            key = (output, *reduced_inputs)
            matrix[target_row][source_index[key]] += (-1) ** (position + 1)

        # (-1)^(n+1) f(a_1,...,a_n) a_{n+1}
        right_input = inputs[-1]
        if product(output, right_input) == output:
            key = (output, *inputs[:-1])
            matrix[target_row][source_index[key]] += (-1) ** (degree + 1)

    return matrix


def rref(matrix: list[list[Fraction]]) -> tuple[list[list[Fraction]], list[int]]:
    rows = [row[:] for row in matrix]
    if not rows:
        return rows, []
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(len(rows[0])):
        pivot = next(
            (row for row in range(pivot_row, len(rows)) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = rows[pivot_row][column]
        rows[pivot_row] = [value / scale for value in rows[pivot_row]]
        for row in range(len(rows)):
            if row == pivot_row or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(rows[row], rows[pivot_row])
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return rows, pivot_columns


def matrix_rank(matrix: list[list[Fraction]]) -> int:
    return len(rref(matrix)[1])


def nullspace(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    reduced, pivots = rref(matrix)
    column_count = len(matrix[0]) if matrix else 0
    free_columns = [column for column in range(column_count) if column not in pivots]
    vectors: list[list[Fraction]] = []
    for free in free_columns:
        vector = [Fraction(0) for _ in range(column_count)]
        vector[free] = 1
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row][free]
        vectors.append(vector)
    return vectors


def matvec(matrix: list[list[Fraction]], vector: list[Fraction]) -> list[Fraction]:
    return [sum(value * weight for value, weight in zip(row, vector)) for row in matrix]


def matmul(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    if not left or not right:
        return []
    right_columns = list(zip(*right))
    return [
        [sum(a * b for a, b in zip(row, column)) for column in right_columns]
        for row in left
    ]


def contracting_homotopy_2(cocycle: list[Fraction]) -> list[Fraction]:
    """h(phi)(a)=sum_i e_i phi(e_i,a) in the primitive basis."""

    c1 = cochain_basis(1)
    c2_index = {basis: index for index, basis in enumerate(cochain_basis(2))}
    return [
        cocycle[c2_index[(output, output, input_sector)]]
        for output, input_sector in c1
    ]


def build_oracle() -> dict:
    delta1 = hochschild_matrix(1)
    delta2 = hochschild_matrix(2)
    rank_delta1 = matrix_rank(delta1)
    rank_delta2 = matrix_rank(delta2)
    cocycles = nullspace(delta2)
    composition = matmul(delta2, delta1)
    contracted = [
        matvec(delta1, contracting_homotopy_2(cocycle)) == cocycle
        for cocycle in cocycles
    ]

    c1_dimension = len(cochain_basis(1))
    c2_dimension = len(cochain_basis(2))
    c3_dimension = len(cochain_basis(3))
    z2_dimension = c2_dimension - rank_delta2
    b2_dimension = rank_delta1

    return {
        "schema": "matching-one.p233-relative-source-hochschild-rigidity.v1",
        "issues": [54, 114, 233, 252],
        "algebra": {
            "presentation": "A=Q[q]/(q^3-q)",
            "crt_isomorphism": "A is isomorphic to Q x Q x Q by evaluation at q=-1,0,1",
            "primitive_basis": ["e_minus", "e_zero", "e_plus"],
            "multiplication": "e_i e_j=delta_ij e_i",
            "unit": "e_minus+e_zero+e_plus",
            "separability_idempotent": "E=sum_i e_i tensor e_i",
            "separable": True,
        },
        "exact_complex": {
            "dimensions": {"C1": c1_dimension, "C2": c2_dimension, "C3": c3_dimension},
            "rank_delta1": rank_delta1,
            "rank_delta2": rank_delta2,
            "dimension_Z2": z2_dimension,
            "dimension_B2": b2_dimension,
            "dimension_HH2": z2_dimension - b2_dimension,
            "delta2_delta1_is_zero": all(
                value == 0 for row in composition for value in row
            ),
            "nullspace_basis_size": len(cocycles),
            "every_Z2_basis_vector_contracts_to_a_coboundary": all(contracted),
        },
        "explicit_contraction": {
            "definition": "for a 2-cocycle phi, h(phi)(a)=sum_i e_i phi(e_i,a)",
            "identity": "delta h(phi)=phi",
            "meaning": "every first-order associative multiplication deformation is removed by a first-order basis change",
        },
        "research_consequence": {
            "exact_no_go": "The three-sector relative-source algebra itself carries no nontrivial Hochschild first-order deformation class: HH^2(A,A)=0.",
            "stronger_standard_fact": "Because A is finite separable in characteristic zero, HH^n(A,M)=0 for every n>0 and every A-bimodule M.",
            "issue_233": "A nontrivial matching tangent cannot live only in the multiplication of the scalar sector algebra generated by q. It must involve a larger non-semisimple defect/partition algebra, a singular Q-dependent projector, module-extension data, or the doubled-space intertwiner Ad_J.",
            "compatible_existing_result": "The matching-odd D in commit 5c9d8d0 remains a real pull-through/intertwiner obstruction; this certificate says it is not a nontrivial deformation of A=span{1,q,q^2} itself.",
            "next_exact_target": "Compute the derivative class only after adjoining connectivity-resolving junctions/defects; first quotient out the explicit separable three-sector subalgebra, then test the remaining non-semisimple radical/module block.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(build_oracle(), indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
