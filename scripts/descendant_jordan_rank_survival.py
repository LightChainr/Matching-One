#!/usr/bin/env python3
"""Certify when an intertwining descendant preserves a finite Jordan chain."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Dict, List, Sequence


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "descendant_jordan_rank_survival_contract.json"
Matrix = List[List[Fraction]]
Vector = List[Fraction]


def identity(size: int) -> Matrix:
    return [
        [Fraction(int(row == column)) for column in range(size)]
        for row in range(size)
    ]


def jordan_nilpotent(rank: int) -> Matrix:
    if rank < 1:
        raise ValueError("Jordan rank must be positive")
    return [
        [Fraction(int(column == row + 1)) for column in range(rank)]
        for row in range(rank)
    ]


def add(left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]) -> Matrix:
    if len(left) != len(right) or any(len(a) != len(b) for a, b in zip(left, right)):
        raise ValueError("matrix dimensions must agree")
    return [[a + b for a, b in zip(left_row, right_row)] for left_row, right_row in zip(left, right)]


def scale(factor: Fraction, matrix: Sequence[Sequence[Fraction]]) -> Matrix:
    return [[factor * value for value in row] for row in matrix]


def multiply(left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]) -> Matrix:
    if not left or not right or len(left[0]) != len(right):
        raise ValueError("matrix dimensions do not compose")
    return [
        [
            sum((left[row][inner] * right[inner][column] for inner in range(len(right))), Fraction())
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def apply(matrix: Sequence[Sequence[Fraction]], vector: Sequence[Fraction]) -> Vector:
    if not matrix or len(matrix[0]) != len(vector):
        raise ValueError("matrix and vector dimensions do not agree")
    return [sum((value * vector[index] for index, value in enumerate(row)), Fraction()) for row in matrix]


def dilatation(rank: int, eigenvalue: Fraction) -> Matrix:
    return add(scale(eigenvalue, identity(rank)), jordan_nilpotent(rank))


def descendant_intertwines(
    source: Sequence[Sequence[Fraction]],
    target: Sequence[Sequence[Fraction]],
    descendant: Sequence[Sequence[Fraction]],
    level: int,
) -> bool:
    """Check D_target A = A D_source + level A exactly."""

    return multiply(target, descendant) == add(
        multiply(descendant, source), scale(Fraction(level), descendant)
    )


def basis_vector(size: int, index: int) -> Vector:
    return [Fraction(int(position == index)) for position in range(size)]


def jordan_chain_length(nilpotent: Sequence[Sequence[Fraction]], vector: Sequence[Fraction]) -> int:
    """Return the number of nonzero vectors before repeated N action vanishes."""

    current = list(vector)
    length = 0
    for _ in range(len(nilpotent) + 1):
        if all(value == 0 for value in current):
            return length
        length += 1
        current = apply(nilpotent, current)
    raise AssertionError("matrix advertised as nilpotent did not annihilate the vector")


def bottom_survives(descendant: Sequence[Sequence[Fraction]]) -> bool:
    bottom = basis_vector(len(descendant[0]), 0)
    return any(value != 0 for value in apply(descendant, bottom))


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def matrix_text(matrix: Sequence[Sequence[Fraction]]) -> List[List[str]]:
    return [[fraction_text(value) for value in row] for row in matrix]


def build_contract() -> Dict[str, object]:
    parent_rank = 2
    parent_eigenvalue = Fraction(5, 4)
    level = 4
    parent_nilpotent = jordan_nilpotent(parent_rank)
    source = dilatation(parent_rank, parent_eigenvalue)
    target = dilatation(parent_rank, parent_eigenvalue + level)
    q4 = identity(parent_rank)
    q4_top = apply(q4, basis_vector(parent_rank, parent_rank - 1))

    control_rank = 3
    control_level = 2
    control_nilpotent = jordan_nilpotent(control_rank)
    control_source = dilatation(control_rank, Fraction(7, 3))
    control_target = dilatation(control_rank, Fraction(13, 3))
    collapsing = control_nilpotent
    collapsing_top = apply(collapsing, basis_vector(control_rank, control_rank - 1))

    return {
        "schema": "matching-one/descendant-jordan-rank-survival/v1",
        "status": "valid_exact_finite_jordan_chain_certificate",
        "parent_issue": "remain open",
        "theorem_scope": "finite equal-rank Jordan chains with an exact homogeneous intertwiner",
        "full_rank_condition": "the image of the bottom vector is nonzero",
        "proof_identity": "N_target^(r-1) A top = A N_source^(r-1) top = A bottom",
        "q4_case": {
            "parent_rank": parent_rank,
            "parent_eigenvalue": fraction_text(parent_eigenvalue),
            "descendant_level": level,
            "descendant_eigenvalue": fraction_text(parent_eigenvalue + level),
            "source_dilatation": matrix_text(source),
            "target_dilatation": matrix_text(target),
            "descendant_map": matrix_text(q4),
            "intertwining_identity_holds": descendant_intertwines(source, target, q4, level),
            "bottom_image_nonzero": bottom_survives(q4),
            "image_chain_rank": jordan_chain_length(parent_nilpotent, q4_top),
            "ordinary_q4_gram_norm": 4930,
        },
        "commuting_collapse_control": {
            "parent_rank": control_rank,
            "descendant_level": control_level,
            "descendant_map": "N",
            "intertwining_identity_holds": descendant_intertwines(
                control_source, control_target, collapsing, control_level
            ),
            "bottom_image_nonzero": bottom_survives(collapsing),
            "image_chain_rank": jordan_chain_length(control_nilpotent, collapsing_top),
            "conclusion": "commutation alone does not preserve Jordan rank",
        },
        "claim_boundary": {
            "proves_lattice_overlap": False,
            "fixes_logarithmic_coefficient": False,
            "derives_torus_ward_response": False,
            "identifies_p4_sprime_readout": False,
        },
        "uses_production_data": False,
    }


def validate_contract(path: Path = CONTRACT) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        frozen = json.load(handle)
    actual = build_contract()
    if frozen != actual:
        raise AssertionError("checked-in descendant Jordan contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
