#!/usr/bin/env python3
"""Exact polynomial commutant and image-chain valuation for one Jordan block."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence

try:
    from scripts.descendant_jordan_rank_survival import (
        Matrix,
        add,
        apply,
        basis_vector,
        identity,
        jordan_chain_length,
        jordan_nilpotent,
        multiply,
        scale,
    )
except ModuleNotFoundError:  # Direct `python3 scripts/...` execution.
    from descendant_jordan_rank_survival import (  # type: ignore[no-redef]
        Matrix,
        add,
        apply,
        basis_vector,
        identity,
        jordan_chain_length,
        jordan_nilpotent,
        multiply,
        scale,
    )


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "jordan_commutant_valuation_contract.json"


def zero_matrix(rank: int) -> Matrix:
    return [[Fraction() for _ in range(rank)] for _ in range(rank)]


def polynomial_in_nilpotent(rank: int, coefficients: Sequence[Fraction]) -> Matrix:
    if rank < 1 or len(coefficients) > rank:
        raise ValueError("coefficient list must fit a positive Jordan rank")
    nilpotent = jordan_nilpotent(rank)
    power = identity(rank)
    result = zero_matrix(rank)
    for coefficient in coefficients:
        result = add(result, scale(coefficient, power))
        power = multiply(power, nilpotent)
    return result


def coefficient_valuation(coefficients: Sequence[Fraction]) -> Optional[int]:
    return next((index for index, value in enumerate(coefficients) if value), None)


def commutant_coefficients(matrix: Sequence[Sequence[Fraction]]) -> List[Fraction]:
    rank = len(matrix)
    if rank < 1 or any(len(row) != rank for row in matrix):
        raise ValueError("commutant matrix must be nonempty and square")
    nilpotent = jordan_nilpotent(rank)
    if multiply(matrix, nilpotent) != multiply(nilpotent, matrix):
        raise ValueError("matrix does not commute with the Jordan nilpotent")
    coefficients = list(matrix[0])
    if polynomial_in_nilpotent(rank, coefficients) != [list(row) for row in matrix]:
        raise AssertionError("single-block commutant was not upper Toeplitz")
    return coefficients


def image_chain_rank(rank: int, coefficients: Sequence[Fraction]) -> int:
    matrix = polynomial_in_nilpotent(rank, coefficients)
    top_image = apply(matrix, basis_vector(rank, rank - 1))
    return jordan_chain_length(jordan_nilpotent(rank), top_image)


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def matrix_text(matrix: Sequence[Sequence[Fraction]]) -> List[List[str]]:
    return [[fraction_text(value) for value in row] for row in matrix]


def control(rank: int, coefficients: Sequence[Fraction]) -> Dict[str, object]:
    matrix = polynomial_in_nilpotent(rank, coefficients)
    valuation = coefficient_valuation(coefficients)
    expected_rank = 0 if valuation is None else rank - valuation
    actual_rank = image_chain_rank(rank, coefficients)
    if actual_rank != expected_rank:
        raise AssertionError("polynomial valuation did not match image-chain rank")
    return {
        "rank": rank,
        "coefficients": [fraction_text(value) for value in coefficients],
        "valuation": valuation,
        "matrix": matrix_text(matrix),
        "commutant_recovery": [
            fraction_text(value) for value in commutant_coefficients(matrix)
        ],
        "image_chain_rank": actual_rank,
    }


def build_contract() -> Dict[str, object]:
    return {
        "schema": "matching-one/jordan-commutant-valuation/v1",
        "status": "valid_exact_single_block_commutant_certificate",
        "parent_issue": "remain open",
        "commutant_theorem": "A commutes with one rank-r Jordan nilpotent iff A=sum_(k=0)^(r-1) a_k N^k",
        "valuation_rule": "if nu=min{k:a_k!=0}, the image of the top vector has chain rank r-nu; the zero map has rank zero",
        "rank_five_controls": [
            control(5, [Fraction(2), Fraction(-1), Fraction(3), Fraction(), Fraction()]),
            control(5, [Fraction(), Fraction(3), Fraction(-2), Fraction(), Fraction()]),
            control(5, [Fraction(), Fraction(), Fraction(), Fraction(5), Fraction(1)]),
            control(5, [Fraction()] * 5),
        ],
        "thermal_q4_label_control": {
            "rank": 2,
            "coefficients": ["1", "0"],
            "valuation": 0,
            "image_chain_rank": image_chain_rank(2, [Fraction(1), Fraction()]),
            "ordinary_q4_gram_norm": 4930,
        },
        "claim_boundary": {
            "adds_new_virasoro_quotient_calculation": False,
            "proves_lattice_overlap": False,
            "fixes_logarithmic_coefficient": False,
            "derives_torus_ward_response": False,
            "identifies_p4_sprime_readout": False,
        },
    }


def validate_contract(path: Path = CONTRACT) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        frozen = json.load(handle)
    actual = build_contract()
    if frozen != actual:
        raise AssertionError("checked-in Jordan commutant contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
