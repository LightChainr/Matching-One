#!/usr/bin/env python3
"""Exact third-registry separation of H4 and H12 landing responses."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Dict, List, Sequence


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "landing_registry_spin_separation_contract.json"


def cosine_alias_at_rational_angle(
    left_spin: int, right_spin: int, numerator: int, denominator: int
) -> bool:
    """Check cos(left*pi*n/d)=cos(right*pi*n/d) using modular integers."""

    if denominator <= 0:
        raise ValueError("angle denominator must be positive")
    left = left_spin * numerator
    right = right_spin * numerator
    modulus = 2 * denominator
    return (left - right) % modulus == 0 or (left + right) % modulus == 0


def matrix_rank(columns: Sequence[Sequence[Fraction]]) -> int:
    if not columns:
        return 0
    height = len(columns[0])
    if any(len(column) != height for column in columns):
        raise ValueError("columns must have equal height")
    rows = [[column[row] for column in columns] for row in range(height)]
    rank = 0
    for column in range(len(columns)):
        pivot = next((row for row in range(rank, height) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        divisor = rows[rank][column]
        rows[rank] = [value / divisor for value in rows[rank]]
        for row in range(height):
            if row == rank or rows[row][column] == 0:
                continue
            factor = rows[row][column]
            rows[row] = [value - factor * base for value, base in zip(rows[row], rows[rank])]
        rank += 1
        if rank == height:
            break
    return rank


def response_columns(include_third_registry: bool) -> List[List[Fraction]]:
    """Return H4/H12 cosine columns at 0, pi/4, and optionally pi/12."""

    h4 = [Fraction(1), Fraction(-1)]
    h12 = [Fraction(1), Fraction(-1)]
    if include_third_registry:
        h4.append(Fraction(1, 2))
        h12.append(Fraction(-1))
    return [h4, h12]


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def matrix_text(columns: Sequence[Sequence[Fraction]]) -> List[List[str]]:
    return [
        [fraction_text(columns[column][row]) for column in range(len(columns))]
        for row in range(len(columns[0]))
    ]


def build_contract() -> Dict[str, object]:
    original = response_columns(False)
    extended = response_columns(True)
    alias_grid = [
        numerator
        for numerator in range(13)
        if cosine_alias_at_rational_angle(4, 12, numerator, 24)
    ]
    return {
        "schema": "matching-one/landing-registry-spin-separation/v1",
        "status": "valid_exact_third_registry_separation_certificate",
        "parent_issue": "remain open",
        "difference_identity": "cos(4 theta)-cos(12 theta)=2 sin(8 theta) sin(4 theta)",
        "alias_criterion": "theta is an integer multiple of pi/8",
        "pi_over_24_grid_alias_numerators_on_0_to_pi_over_2": alias_grid,
        "original_axis_diagonal": {
            "angles": ["0", "pi/4"],
            "columns": ["H4", "H12"],
            "response_matrix": matrix_text(original),
            "rank": matrix_rank(original),
        },
        "extended_registry": {
            "angles": ["0", "pi/4", "pi/12"],
            "columns": ["H4", "H12"],
            "response_matrix": matrix_text(extended),
            "rank": matrix_rank(extended),
            "separating_minor": "-3/2",
            "conclusion": "the pi/12 cosine sample separates H4 and H12 exactly",
        },
        "claim_boundary": {
            "proves_lattice_registry_feasibility": False,
            "runs_new_landing_simulation": False,
            "derives_arm_exponent": False,
            "derives_ope_coefficient": False,
            "identifies_x_21_over_4": False,
        },
    }


def validate_contract(path: Path = CONTRACT) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        frozen = json.load(handle)
    actual = build_contract()
    if frozen != actual:
        raise AssertionError("checked-in landing registry separation contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
