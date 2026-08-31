#!/usr/bin/env python3
"""Exact spin aliases of the axis/diagonal four-arm landing registry."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Dict, List, Sequence


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "landing_registry_spin_alias_contract.json"


def validate_spin(spin: int) -> None:
    if spin < 0 or spin % 4:
        raise ValueError("the C4-invariant landing registry requires spin 4k")


def cosine_response(spin: int) -> List[int]:
    """Return cos(s theta) at theta=(0, pi/4), exactly for s=4k."""

    validate_spin(spin)
    return [1, -1 if (spin // 4) % 2 else 1]


def sine_response(spin: int) -> List[int]:
    """Return sin(s theta) at theta=(0, pi/4), exactly for s=4k."""

    validate_spin(spin)
    return [0, 0]


def axis_minus_diagonal(spin: int) -> int:
    response = cosine_response(spin)
    return response[0] - response[1]


def matrix_rank(columns: Sequence[Sequence[int]]) -> int:
    if not columns:
        return 0
    height = len(columns[0])
    if any(len(column) != height for column in columns):
        raise ValueError("columns must have equal height")
    rows = [[Fraction(column[row]) for column in columns] for row in range(height)]
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


def build_contract(max_spin: int = 32) -> Dict[str, object]:
    validate_spin(max_spin)
    spins = list(range(0, max_spin + 1, 4))
    response_classes: Dict[str, List[int]] = {
        "axis_equals_diagonal": [],
        "axis_opposes_diagonal": [],
    }
    table = []
    for spin in spins:
        response = cosine_response(spin)
        contrast = axis_minus_diagonal(spin)
        key = "axis_equals_diagonal" if contrast == 0 else "axis_opposes_diagonal"
        response_classes[key].append(spin)
        table.append(
            {
                "spin": spin,
                "cosine_response": response,
                "sine_response": sine_response(spin),
                "axis_minus_diagonal": contrast,
                "alias_class": key,
            }
        )

    h4_h12 = [cosine_response(4), cosine_response(12)]
    return {
        "schema": "matching-one/landing-registry-spin-alias/v1",
        "status": "valid_exact_two_registry_alias_certificate",
        "parent_issue": "remain open",
        "registry_angles": ["0", "pi/4"],
        "observable": "axis minus diagonal",
        "response_formula": "(cos(0), cos(s*pi/4)) = (1, (-1)^(s/4)) for s divisible by 4",
        "selection_rule": {
            "selected": "s congruent to 4 modulo 8",
            "annihilated": "s congruent to 0 modulo 8",
        },
        "low_spin_maximum": max_spin,
        "low_spin_table": table,
        "complete_low_spin_alias_classes": response_classes,
        "h4_h12_no_go": {
            "columns": ["H4", "H12"],
            "response_matrix": [
                [h4_h12[column][row] for column in range(2)] for row in range(2)
            ],
            "rank": matrix_rank(h4_h12),
            "determinant": 0,
            "conclusion": "no linear statistic of these two registry means distinguishes H4 from H12",
        },
        "sine_quadrature_available": False,
        "claim_boundary": {
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
        raise AssertionError("checked-in landing registry alias contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
