#!/usr/bin/env python3
"""Exact additive-versus-geometric Hankel calibration for Issue 400."""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import comb
import json
from pathlib import Path
from typing import Any, Mapping, Sequence, Union


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "additive_vs_geometric_hankel_certificate.json"
SCHEMA = "matching-one/additive-vs-geometric-hankel/v1"
ExactInput = Union[int, str, Fraction]
Matrix = tuple[tuple[Fraction, ...], ...]


def exact_fraction(value: ExactInput, *, field: str) -> Fraction:
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError(f"{field} must be exact; floats and booleans are forbidden")
    try:
        return value if isinstance(value, Fraction) else Fraction(value)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"invalid exact value for {field}") from exc


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}" if value.denominator != 1 else str(value.numerator)


def exact_sequence(values: Sequence[ExactInput], *, field: str = "sequence") -> tuple[Fraction, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values:
        raise ValueError(f"{field} must be a nonempty sequence")
    return tuple(exact_fraction(value, field=f"{field}[{index}]") for index, value in enumerate(values))


def hankel(values: Sequence[ExactInput], size: int, *, offset: int = 0) -> Matrix:
    sequence = exact_sequence(values)
    if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
        raise ValueError("size must be a positive integer")
    if isinstance(offset, bool) or not isinstance(offset, int) or offset < 0:
        raise ValueError("offset must be a nonnegative integer")
    required = offset + 2 * size - 1
    if len(sequence) < required:
        raise ValueError(f"sequence needs at least {required} entries")
    return tuple(tuple(sequence[offset + row + column] for column in range(size)) for row in range(size))


def exact_rank(matrix: Sequence[Sequence[ExactInput]]) -> int:
    if isinstance(matrix, (str, bytes)) or not isinstance(matrix, Sequence) or not matrix:
        raise ValueError("matrix must be a nonempty row sequence")
    rows = [list(exact_sequence(row, field=f"matrix[{index}]")) for index, row in enumerate(matrix)]
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("matrix must be rectangular")
    rank = 0
    for column in range(width):
        pivot = next((row for row in range(rank, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][column]
        rows[rank] = [value / pivot_value for value in rows[rank]]
        for row in range(len(rows)):
            if row != rank and rows[row][column]:
                factor = rows[row][column]
                rows[row] = [left - factor * right for left, right in zip(rows[row], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def exact_determinant(matrix: Sequence[Sequence[ExactInput]]) -> Fraction:
    rows = [list(exact_sequence(row, field=f"matrix[{index}]")) for index, row in enumerate(matrix)]
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("determinant requires a nonempty square matrix")
    determinant = Fraction(1)
    for column in range(len(rows)):
        pivot = next((row for row in range(column, len(rows)) if rows[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            determinant = -determinant
        pivot_value = rows[column][column]
        determinant *= pivot_value
        for row in range(column + 1, len(rows)):
            factor = rows[row][column] / pivot_value
            for inner in range(column + 1, len(rows)):
                rows[row][inner] -= factor * rows[column][inner]
    return determinant


def hilbert_determinant(size: int) -> Fraction:
    if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
        raise ValueError("size must be a positive integer")
    denominator = 1
    for index in range(size):
        denominator *= (2 * index + 1) * comb(2 * index, index) ** 2
    return Fraction(1, denominator)


def additive_power_sequence(length: int) -> tuple[Fraction, ...]:
    if isinstance(length, bool) or not isinstance(length, int) or length <= 0:
        raise ValueError("length must be a positive integer")
    return tuple(Fraction(1, index + 1) for index in range(length))


def geometric_power_sequence(length: int, *, ratio: ExactInput = Fraction(1, 2)) -> tuple[Fraction, ...]:
    if isinstance(length, bool) or not isinstance(length, int) or length <= 0:
        raise ValueError("length must be a positive integer")
    q = exact_fraction(ratio, field="ratio")
    if q == 0:
        raise ValueError("ratio must be nonzero")
    return tuple(q**index for index in range(length))


def logarithmic_partner_sequence(length: int, *, ratio: ExactInput = Fraction(1, 2)) -> tuple[Fraction, ...]:
    q = exact_fraction(ratio, field="ratio")
    return tuple(Fraction(index + 1) * q**index for index in range(length))


def exponential_sum_sequence(length: int, ratios: Sequence[ExactInput]) -> tuple[Fraction, ...]:
    qs = exact_sequence(ratios, field="ratios")
    if len(set(qs)) != len(qs) or any(q == 0 for q in qs):
        raise ValueError("ratios must be distinct and nonzero")
    return tuple(sum(q**index for q in qs) for index in range(length))


def ranks_through(values: Sequence[ExactInput], maximum_size: int) -> list[int]:
    return [exact_rank(hankel(values, size)) for size in range(1, maximum_size + 1)]


def build_certificate(maximum_size: int = 8) -> dict[str, Any]:
    if isinstance(maximum_size, bool) or not isinstance(maximum_size, int) or maximum_size < 3:
        raise ValueError("maximum_size must be an integer at least 3")
    length = 2 * maximum_size - 1
    additive = additive_power_sequence(max(length, 33))
    geometric = geometric_power_sequence(length)
    partner = logarithmic_partner_sequence(length)
    three_mode = exponential_sum_sequence(length, (Fraction(1, 2), Fraction(1, 3), Fraction(1, 5)))

    additive_blocks = []
    for size in range(1, maximum_size + 1):
        block = hankel(additive, size)
        determinant = exact_determinant(block)
        formula = hilbert_determinant(size)
        additive_blocks.append({
            "size": size,
            "rank": exact_rank(block),
            "determinant": fraction_text(determinant),
            "formula": fraction_text(formula),
            "formula_matches": determinant == formula,
        })

    recurrence_residuals = [
        (index + 2) * additive[index + 1] - (index + 1) * additive[index]
        for index in range(32)
    ]
    result = {
        "schema": SCHEMA,
        "issue": 400,
        "status": "exact_additive_vs_geometric_hankel_control",
        "maximum_hankel_size": maximum_size,
        "additive_grid": {
            "sequence": "g_n=1/(n+1)",
            "blocks": additive_blocks,
            "ranks": [entry["rank"] for entry in additive_blocks],
            "variable_coefficient_recurrence": "(n+2)g_(n+1)-(n+1)g_n=0",
            "recurrence_indices": [0, 31],
            "recurrence_residuals_all_zero": all(value == 0 for value in recurrence_residuals),
            "constant_recurrence_obstruction": (
                "full size-s Hankel rank for every s=1..8; no nonzero constant-coefficient "
                "recurrence of order <s closes the corresponding checked support window"
            ),
        },
        "geometric_grid": {
            "power": {
                "sequence": "(1/2)^n",
                "ranks": ranks_through(geometric, maximum_size),
                "eventual_rank": 1,
            },
            "logarithmic_partner": {
                "sequence": "(n+1)(1/2)^n",
                "ranks": ranks_through(partner, maximum_size),
                "eventual_rank": 2,
            },
            "three_mode_control": {
                "sequence": "(1/2)^n+(1/3)^n+(1/5)^n",
                "ranks": ranks_through(three_mode, maximum_size),
                "eventual_rank": 3,
            },
        },
        "exact_checks": {
            "all_additive_blocks_full_rank": all(entry["rank"] == entry["size"] for entry in additive_blocks),
            "all_hilbert_determinants_match": all(entry["formula_matches"] for entry in additive_blocks),
            "variable_recurrence_holds": all(value == 0 for value in recurrence_residuals),
            "geometric_power_rank_one": ranks_through(geometric, maximum_size) == [1] * maximum_size,
            "log_partner_rank_two": ranks_through(partner, maximum_size) == [1] + [2] * (maximum_size - 1),
            "three_mode_rank_three": ranks_through(three_mode, maximum_size) == [1, 2] + [3] * (maximum_size - 2),
        },
        "claim_boundary": {
            "included": "finite exact rational Hankel ranks and recurrence controls through size 8",
            "excluded": "P250 production data, physical field counts, covariance, held-out predictions, or continuum interpretation",
            "parent_issue": "remain open",
        },
    }
    validate_certificate(result)
    return result


def validate_certificate(payload: Mapping[str, Any]) -> None:
    if payload.get("schema") != SCHEMA or payload.get("issue") != 400:
        raise ValueError("certificate schema or issue mismatch")
    checks = payload.get("exact_checks")
    if not isinstance(checks, Mapping) or set(checks.values()) != {True}:
        raise ValueError("all exact checks must be true")
    boundary = payload.get("claim_boundary")
    if not isinstance(boundary, Mapping) or boundary.get("parent_issue") != "remain open":
        raise ValueError("parent issue boundary must remain open")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = build_certificate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
