#!/usr/bin/env python3
"""Exact matching-parity and identifiability audit for the four P4 channels."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from math import factorial
from pathlib import Path
from typing import Any, Sequence


@dataclass(frozen=True)
class Field:
    key: str
    label: str
    dimension: Fraction
    matching_parity: int


FIELDS = (
    Field("I", "matching-even identity-family spin 4", Fraction(4), 1),
    Field("T", "matching-odd thermal-family spin 4", Fraction(21, 4), -1),
)
THERMAL_EXPONENT = Fraction(3, 4)
PRIMARY_CHANNELS = (("P4_S", "S", 0), ("P4_D", "D", 0),
                    ("P4_S_prime", "S", 1), ("P4_D_prime", "D", 1))


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def selection_coefficient(field: Field, combination: str, derivative_order: int) -> Fraction:
    if combination not in ("S", "D"):
        raise ValueError("combination must be S or D")
    if derivative_order < 0:
        raise ValueError("derivative order must be nonnegative")
    sign = 1 if combination == "S" else -1
    return Fraction(1 + sign * field.matching_parity * ((-1) ** derivative_order), 2)


def length_exponent(field: Field, derivative_order: int) -> Fraction:
    return Fraction(2) - field.dimension + derivative_order * THERMAL_EXPONENT


def matrix_rank(matrix: Sequence[Sequence[Fraction]]) -> int:
    rows = [list(row) for row in matrix]
    if not rows:
        return 0
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("ragged matrix")
    rank = 0
    for column in range(width):
        pivot = next((row for row in range(rank, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][column]
        rows[rank] = [value / scale for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank or not rows[row][column]:
                continue
            scale = rows[row][column]
            rows[row] = [a - scale * b for a, b in zip(rows[row], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def build_artifact() -> dict[str, Any]:
    derivative_table = []
    for field in FIELDS:
        for derivative_order in range(4):
            for combination in ("S", "D"):
                coefficient = selection_coefficient(field, combination, derivative_order)
                exponent_l = length_exponent(field, derivative_order)
                derivative_table.append({
                    "field": field.key,
                    "combination": combination,
                    "derivative_order": derivative_order,
                    "allowed": bool(coefficient),
                    "selection_coefficient": fraction_text(coefficient),
                    "derivative_factor": factorial(derivative_order) if coefficient else 0,
                    "L_exponent": fraction_text(exponent_l),
                    "N_exponent": fraction_text(exponent_l / 2),
                })

    columns = ["I:f0", "I:f1", "T:f0", "T:f1"]
    matrix = []
    primary_rows = []
    for channel, combination, derivative_order in PRIMARY_CHANNELS:
        row = []
        contributions = []
        for column in columns:
            field_key, coefficient_key = column.split(":")
            coefficient_order = int(coefficient_key[1:])
            field = next(item for item in FIELDS if item.key == field_key)
            value = Fraction(0)
            if coefficient_order == derivative_order:
                value = selection_coefficient(field, combination, derivative_order) * factorial(derivative_order)
            row.append(value)
            if value:
                exponent_l = length_exponent(field, derivative_order)
                contributions.append({
                    "column": column,
                    "coefficient": fraction_text(value),
                    "L_exponent": fraction_text(exponent_l),
                    "N_exponent": fraction_text(exponent_l / 2),
                })
        matrix.append(row)
        primary_rows.append({"channel": channel, "entries": [fraction_text(value) for value in row],
                             "contributions": contributions})

    rank = matrix_rank(matrix)
    assert rank == 4
    return {
        "schema": "matching-one/operator-mixing-parity-oracle/v1",
        "issue": 125,
        "thermal_exponent_y_t": fraction_text(THERMAL_EXPONENT),
        "fields": [
            {"key": field.key, "label": field.label,
             "dimension": fraction_text(field.dimension), "matching_parity": field.matching_parity}
            for field in FIELDS
        ],
        "selection_rule": {
            "S": "allowed iff (-1)^n = eta",
            "D": "allowed iff (-1)^n = -eta",
        },
        "derivative_table": derivative_table,
        "primary_map": {
            "rows": [channel for channel, _, _ in PRIMARY_CHANNELS],
            "columns": columns,
            "matrix": [[fraction_text(value) for value in row] for row in matrix],
            "rank": rank,
        },
        "identifiability": {
            "independent_coefficients": 4,
            "observables": 4,
            "conclusion": "the symmetry-allowed four-channel map is full rank",
            "not_implied_by_symmetry": ["f_I1/f_I0", "f_T1/f_T0"],
            "two_amplitude_reduction_requires": (
                "two external relations fixing the derivative-to-center Taylor ratios; "
                "parity alone supplies only structural zeros"
            ),
        },
        "primary_rows": primary_rows,
        "boundary": (
            "This exact audit fixes selection zeros and powers only. It does not fit amplitudes, "
            "identify fields, or supply the missing Taylor-ratio dynamics."
        ),
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    lines = [
        "# Four-channel operator-mixing parity audit", "",
        "The matching involution fixes an exact zero pattern, but it does not reduce the four", 
        "allowed Taylor coefficients to two amplitudes.", "",
        "| channel | nonzero coefficient | L power | N power |", "|---|---|---:|---:|",
    ]
    for row in artifact["primary_rows"]:
        contribution = row["contributions"][0]
        lines.append("| %s | `%s` | `%s` | `%s` |" % (
            row["channel"], contribution["column"], contribution["L_exponent"],
            contribution["N_exponent"]))
    lines.extend([
        "", "## Exact identifiability conclusion", "",
        "With columns `[I:f0, I:f1, T:f0, T:f1]`, the structural matrix has rank `4`.",
        "Therefore parity supplies zeros and channel assignments, not the ratios `f_I1/f_I0` or",
        "`f_T1/f_T0`. A genuine two-amplitude joint prediction needs those two relations from an",
        "independent dynamical calculation or frozen training data.", "", "## Boundary", "",
        artifact["boundary"], "",
    ])
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = build_artifact()
    rendered = (json.dumps(artifact, indent=2, sort_keys=True) + "\n"
                if args.format == "json" else render_markdown(artifact))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
