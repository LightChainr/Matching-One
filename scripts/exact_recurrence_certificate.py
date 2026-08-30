#!/usr/bin/env python3
"""Compile the first uniquely supported exact recurrence of a finite sequence."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "results/model-certificates/synthetic/m2d-vs-m2j/latest.json"
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/finite-recurrence/latest.json"
SCHEMA = "matching-one/exact-finite-recurrence-certificate/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def solve_linear_system(rows: Sequence[Sequence[Fraction]], variable_count: int) -> tuple[str, list[Fraction] | None, int]:
    _require(variable_count >= 1, "variable count must be positive")
    _require(rows and all(len(row) == variable_count + 1 for row in rows), "augmented row shape drift")
    work = [list(row) for row in rows]
    pivot_columns = []
    pivot_row = 0
    for column in range(variable_count):
        selected = next((row for row in range(pivot_row, len(work)) if work[row][column]), None)
        if selected is None:
            continue
        work[pivot_row], work[selected] = work[selected], work[pivot_row]
        pivot = work[pivot_row][column]
        work[pivot_row] = [value / pivot for value in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [left - scale * right for left, right in zip(work[row], work[pivot_row])]
        pivot_columns.append(column)
        pivot_row += 1
    for row in work:
        if not any(row[:variable_count]) and row[variable_count]:
            return "inconsistent", None, len(pivot_columns)
    if len(pivot_columns) < variable_count:
        return "underdetermined", None, len(pivot_columns)
    solution = [Fraction() for _ in range(variable_count)]
    for row, column in enumerate(pivot_columns):
        solution[column] = work[row][variable_count]
    return "supported_unique", solution, len(pivot_columns)


def recurrence_candidate(sequence: Sequence[Fraction], order: int) -> dict[str, Any]:
    _require(order >= 1, "recurrence order must be positive")
    _require(len(sequence) >= 2 * order, "insufficient sequence length for a unique recurrence gate")
    equations = [
        list(sequence[start : start + order]) + [sequence[start + order]]
        for start in range(len(sequence) - order)
    ]
    status, coefficients, coefficient_rank = solve_linear_system(equations, order)
    residuals: list[Fraction] = []
    if coefficients is not None:
        residuals = [
            sequence[start + order]
            - sum((coefficients[index] * sequence[start + index] for index in range(order)), Fraction())
            for start in range(len(sequence) - order)
        ]
        _require(not any(residuals), "supported recurrence has a nonzero residual")
    return {
        "order": order,
        "equation_count": len(equations),
        "coefficient_rank": coefficient_rank,
        "status": status,
        "coefficients_low_to_high": None if coefficients is None else [str(value) for value in coefficients],
        "residuals": [str(value) for value in residuals],
    }


def compile_first_recurrence(sequence: Sequence[Fraction], maximum_order: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    _require(maximum_order >= 1, "maximum order must be positive")
    attempts = []
    selected = None
    for order in range(1, maximum_order + 1):
        candidate = recurrence_candidate(sequence, order)
        attempts.append(candidate)
        if candidate["status"] == "supported_unique":
            selected = candidate
            break
    _require(selected is not None, "no uniquely supported recurrence within the declared order cap")
    return attempts, selected


def build_result(source_path: Path = DEFAULT_SOURCE) -> dict[str, Any]:
    source = json.loads(source_path.read_text(encoding="utf-8"))
    sequence = [Fraction(value) for value in source["synthetic_input"]["moments"]]
    attempts, selected = compile_first_recurrence(sequence, maximum_order=2)
    coefficients = [Fraction(value) for value in selected["coefficients_low_to_high"]]
    characteristic = [-value for value in coefficients] + [Fraction(1)]
    discriminant = characteristic[1] ** 2 - 4 * characteristic[2] * characteristic[0]
    _require(characteristic == [Fraction(1), Fraction(-2), Fraction(1)], "characteristic polynomial drift")
    _require(discriminant == 0, "frozen repeated-root control drift")
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact",
        "source": {
            "path": str(source_path.relative_to(ROOT)),
            "sha256": _sha256_file(source_path),
            "dependency_group": source["synthetic_input"]["dependency_group"],
        },
        "sequence": [str(value) for value in sequence],
        "maximum_order": 2,
        "attempts": attempts,
        "selected_recurrence": selected,
        "characteristic_polynomial_coefficients_low_to_high": [str(value) for value in characteristic],
        "quadratic_discriminant": str(discriminant),
        "summary": {
            "first_uniquely_supported_order": selected["order"],
            "lower_orders_rejected": [row["order"] for row in attempts if row["status"] == "inconsistent"],
            "status": "exact_finite_recurrence_verified",
        },
        "claim_boundary": {
            "included": "the first unique constant-coefficient recurrence supported by every declared finite-sequence equation through order two",
            "excluded": "infinite continuation, flat extension, diagonalizability by itself, noisy recurrence estimation, or physical transfer spectrum",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], source_path: Path = DEFAULT_SOURCE) -> Mapping[str, Any]:
    expected = build_result(source_path)
    _require(result == expected, "finite-recurrence certificate does not exactly reproduce")
    return {
        "schema": result["schema"],
        "status": "valid_exact_finite_recurrence_certificate",
        "selected_order": result["summary"]["first_uniquely_supported_order"],
        "lower_orders_rejected": result["summary"]["lower_orders_rejected"],
        "discriminant": result["quadratic_discriminant"],
        "source_sha256": result["source"]["sha256"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        value = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_result(value, args.source), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_result(args.source), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
