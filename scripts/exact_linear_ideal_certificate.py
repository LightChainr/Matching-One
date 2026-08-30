#!/usr/bin/env python3
"""Verify supplied rational polynomial ideal-infeasibility witnesses."""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import reduce
import hashlib
import json
from math import gcd
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / "results/model-certificates/synthetic/semantic-zero-contradiction/latest.json"
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/linear-ideal/latest.json"
SCHEMA = "matching-one/exact-linear-ideal-certificate/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _lcm(left: int, right: int) -> int:
    return abs(left * right) // gcd(left, right) if left and right else 0


def normalize(polynomial: Sequence[Fraction]) -> list[Fraction]:
    values = list(polynomial) or [Fraction()]
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return values


def add(left: Sequence[Fraction], right: Sequence[Fraction]) -> list[Fraction]:
    size = max(len(left), len(right))
    return normalize(
        [
            (left[index] if index < len(left) else Fraction())
            + (right[index] if index < len(right) else Fraction())
            for index in range(size)
        ]
    )


def multiply(left: Sequence[Fraction], right: Sequence[Fraction]) -> list[Fraction]:
    result = [Fraction()] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return normalize(result)


def verify_ideal_witness(
    constraints_source: Sequence[Sequence[Any]], multipliers_source: Sequence[Sequence[Any]]
) -> Mapping[str, Any]:
    _require(len(constraints_source) == len(multipliers_source) and constraints_source, "constraint/multiplier count mismatch")
    constraints = [[Fraction(value) for value in polynomial] for polynomial in constraints_source]
    multipliers = [[Fraction(value) for value in polynomial] for polynomial in multipliers_source]
    terms = [multiply(multiplier, constraint) for multiplier, constraint in zip(multipliers, constraints)]
    result = reduce(add, terms, [Fraction()])
    _require(result == [Fraction(1)], "supplied ideal witness does not produce constant one")

    multiplier_coefficients = [value for polynomial in multipliers for value in polynomial]
    denominator_lcm = reduce(_lcm, (value.denominator for value in multiplier_coefficients), 1)
    cleared = [int(value * denominator_lcm) for value in multiplier_coefficients]
    content = reduce(gcd, (abs(value) for value in cleared), 0)
    return {
        "term_products": [[str(value) for value in term] for term in terms],
        "result_coefficients_low_to_high": [str(value) for value in result],
        "multiplier_denominator_lcm": denominator_lcm,
        "cleared_multiplier_content": content,
        "primitive_after_common_denominator": content == 1,
        "status": "exact_ideal_contains_one",
    }


def build_result(source_path: Path = DEFAULT_SOURCE) -> dict[str, Any]:
    source = json.loads(source_path.read_text(encoding="utf-8"))
    system = source["polynomial_system"]
    witness = source["primitive_bezout_infeasibility_witness"]
    constraints = [
        system["semantic_zero_coefficients_constant_first"],
        system["observed_row_coefficients_constant_first"],
    ]
    multipliers = [[value] for value in witness["multipliers"]]
    verification = verify_ideal_witness(constraints, multipliers)
    _require(verification["primitive_after_common_denominator"], "frozen witness lost primitive normalization")
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact",
        "source": {
            "path": str(source_path.relative_to(ROOT)),
            "sha256": _sha256_file(source_path),
            "dependency_group": source["synthetic_input"]["dependency_group"],
        },
        "variable_order": [system["variable"]],
        "coefficient_order": "univariate low degree to high degree",
        "constraints": constraints,
        "multipliers": multipliers,
        "verification": verification,
        "solver_invoked": False,
        "claim_boundary": {
            "included": "exact verification that the supplied rational univariate ideal witness combines the declared generators to one",
            "excluded": "witness search, Groebner or SOS completeness, multivariate monomial parsing, noisy constraints, selection-rule derivation, or model validation",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], source_path: Path = DEFAULT_SOURCE) -> Mapping[str, Any]:
    expected = build_result(source_path)
    _require(result == expected, "linear-ideal certificate does not exactly reproduce")
    verification = result["verification"]
    return {
        "schema": result["schema"],
        "status": "valid_exact_linear_ideal_certificate",
        "constraint_count": len(result["constraints"]),
        "result": verification["result_coefficients_low_to_high"],
        "primitive": verification["primitive_after_common_denominator"],
        "solver_invoked": result["solver_invoked"],
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
