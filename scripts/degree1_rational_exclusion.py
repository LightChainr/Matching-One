#!/usr/bin/env python3
"""Exhaust the frozen degree-one integer-polynomial search exactly.

This is one bounded slice of Issue 1.  It does not run PSLQ or inspect the
degree-two through degree-four search space.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "pslq_search_contract.json"
DEFAULT_OUTPUT = ROOT / "results" / "pslq-degree1-rational-exclusion" / "latest.json"
SCHEMA = "matching-one/degree1-rational-exclusion/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _fraction(text: Any, label: str) -> Fraction:
    _require(isinstance(text, str) and text.strip() == text, "%s must be an exact string" % label)
    try:
        return Fraction(text)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError("%s is not an exact decimal/rational" % label) from exc


def _canonical_fraction(value: Fraction) -> str:
    return "%d/%d" % (value.numerator, value.denominator)


def primitive_degree_one_coefficients(height: int) -> tuple[tuple[int, int], ...]:
    """Return one sign-normalized representative of every primitive a0+a1*p."""

    _require(type(height) is int and height >= 1, "height must be a positive integer")
    return tuple(
        (a0, a1)
        for a1 in range(1, height + 1)
        for a0 in range(-height, height + 1)
        if gcd(abs(a0), a1) == 1
    )


def _residual_interval(a0: int, a1: int, lower: Fraction, upper: Fraction) -> tuple[Fraction, Fraction]:
    return a0 + a1 * lower, a0 + a1 * upper


def _distance_from_zero(lower: Fraction, upper: Fraction) -> Fraction:
    if lower <= 0 <= upper:
        return Fraction(0)
    return min(abs(lower), abs(upper))


def build_result(contract_path: Path = DEFAULT_CONTRACT) -> dict[str, Any]:
    contract_bytes = contract_path.read_bytes()
    contract = json.loads(contract_bytes)
    _require(contract.get("schema") == "matching-one/pslq-search-contract/v1", "wrong contract schema")
    stage = contract.get("search_stages", {}).get("algebraic_polynomial", {})
    _require(stage.get("degree_min") == 1, "degree-one search is not frozen")
    _require(stage.get("primitive_coefficients_only") is True, "primitive coefficient policy required")
    _require(stage.get("search_each_interval_separately") is True, "separate interval searches required")
    height = stage.get("coefficient_height_max")
    coefficients = primitive_degree_one_coefficients(height)

    interval_results = []
    for row in contract.get("intervals", []):
        lower = _fraction(row.get("lower"), "%s lower" % row.get("id"))
        upper = _fraction(row.get("upper"), "%s upper" % row.get("id"))
        _require(lower < upper, "interval must be nonempty")
        containing = []
        best = None
        for a0, a1 in coefficients:
            residual_lower, residual_upper = _residual_interval(a0, a1, lower, upper)
            distance = _distance_from_zero(residual_lower, residual_upper)
            candidate = (distance, a1, abs(a0), a0, residual_lower, residual_upper)
            if distance == 0:
                containing.append((a0, a1))
            if best is None or candidate[:4] < best[:4]:
                best = candidate
        assert best is not None
        distance, a1, _, a0, residual_lower, residual_upper = best
        interval_results.append(
            {
                "interval_id": row["id"],
                "source_id": row["source_id"],
                "lower": row["lower"],
                "upper": row["upper"],
                "primitive_polynomials_checked": len(coefficients),
                "zero_containing_residuals": len(containing),
                "excluded": not containing,
                "closest_polynomial": {
                    "coefficients_ascending": [a0, a1],
                    "height": max(abs(a0), a1),
                    "root": _canonical_fraction(Fraction(-a0, a1)),
                    "slope_abs": a1,
                    "root_condition_to_additive_residual": _canonical_fraction(Fraction(1, a1)),
                    "residual_interval": [
                        _canonical_fraction(residual_lower),
                        _canonical_fraction(residual_upper),
                    ],
                    "minimum_absolute_residual": _canonical_fraction(distance),
                },
            }
        )

    provenance = contract.get("provenance", {})
    provenance_path = ROOT / provenance["path"]
    provenance_digest = hashlib.sha256(provenance_path.read_bytes()).hexdigest()
    _require(provenance_digest == provenance.get("sha256"), "contract provenance digest drift")
    return {
        "schema": SCHEMA,
        "issue": 1,
        "status": "degree1_exact_exclusion_complete",
        "contract": {
            "path": str(contract_path.relative_to(ROOT)),
            "sha256": hashlib.sha256(contract_bytes).hexdigest(),
            "provenance_sha256": provenance_digest,
        },
        "search": {
            "degree": 1,
            "coefficient_height_max": height,
            "primitive_coefficients_only": True,
            "sign_normalization": "a1_positive",
            "arithmetic": "fractions.Fraction exact rational endpoints",
            "unique_polynomials_checked_per_interval": len(coefficients),
        },
        "interval_results": interval_results,
        "conclusion": {
            "all_method_intervals_excluded": all(row["excluded"] for row in interval_results),
            "meaning": "no primitive degree-one integer polynomial in the frozen height bound has a root in any declared method-specific interval",
            "does_not_imply": [
                "exclusion at algebraic degree two through four",
                "exclusion of the frozen standard-constant bases",
                "transcendence or absence of a more complex exact representation",
            ],
        },
        "claim_boundary": {
            "included": "exact exhaustive degree-one integer-polynomial exclusion at frozen height 100",
            "excluded": "degrees 2-4, PSLQ constant bases, near-hit promotion, synthetic calibration, or closed-form claims",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], contract_path: Path = DEFAULT_CONTRACT) -> dict[str, Any]:
    expected = build_result(contract_path)
    _require(result == expected, "result does not exactly reproduce the frozen search")
    _require(result.get("claim_boundary", {}).get("parent_issue") == "remain open", "parent boundary drift")
    return {
        "schema": result["schema"],
        "status": "valid_exact_degree1_exclusion",
        "interval_count": len(result["interval_results"]),
        "polynomials_per_interval": result["search"]["unique_polynomials_checked_per_interval"],
        "all_method_intervals_excluded": result["conclusion"]["all_method_intervals_excluded"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate is not None:
        result = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_result(result, args.contract), indent=2, sort_keys=True))
        return 0
    result = build_result(args.contract)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = args.output
    if output is None:
        print(rendered, end="")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
