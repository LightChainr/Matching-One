#!/usr/bin/env python3
"""Exhaust the frozen degree-two integer-polynomial search exactly."""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from math import gcd
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.exact_polynomial_root_certificate import (
        evaluate,
        isolate_roots,
        open_root_count,
        sturm_sequence,
    )
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from exact_polynomial_root_certificate import evaluate, isolate_roots, open_root_count, sturm_sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "pslq_search_contract.json"
DEFAULT_OUTPUT = ROOT / "results" / "pslq-degree2-polynomial-exclusion" / "latest.json"
SCHEMA = "matching-one/degree2-polynomial-exclusion/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _fraction_text(value: Fraction) -> str:
    return "%d/%d" % (value.numerator, value.denominator)


def _parse_interval(row: Mapping[str, Any]) -> tuple[Fraction, Fraction]:
    try:
        lower = Fraction(row["lower"])
        upper = Fraction(row["upper"])
    except (KeyError, ValueError, ZeroDivisionError) as exc:
        raise ValueError("interval endpoints must be exact decimal strings") from exc
    _require(lower < upper, "interval must be nonempty")
    return lower, upper


def primitive_quadratic_count(height: int) -> int:
    _require(type(height) is int and height >= 1, "height must be a positive integer")
    return sum(
        gcd(gcd(abs(a0), abs(a1)), a2) == 1
        for a2 in range(1, height + 1)
        for a1 in range(-height, height + 1)
        for a0 in range(-height, height + 1)
    )


def _scaled_evaluate(a0: int, a1: int, a2: int, point: Fraction) -> int:
    numerator, denominator = point.numerator, point.denominator
    return a0 * denominator * denominator + a1 * numerator * denominator + a2 * numerator * numerator


def _vertex_inside(a1: int, a2: int, lower: Fraction, upper: Fraction) -> bool:
    # lower <= -a1/(2*a2) <= upper, with a2 positive.
    scale = 2 * a2
    return lower.numerator * scale <= -a1 * lower.denominator and -a1 * upper.denominator <= upper.numerator * scale


def _has_root_in_interval(
    a0: int,
    a1: int,
    a2: int,
    lower: Fraction,
    upper: Fraction,
    lower_value: int,
    upper_value: int,
) -> bool:
    if lower_value == 0 or upper_value == 0 or (lower_value < 0 < upper_value) or (upper_value < 0 < lower_value):
        return True
    if lower_value > 0 and upper_value > 0 and _vertex_inside(a1, a2, lower, upper):
        return 4 * a0 * a2 - a1 * a1 <= 0
    return False


def _minimum_residual(
    a0: int,
    a1: int,
    a2: int,
    lower: Fraction,
    upper: Fraction,
    lower_value: int,
    upper_value: int,
) -> Fraction:
    candidates = [
        Fraction(abs(lower_value), lower.denominator * lower.denominator),
        Fraction(abs(upper_value), upper.denominator * upper.denominator),
    ]
    if _vertex_inside(a1, a2, lower, upper):
        candidates.append(Fraction(abs(4 * a0 * a2 - a1 * a1), 4 * a2))
    return min(candidates)


def _minimum_residual_pair(
    a0: int,
    a1: int,
    a2: int,
    lower: Fraction,
    upper: Fraction,
    lower_value: int,
    upper_value: int,
) -> tuple[int, int]:
    """Return an unreduced exact nonnegative numerator/denominator pair."""

    candidates = [
        (abs(lower_value), lower.denominator * lower.denominator),
        (abs(upper_value), upper.denominator * upper.denominator),
    ]
    if _vertex_inside(a1, a2, lower, upper):
        candidates.append((abs(4 * a0 * a2 - a1 * a1), 4 * a2))
    winner = candidates[0]
    for candidate in candidates[1:]:
        if candidate[0] * winner[1] < winner[0] * candidate[1]:
            winner = candidate
    return winner


def _polynomial_range(a0: int, a1: int, a2: int, lower: Fraction, upper: Fraction) -> tuple[Fraction, Fraction]:
    values = [
        Fraction(_scaled_evaluate(a0, a1, a2, lower), lower.denominator**2),
        Fraction(_scaled_evaluate(a0, a1, a2, upper), upper.denominator**2),
    ]
    if _vertex_inside(a1, a2, lower, upper):
        values.append(Fraction(4 * a0 * a2 - a1 * a1, 4 * a2))
    return min(values), max(values)


def _nearest_root_certificate(coefficients: Sequence[int], lower: Fraction, upper: Fraction) -> Mapping[str, Any]:
    polynomial = [Fraction(value) for value in coefficients]
    sequence = sturm_sequence(polynomial)
    _require(evaluate(polynomial, lower) != 0 and evaluate(polynomial, upper) != 0, "near-case endpoint unexpectedly is a root")
    interval_count = open_root_count(sequence, lower, upper)
    roots = isolate_roots(polynomial, Fraction(0), Fraction(1), bits=60)
    separated = []
    for lo, hi in roots:
        if hi < lower:
            separated.append((lower - hi, lo, hi, "below"))
        elif lo > upper:
            separated.append((lo - upper, lo, hi, "above"))
    nearest = min(separated) if separated else None
    return {
        "sturm_open_root_count_in_method_interval": interval_count,
        "unit_interval_root_brackets": [[_fraction_text(lo), _fraction_text(hi)] for lo, hi in roots],
        "nearest_root_separation_lower_bound": None if nearest is None else _fraction_text(nearest[0]),
        "nearest_root_side": None if nearest is None else nearest[3],
        "isolation_bits": 60,
    }


def run_search(contract: Mapping[str, Any]) -> tuple[int, list[dict[str, Any]]]:
    stage = contract.get("search_stages", {}).get("algebraic_polynomial", {})
    _require(stage.get("degree_min") == 1 and stage.get("degree_max") >= 2, "quadratic stage is not frozen")
    _require(stage.get("primitive_coefficients_only") is True, "primitive coefficients are required")
    _require(stage.get("search_each_interval_separately") is True, "separate interval searches are required")
    height = stage.get("coefficient_height_max")
    _require(type(height) is int and 1 <= height <= 100, "unexpected coefficient height")
    interval_rows = contract.get("intervals", [])
    intervals = [_parse_interval(row) for row in interval_rows]
    best: list[tuple[int, int, tuple[int, int, int]] | None] = [None] * len(intervals)
    containing = [0] * len(intervals)
    checked = 0

    for a2 in range(1, height + 1):
        for a1 in range(-height, height + 1):
            for a0 in range(-height, height + 1):
                if gcd(gcd(abs(a0), abs(a1)), a2) != 1:
                    continue
                checked += 1
                for index, (lower, upper) in enumerate(intervals):
                    lower_value = _scaled_evaluate(a0, a1, a2, lower)
                    upper_value = _scaled_evaluate(a0, a1, a2, upper)
                    if _has_root_in_interval(a0, a1, a2, lower, upper, lower_value, upper_value):
                        containing[index] += 1
                    residual_num, residual_den = _minimum_residual_pair(
                        a0, a1, a2, lower, upper, lower_value, upper_value
                    )
                    current = best[index]
                    coefficients = (a0, a1, a2)
                    if (
                        current is None
                        or residual_num * current[1] < current[0] * residual_den
                        or (
                            residual_num * current[1] == current[0] * residual_den
                            and coefficients < current[2]
                        )
                    ):
                        best[index] = (residual_num, residual_den, coefficients)

    results = []
    for row, (lower, upper), winner, root_count in zip(interval_rows, intervals, best, containing):
        assert winner is not None
        residual_num, residual_den, coefficients = winner
        residual = Fraction(residual_num, residual_den)
        a0, a1, a2 = coefficients
        range_lower, range_upper = _polynomial_range(a0, a1, a2, lower, upper)
        sturm = _nearest_root_certificate(coefficients, lower, upper)
        _require(sturm["sturm_open_root_count_in_method_interval"] == 0, "closest polynomial has an interval root")
        results.append(
            {
                "interval_id": row["id"],
                "source_id": row["source_id"],
                "lower": row["lower"],
                "upper": row["upper"],
                "primitive_quadratics_checked": checked,
                "root_containing_polynomials": root_count,
                "excluded": root_count == 0,
                "closest_polynomial": {
                    "coefficients_ascending": list(coefficients),
                    "height": max(abs(value) for value in coefficients),
                    "discriminant": a1 * a1 - 4 * a0 * a2,
                    "minimum_absolute_residual": _fraction_text(residual),
                    "polynomial_range": [_fraction_text(range_lower), _fraction_text(range_upper)],
                    "independent_sturm_certificate": sturm,
                },
            }
        )
    return checked, results


@lru_cache(maxsize=4)
def build_result(contract_path: Path = DEFAULT_CONTRACT) -> dict[str, Any]:
    contract_bytes = contract_path.read_bytes()
    contract = json.loads(contract_bytes)
    _require(contract.get("schema") == "matching-one/pslq-search-contract/v1", "wrong contract schema")
    provenance = contract.get("provenance", {})
    provenance_path = ROOT / provenance["path"]
    provenance_digest = hashlib.sha256(provenance_path.read_bytes()).hexdigest()
    _require(provenance_digest == provenance.get("sha256"), "provenance digest drift")
    checked, interval_results = run_search(contract)
    return {
        "schema": SCHEMA,
        "issue": 1,
        "status": "degree2_exact_exclusion_complete",
        "contract": {
            "path": str(contract_path.relative_to(ROOT)),
            "sha256": hashlib.sha256(contract_bytes).hexdigest(),
            "provenance_sha256": provenance_digest,
        },
        "search": {
            "degree": 2,
            "coefficient_height_max": contract["search_stages"]["algebraic_polynomial"]["coefficient_height_max"],
            "primitive_coefficients_only": True,
            "sign_normalization": "a2_positive",
            "arithmetic": "exact integer-scaled rational endpoint and vertex tests",
            "unique_polynomials_checked_per_interval": checked,
        },
        "interval_results": interval_results,
        "conclusion": {
            "all_method_intervals_excluded": all(row["excluded"] for row in interval_results),
            "meaning": "no primitive quadratic in the frozen height bound has a real root in any declared method-specific interval",
            "does_not_imply": [
                "exclusion at algebraic degree three or four",
                "exclusion of either frozen constant library",
                "transcendence or absence of a more complex exact representation",
            ],
        },
        "claim_boundary": {
            "included": "exact exhaustive degree-two integer-polynomial exclusion at frozen height 100",
            "excluded": "degrees 3-4, constant-basis PSLQ, near-hit promotion, false-positive calibration, or closed-form claims",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], contract_path: Path = DEFAULT_CONTRACT) -> Mapping[str, Any]:
    expected = build_result() if contract_path == DEFAULT_CONTRACT else build_result(contract_path)
    _require(result == expected, "result does not exactly reproduce the frozen quadratic search")
    return {
        "schema": result["schema"],
        "status": "valid_exact_degree2_exclusion",
        "interval_count": len(result["interval_results"]),
        "quadratics_per_interval": result["search"]["unique_polynomials_checked_per_interval"],
        "all_method_intervals_excluded": result["conclusion"]["all_method_intervals_excluded"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        result = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_result(result, args.contract), indent=2, sort_keys=True))
        return 0
    result = build_result(args.contract)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
