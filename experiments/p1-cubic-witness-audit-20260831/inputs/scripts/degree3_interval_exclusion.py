#!/usr/bin/env python3
"""Exhaust one frozen degree-three integer-polynomial interval exactly."""

from __future__ import annotations

import argparse
from bisect import bisect_left, bisect_right
from fractions import Fraction
from functools import lru_cache
from math import gcd
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.exact_polynomial_root_certificate import isolate_roots, open_root_count, sturm_sequence
except ModuleNotFoundError:
    from exact_polynomial_root_certificate import isolate_roots, open_root_count, sturm_sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "pslq_search_contract.json"
SCHEMA = "matching-one/degree3-interval-exclusion/v1"
HEIGHT = 100


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def output_path(interval_id: str) -> Path:
    return ROOT / "results" / f"pslq-degree3-{interval_id}" / "latest.json"


def _scaled_cubic_tail(a1: int, a2: int, a3: int, point: Fraction) -> int:
    n, d = point.numerator, point.denominator
    return a1 * n * d * d + a2 * n * n * d + a3 * n * n * n


def _scaled_derivative(a1: int, a2: int, a3: int, point: Fraction) -> int:
    n, d = point.numerator, point.denominator
    return a1 * d * d + 2 * a2 * n * d + 3 * a3 * n * n


def _derivative_has_interval_root(
    a1: int, a2: int, a3: int, lower: Fraction, upper: Fraction
) -> bool:
    left = _scaled_derivative(a1, a2, a3, lower)
    right = _scaled_derivative(a1, a2, a3, upper)
    if left == 0 or right == 0 or left * right < 0:
        return True
    # The derivative is an upward quadratic.  Only two roots between positive
    # endpoints need an additional vertex check.
    vertex_inside = lower * (3 * a3) <= -a2 <= upper * (3 * a3)
    return left > 0 and right > 0 and vertex_inside and 3 * a3 * a1 - a2 * a2 <= 0


@lru_cache(maxsize=2)
def _primitive_a0_values(height: int = HEIGHT) -> dict[int, tuple[int, ...]]:
    return {
        divisor: tuple(a0 for a0 in range(-height, height + 1) if gcd(abs(a0), divisor) == 1)
        for divisor in range(1, height + 1)
    }


def _nearest_root_certificate(
    coefficients: Sequence[int], lower: Fraction, upper: Fraction
) -> Mapping[str, Any]:
    polynomial = [Fraction(value) for value in coefficients]
    sequence = sturm_sequence(polynomial)
    count = open_root_count(sequence, lower, upper)
    roots = isolate_roots(polynomial, Fraction(0), Fraction(1), bits=80)
    separated = []
    for lo, hi in roots:
        if hi < lower:
            separated.append((lower - hi, "below", lo, hi))
        elif lo > upper:
            separated.append((lo - upper, "above", lo, hi))
    nearest = min(separated) if separated else None
    return {
        "sturm_open_root_count_in_method_interval": count,
        "unit_interval_root_brackets": [[_text(lo), _text(hi)] for lo, hi in roots],
        "nearest_root_separation_lower_bound": None if nearest is None else _text(nearest[0]),
        "nearest_root_side": None if nearest is None else nearest[1],
        "isolation_bits": 80,
    }


def run_search(interval: Mapping[str, Any], height: int = HEIGHT) -> dict[str, Any]:
    lower, upper = Fraction(interval["lower"]), Fraction(interval["upper"])
    _require(lower < upper, "interval must be nonempty")
    allowed = _primitive_a0_values(height)
    checked = 0
    containing = 0
    stationary_fibers = 0
    best: tuple[int, int, tuple[int, int, int, int]] | None = None

    for a3 in range(1, height + 1):
        for a2 in range(-height, height + 1):
            for a1 in range(-height, height + 1):
                if _derivative_has_interval_root(a1, a2, a3, lower, upper):
                    stationary_fibers += 1
                    continue

                lower_tail = _scaled_cubic_tail(a1, a2, a3, lower)
                upper_tail = _scaled_cubic_tail(a1, a2, a3, upper)
                lower_den = lower.denominator**3
                upper_den = upper.denominator**3
                if lower_tail * upper_den <= upper_tail * lower_den:
                    qmin_num, qmin_den = lower_tail, lower_den
                    qmax_num, qmax_den = upper_tail, upper_den
                else:
                    qmin_num, qmin_den = upper_tail, upper_den
                    qmax_num, qmax_den = lower_tail, lower_den

                divisor = gcd(gcd(abs(a1), abs(a2)), a3)
                a0_values = allowed[divisor]
                checked += len(a0_values)

                first_root_a0 = _ceil_div(-qmax_num, qmax_den)
                last_root_a0 = (-qmin_num) // qmin_den
                if first_root_a0 <= height and last_root_a0 >= -height:
                    containing += bisect_right(a0_values, min(height, last_root_a0)) - bisect_left(
                        a0_values, max(-height, first_root_a0)
                    )

                target = (-qmin_num) // qmin_den
                insertion = bisect_left(a0_values, target)
                for index in range(max(0, insertion - 2), min(len(a0_values), insertion + 3)):
                    a0 = a0_values[index]
                    low_value = a0 * qmin_den + qmin_num
                    high_value = a0 * qmax_den + qmax_num
                    if low_value <= 0 <= high_value:
                        residual_num, residual_den = 0, 1
                    elif low_value > 0:
                        residual_num, residual_den = low_value, qmin_den
                    else:
                        residual_num, residual_den = -high_value, qmax_den
                    coefficients = (a0, a1, a2, a3)
                    if (
                        best is None
                        or residual_num * best[1] < best[0] * residual_den
                        or (
                            residual_num * best[1] == best[0] * residual_den
                            and coefficients < best[2]
                        )
                    ):
                        best = (residual_num, residual_den, coefficients)

    _require(stationary_fibers == 0, "a derivative stationary point requires a separate exact range path")
    assert best is not None
    residual = Fraction(best[0], best[1])
    coefficients = best[2]
    lo_value = sum(Fraction(value) * lower**power for power, value in enumerate(coefficients))
    hi_value = sum(Fraction(value) * upper**power for power, value in enumerate(coefficients))
    certificate = _nearest_root_certificate(coefficients, lower, upper)
    _require(certificate["sturm_open_root_count_in_method_interval"] == 0, "closest cubic has an interval root")
    return {
        "interval_id": interval["id"],
        "source_id": interval["source_id"],
        "lower": interval["lower"],
        "upper": interval["upper"],
        "primitive_cubics_checked": checked,
        "coefficient_fibers_checked": height * (2 * height + 1) ** 2,
        "derivative_stationary_fibers": stationary_fibers,
        "root_containing_polynomials": containing,
        "excluded": containing == 0,
        "closest_polynomial": {
            "coefficients_ascending": list(coefficients),
            "height": max(abs(value) for value in coefficients),
            "minimum_absolute_residual": _text(residual),
            "polynomial_endpoint_values": [_text(lo_value), _text(hi_value)],
            "independent_sturm_certificate": certificate,
        },
    }


@lru_cache(maxsize=8)
def build_result(interval_id: str, contract_path: Path = DEFAULT_CONTRACT) -> dict[str, Any]:
    raw = contract_path.read_bytes()
    contract = json.loads(raw)
    stage = contract["search_stages"]["algebraic_polynomial"]
    _require(stage["degree_min"] == 1 and stage["degree_max"] >= 3, "cubic stage is not frozen")
    _require(stage["coefficient_height_max"] == HEIGHT, "cubic height drift")
    _require(stage["primitive_coefficients_only"] is True, "primitive coefficients are required")
    rows = [row for row in contract["intervals"] if row["id"] == interval_id]
    _require(len(rows) == 1, "interval id is not uniquely frozen")
    provenance = contract["provenance"]
    provenance_digest = hashlib.sha256((ROOT / provenance["path"]).read_bytes()).hexdigest()
    _require(provenance_digest == provenance["sha256"], "provenance digest drift")
    result = run_search(rows[0])
    return {
        "schema": SCHEMA,
        "issue": 1,
        "status": "degree3_interval_exact_exclusion_complete",
        "contract_sha256": hashlib.sha256(raw).hexdigest(),
        "provenance_sha256": provenance_digest,
        "search": {
            "degree": 3,
            "coefficient_height_max": HEIGHT,
            "primitive_coefficients_only": True,
            "sign_normalization": "a3_positive",
            "arithmetic": "exact integer-scaled rational endpoint ranges with exact derivative monotonicity",
        },
        "interval_result": result,
        "claim_boundary": {
            "included": f"degree-3 height-100 exclusion on {interval_id} only",
            "excluded": "other method intervals, degree 4, library expansion, near-hit promotion, p-values, closed forms, or transcendence",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], interval_id: str) -> Mapping[str, Any]:
    expected = build_result(interval_id)
    _require(result == expected, "degree-3 interval result does not exactly reproduce")
    row = expected["interval_result"]
    return {
        "schema": SCHEMA,
        "status": "valid",
        "interval_id": interval_id,
        "primitive_cubics_checked": row["primitive_cubics_checked"],
        "excluded": row["excluded"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("interval_id")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_result(json.loads(args.validate.read_text()), args.interval_id), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_result(args.interval_id), indent=2, sort_keys=True) + "\n"
    destination = args.output
    if destination:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(rendered)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
