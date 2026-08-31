#!/usr/bin/env python3
"""Verify supplied rational interval arithmetic and algebraic enclosures exactly."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/interval-enclosure/latest.json"
SCHEMA = "matching-one/exact-rational-interval-enclosure/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def interval(values: Sequence[Any]) -> tuple[Fraction, Fraction]:
    _require(isinstance(values, list) and len(values) == 2, "interval must have two endpoints")
    lower, upper = (Fraction(value) for value in values)
    _require(lower <= upper, "interval endpoints are reversed")
    return lower, upper


def add(left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return left[0] + right[0], left[1] + right[1]


def multiply(left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    products = [left_endpoint * right_endpoint for left_endpoint in left for right_endpoint in right]
    return min(products), max(products)


def _render(value: tuple[Fraction, Fraction]) -> list[str]:
    return [str(endpoint) for endpoint in value]


def verify_enclosure(descriptor: Mapping[str, Any]) -> Mapping[str, Any]:
    _require(set(descriptor) == {"left", "right", "declared_sum", "declared_product", "sqrt_radicand", "sqrt_interval"}, "descriptor fields drift")
    left = interval(descriptor["left"])
    right = interval(descriptor["right"])
    declared_sum = interval(descriptor["declared_sum"])
    declared_product = interval(descriptor["declared_product"])
    actual_sum = add(left, right)
    actual_product = multiply(left, right)
    _require(actual_sum == declared_sum, "declared interval sum is not exact")
    _require(actual_product == declared_product, "declared interval product is not exact")
    radicand = Fraction(descriptor["sqrt_radicand"])
    sqrt_interval = interval(descriptor["sqrt_interval"])
    _require(radicand >= 0 and sqrt_interval[0] >= 0, "square-root enclosure must be nonnegative")
    lower_square = sqrt_interval[0] ** 2
    upper_square = sqrt_interval[1] ** 2
    _require(lower_square <= radicand <= upper_square, "square-root interval does not enclose radicand")
    return {
        "sum": _render(actual_sum),
        "product": _render(actual_product),
        "sqrt_interval": _render(sqrt_interval),
        "sqrt_lower_square": str(lower_square),
        "sqrt_radicand": str(radicand),
        "sqrt_upper_square": str(upper_square),
        "all_checks_exact": True,
        "status": "exact_rational_interval_enclosure_verified",
    }


def frozen_descriptor() -> dict[str, Any]:
    return {
        "left": ["1/3", "1/2"],
        "right": ["-1/4", "2/3"],
        "declared_sum": ["1/12", "7/6"],
        "declared_product": ["-1/8", "1/3"],
        "sqrt_radicand": "2",
        "sqrt_interval": ["707/500", "283/200"],
    }


def build_result() -> dict[str, Any]:
    descriptor = frozen_descriptor()
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact_control",
        "dependency_group": "synthetic-exact-rational-intervals",
        "descriptor": descriptor,
        "verification": verify_enclosure(descriptor),
        "claim_boundary": {
            "included": "exact rational verification of supplied interval addition, multiplication, and one positive square-root enclosure",
            "excluded": "automatic outward rounding, transcendental constants, general algebraic-number isolation, measured-input enclosure, or statistical coverage",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    _require(result == expected, "interval certificate does not exactly reproduce")
    return {"schema": result["schema"], "status": "valid_exact_interval_enclosure", "sqrt_interval": result["verification"]["sqrt_interval"]}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_result(json.loads(args.validate.read_text(encoding="utf-8"))), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_result(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
