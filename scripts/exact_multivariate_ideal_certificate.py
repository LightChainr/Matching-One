#!/usr/bin/env python3
"""Verify supplied sparse multivariate rational Bezout identities exactly."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/multivariate-ideal/latest.json"
SCHEMA = "matching-one/exact-multivariate-ideal-certificate/v1"
Polynomial = dict[tuple[int, ...], Fraction]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def parse_polynomial(terms: Sequence[Mapping[str, Any]], variables: Sequence[str]) -> Polynomial:
    _require(isinstance(terms, list), "polynomial terms must be a list")
    result: Polynomial = {}
    for term in terms:
        _require(set(term) == {"coefficient", "powers"}, "term fields drift")
        powers = term["powers"]
        _require(isinstance(powers, dict) and set(powers).issubset(variables), "unknown polynomial variable")
        _require(all(isinstance(value, int) and value >= 0 for value in powers.values()), "polynomial exponents must be nonnegative integers")
        monomial = tuple(powers.get(variable, 0) for variable in variables)
        result[monomial] = result.get(monomial, Fraction()) + Fraction(term["coefficient"])
        if result[monomial] == 0:
            del result[monomial]
    return result


def add(left: Polynomial, right: Polynomial) -> Polynomial:
    result = dict(left)
    for monomial, coefficient in right.items():
        result[monomial] = result.get(monomial, Fraction()) + coefficient
        if result[monomial] == 0:
            del result[monomial]
    return result


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result: Polynomial = {}
    for left_monomial, left_coefficient in left.items():
        for right_monomial, right_coefficient in right.items():
            monomial = tuple(a + b for a, b in zip(left_monomial, right_monomial))
            result[monomial] = result.get(monomial, Fraction()) + left_coefficient * right_coefficient
            if result[monomial] == 0:
                del result[monomial]
    return result


def render(polynomial: Polynomial, variables: Sequence[str]) -> list[dict[str, Any]]:
    return [
        {
            "coefficient": str(polynomial[monomial]),
            "powers": {variable: exponent for variable, exponent in zip(variables, monomial) if exponent},
        }
        for monomial in sorted(polynomial)
    ]


def verify_witness(descriptor: Mapping[str, Any]) -> Mapping[str, Any]:
    _require(set(descriptor) == {"variables", "constraints", "multipliers"}, "descriptor fields drift")
    variables = descriptor["variables"]
    _require(isinstance(variables, list) and variables and all(isinstance(variable, str) and variable.isidentifier() for variable in variables) and len(variables) == len(set(variables)), "variables must be valid, nonempty, and unique")
    constraints = descriptor["constraints"]
    multipliers = descriptor["multipliers"]
    _require(isinstance(constraints, list) and len(constraints) == len(multipliers) and constraints, "constraint/multiplier count mismatch")
    products = []
    total: Polynomial = {}
    for constraint_terms, multiplier_terms in zip(constraints, multipliers):
        product = multiply(parse_polynomial(constraint_terms, variables), parse_polynomial(multiplier_terms, variables))
        products.append(product)
        total = add(total, product)
    one = {(0,) * len(variables): Fraction(1)}
    _require(total == one, "supplied multivariate ideal witness does not produce one")
    return {
        "variable_order": variables,
        "constraint_count": len(constraints),
        "term_products": [render(product, variables) for product in products],
        "result": render(total, variables),
        "status": "exact_multivariate_ideal_contains_one",
    }


def frozen_descriptor() -> dict[str, Any]:
    return {
        "variables": ["x", "y"],
        "constraints": [
            [{"coefficient": "1", "powers": {"x": 1}}],
            [{"coefficient": "1", "powers": {}}, {"coefficient": "-1", "powers": {"x": 1, "y": 1}}],
        ],
        "multipliers": [
            [{"coefficient": "1", "powers": {"y": 1}}],
            [{"coefficient": "1", "powers": {}}],
        ],
    }


def build_result() -> dict[str, Any]:
    descriptor = frozen_descriptor()
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact",
        "dependency_group": "synthetic-exact-multivariate-bezout",
        "descriptor": descriptor,
        "verification": verify_witness(descriptor),
        "solver_invoked": False,
        "claim_boundary": {
            "included": "exact verification that supplied sparse rational multivariate multipliers combine supplied equality generators to one",
            "excluded": "witness search, Groebner bases, ideal membership completeness, inequalities, SOS terms, noisy constraints, or model validation",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    _require(result == expected, "multivariate ideal certificate does not exactly reproduce")
    return {"schema": result["schema"], "status": "valid_exact_multivariate_ideal_certificate", "constraint_count": result["verification"]["constraint_count"]}


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
