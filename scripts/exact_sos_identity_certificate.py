#!/usr/bin/env python3
"""Verify supplied rational SOS and inequality-localizing identities exactly."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/sos-identity/latest.json"
SCHEMA = "matching-one/exact-sos-localizing-identity/v1"
Polynomial = dict[tuple[int, ...], Fraction]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def parse(terms: Sequence[Mapping[str, Any]], variables: Sequence[str]) -> Polynomial:
    _require(isinstance(terms, list), "polynomial terms must be a list")
    result: Polynomial = {}
    for term in terms:
        _require(set(term) == {"coefficient", "powers"}, "term fields drift")
        powers = term["powers"]
        _require(isinstance(powers, dict) and set(powers).issubset(variables), "unknown polynomial variable")
        _require(all(isinstance(value, int) and value >= 0 for value in powers.values()), "invalid polynomial exponent")
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
    for lm, lc in left.items():
        for rm, rc in right.items():
            monomial = tuple(a + b for a, b in zip(lm, rm))
            result[monomial] = result.get(monomial, Fraction()) + lc * rc
            if result[monomial] == 0:
                del result[monomial]
    return result


def render(polynomial: Polynomial, variables: Sequence[str]) -> list[dict[str, Any]]:
    return [{"coefficient": str(polynomial[m]), "powers": {v: e for v, e in zip(variables, m) if e}} for m in sorted(polynomial)]


def sum_of_squares(square_terms: Sequence[Sequence[Mapping[str, Any]]], variables: Sequence[str]) -> Polynomial:
    _require(isinstance(square_terms, list), "square terms must be a list")
    total: Polynomial = {}
    for terms in square_terms:
        polynomial = parse(terms, variables)
        total = add(total, multiply(polynomial, polynomial))
    return total


def verify_identity(descriptor: Mapping[str, Any]) -> Mapping[str, Any]:
    _require(set(descriptor) == {"variables", "target", "free_squares", "inequalities", "localizing_squares"}, "descriptor fields drift")
    variables = descriptor["variables"]
    _require(isinstance(variables, list) and variables and all(isinstance(variable, str) and variable.isidentifier() for variable in variables) and len(variables) == len(set(variables)), "variables must be valid, nonempty, and unique")
    inequalities = descriptor["inequalities"]
    localizers = descriptor["localizing_squares"]
    _require(isinstance(inequalities, list) and len(inequalities) == len(localizers), "inequality/localizer count mismatch")
    total = sum_of_squares(descriptor["free_squares"], variables)
    free_sos = total
    localizing_products = []
    for inequality_terms, square_terms in zip(inequalities, localizers):
        sigma = sum_of_squares(square_terms, variables)
        product = multiply(sigma, parse(inequality_terms, variables))
        localizing_products.append(product)
        total = add(total, product)
    target = parse(descriptor["target"], variables)
    _require(total == target, "supplied SOS/localizing identity does not match target")
    return {
        "variable_order": variables,
        "free_sos": render(free_sos, variables),
        "localizing_products": [render(product, variables) for product in localizing_products],
        "result": render(total, variables),
        "inequality_count": len(inequalities),
        "status": "exact_sos_localizing_identity_verified",
    }


def frozen_descriptor() -> dict[str, Any]:
    one = [{"coefficient": "1", "powers": {}}]
    return {
        "variables": ["x", "y"],
        "target": [{"coefficient": "1", "powers": {}}, {"coefficient": "2", "powers": {"x": 1, "y": 1}}, {"coefficient": "1", "powers": {"y": 2}}],
        "free_squares": [[{"coefficient": "1", "powers": {"x": 1}}, {"coefficient": "1", "powers": {"y": 1}}]],
        "inequalities": [[{"coefficient": "1", "powers": {}}, {"coefficient": "-1", "powers": {"x": 2}}]],
        "localizing_squares": [[one]],
    }


def build_result() -> dict[str, Any]:
    descriptor = frozen_descriptor()
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact",
        "dependency_group": "synthetic-exact-sos-localizer",
        "descriptor": descriptor,
        "verification": verify_identity(descriptor),
        "solver_invoked": False,
        "claim_boundary": {
            "included": "exact polynomial verification of supplied rational square terms and supplied inequality localizers",
            "excluded": "SOS search, SDP solving, domain sufficiency, Archimedean completeness, equality multipliers, numerical reconstruction, or physical positivity",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    _require(result == expected, "SOS identity certificate does not exactly reproduce")
    return {"schema": result["schema"], "status": "valid_exact_sos_localizing_identity", "inequality_count": result["verification"]["inequality_count"]}


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
