#!/usr/bin/env python3
"""Compile rational parameter boxes into canonical polynomial inequalities."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/parameter-bound-problem/latest.json"
SCHEMA = "matching-one/canonical-parameter-bound-problem/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _term(coefficient: Fraction, exponents: Sequence[int]) -> Mapping[str, Any]:
    return {"coefficient": str(coefficient), "exponents": list(exponents)}


def compile_bounds(bounds: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    _require(isinstance(bounds, list) and bounds, "bounds cannot be empty")
    names = [bound.get("name") for bound in bounds]
    _require(all(isinstance(name, str) and name for name in names) and len(names) == len(set(names)), "bound names must be nonempty and unique")
    dimension = len(bounds)
    inequalities = []
    radius_squared = Fraction()
    normalized_bounds = []
    for index, bound in enumerate(bounds):
        _require(set(bound) == {"name", "lower", "upper", "provenance"}, "bound fields drift")
        _require(isinstance(bound["provenance"], str) and bound["provenance"], "bound provenance is required")
        lower = Fraction(bound["lower"])
        upper = Fraction(bound["upper"])
        _require(lower <= upper, f"reversed bounds for {bound['name']}")
        unit = [0] * dimension
        unit[index] = 1
        zero = [0] * dimension
        inequalities.append({"id": f"{bound['name']}_lower", "meaning": f"{bound['name']}-lower >= 0", "terms": [_term(-lower, zero), _term(Fraction(1), unit)]})
        inequalities.append({"id": f"{bound['name']}_upper", "meaning": f"upper-{bound['name']} >= 0", "terms": [_term(upper, zero), _term(Fraction(-1), unit)]})
        coordinate_radius = max(lower * lower, upper * upper)
        radius_squared += coordinate_radius
        normalized_bounds.append({"name": bound["name"], "lower": str(lower), "upper": str(upper), "coordinate_square_bound": str(coordinate_radius), "provenance": bound["provenance"]})
    ball_terms = [_term(radius_squared, [0] * dimension)]
    for index in range(dimension):
        square = [0] * dimension
        square[index] = 2
        ball_terms.append(_term(Fraction(-1), square))
    inequalities.append({"id": "archimedean_ball", "meaning": "R-sum(x_i^2) >= 0 on the declared box", "terms": ball_terms})
    return {
        "coefficient_field": "Q",
        "variable_order": names,
        "bounds": normalized_bounds,
        "inequalities": inequalities,
        "inequality_count": len(inequalities),
        "archimedean_radius_squared": str(radius_squared),
        "derivation": "sum of exact coordinate-wise max(lower^2, upper^2) bounds",
        "status": "canonical_parameter_bound_problem_compiled",
    }


def evaluate(terms: Sequence[Mapping[str, Any]], values: Sequence[Any]) -> Fraction:
    rational_values = [Fraction(value) for value in values]
    total = Fraction()
    for term in terms:
        _require(len(term["exponents"]) == len(rational_values), "evaluation dimension mismatch")
        value = Fraction(term["coefficient"])
        for base, exponent in zip(rational_values, term["exponents"]):
            value *= base ** exponent
        total += value
    return total


def frozen_bounds() -> list[dict[str, str]]:
    return [
        {"name": "x", "lower": "-1", "upper": "2", "provenance": "synthetic compact x chart"},
        {"name": "y", "lower": "-1/2", "upper": "1/3", "provenance": "synthetic compact y chart"},
    ]


def build_result() -> dict[str, Any]:
    bounds = frozen_bounds()
    problem = compile_bounds(bounds)
    assignment = ["1/2", "0"]
    evaluations = [evaluate(item["terms"], assignment) for item in problem["inequalities"]]
    _require(all(value >= 0 for value in evaluations), "frozen assignment violates compiled bounds")
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact_compilation",
        "dependency_group": "synthetic-exact-parameter-bound-problem",
        "input_bounds": bounds,
        "problem": problem,
        "supplied_assignment_check": {"assignment": assignment, "inequality_values": [str(value) for value in evaluations], "all_nonnegative": True},
        "claim_boundary": {
            "included": "deterministic compilation of supplied rational boxes into canonical linear inequalities plus a derived Archimedean ball inequality",
            "excluded": "correlated constraints, optimal compactification, physical-bound derivation, general inequality parsing, SOS relaxation generation, or feasibility solving",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    _require(result == expected, "parameter-bound problem does not exactly reproduce")
    return {"schema": result["schema"], "status": "valid_canonical_parameter_bound_problem", "inequality_count": result["problem"]["inequality_count"], "radius_squared": result["problem"]["archimedean_radius_squared"]}


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
