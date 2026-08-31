#!/usr/bin/env python3
"""Verify supplied rational parameter assignments and declared box bounds."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/parameter-box/latest.json"
SCHEMA = "matching-one/exact-parameter-box-certificate/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def verify_parameter_box(parameters: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    _require(isinstance(parameters, list) and parameters, "parameter list cannot be empty")
    seen = set()
    verified = []
    for parameter in parameters:
        _require(set(parameter) == {"name", "value", "lower", "upper", "provenance"}, "parameter fields drift")
        name = parameter["name"]
        _require(isinstance(name, str) and name and name not in seen, "parameter names must be nonempty and unique")
        seen.add(name)
        _require(isinstance(parameter["provenance"], str) and parameter["provenance"], "bound provenance is required")
        value = Fraction(parameter["value"])
        lower = Fraction(parameter["lower"])
        upper = Fraction(parameter["upper"])
        _require(lower <= upper, f"reversed bounds for {name}")
        _require(lower <= value <= upper, f"parameter outside bounds: {name}")
        verified.append({
            "name": name,
            "value": str(value),
            "lower": str(lower),
            "upper": str(upper),
            "lower_slack": str(value - lower),
            "upper_slack": str(upper - value),
            "active_lower": value == lower,
            "active_upper": value == upper,
            "provenance": parameter["provenance"],
        })
    active = []
    for item in verified:
        if item["active_lower"]:
            active.append(item["name"] + ":lower")
        if item["active_upper"]:
            active.append(item["name"] + ":upper")
    return {"parameter_count": len(verified), "parameters": verified, "active_boundaries": active, "status": "exact_parameter_box_verified"}


def frozen_parameters() -> list[dict[str, str]]:
    return [
        {"name": "x", "value": "1/2", "lower": "-1", "upper": "1", "provenance": "synthetic invariant-coordinate normalization"},
        {"name": "k", "value": "0", "lower": "0", "upper": "2", "provenance": "synthetic nonnegative coupling cap"},
        {"name": "lambda", "value": "3/2", "lower": "1", "upper": "3/2", "provenance": "synthetic frozen transfer interval"},
    ]


def build_result() -> dict[str, Any]:
    parameters = frozen_parameters()
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact_control",
        "dependency_group": "synthetic-exact-parameter-box",
        "parameters": parameters,
        "verification": verify_parameter_box(parameters),
        "claim_boundary": {
            "included": "exact verification of supplied rational parameter values, box bounds, provenance strings, slacks, and active faces",
            "excluded": "derivation or exhaustion of physical bounds, unbounded parameters, correlated semialgebraic domains, model feasibility, or statistical uncertainty",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    _require(result == expected, "parameter-box certificate does not exactly reproduce")
    return {"schema": result["schema"], "status": "valid_exact_parameter_box", "active_boundaries": result["verification"]["active_boundaries"]}


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
