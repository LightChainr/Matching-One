#!/usr/bin/env python3
"""Validate an exact coordinate-free finite jet for the parametric curve M(U)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

from exact_jet_algebra import (
    normalized_odd_invariants,
    parametric_derivatives,
    parse_derivatives,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "parametric_matching_profile_contract.json"
EXPECTED_SCHEMA = "matching-one/parametric-matching-profile/v1"
FORBIDDEN_KEYS = {
    "batches",
    "covariance",
    "estimate",
    "measurements",
    "samples",
    "standard_error",
    "universality",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _walk_forbidden(value: Any, path: str = "$") -> None:
    if isinstance(value, dict):
        bad = sorted(set(value) & FORBIDDEN_KEYS)
        _require(not bad, f"{path} contains forbidden empirical keys: {','.join(bad)}")
        for key, child in value.items():
            _walk_forbidden(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _walk_forbidden(child, f"{path}[{index}]")


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _walk_forbidden(contract)
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 25, "issue must be 25")
    _require(contract.get("status") == "exact_algebra_only", "status drifted")
    order = contract.get("maximum_order")
    _require(order == 5, "maximum order must remain five")

    expected = parse_derivatives(
        contract.get("expected_parametric_derivatives", []), order, "expected parametric jet"
    )
    expected_invariants = {
        key: str(value) for key, value in normalized_odd_invariants(expected).items()
    }
    _require(
        contract.get("expected_normalized_odd_invariants") == expected_invariants,
        "expected normalized invariants drifted",
    )

    representations = contract.get("representations")
    _require(isinstance(representations, list) and len(representations) == 2, "two representations are required")
    seen: set[str] = set()
    rows = []
    for representation in representations:
        identifier = representation.get("id")
        _require(isinstance(identifier, str) and identifier not in seen, "representation id is missing or repeated")
        seen.add(identifier)
        matching = parse_derivatives(
            representation.get("matching_derivatives", []), order, f"{identifier} matching jet"
        )
        reference = parse_derivatives(
            representation.get("reference_derivatives", []), order, f"{identifier} reference jet"
        )
        derived = parametric_derivatives(matching, reference, order)
        _require(derived == expected, f"{identifier} does not reconstruct the frozen M(U) jet")
        rows.append(
            {
                "id": identifier,
                "parametric_derivatives": [str(value) for value in derived],
                "normalized_odd_invariants": {
                    key: str(value) for key, value in normalized_odd_invariants(derived).items()
                },
            }
        )

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_exact_parametric_jet",
        "maximum_order": order,
        "representations_verified": rows,
        "common_parametric_derivatives": [str(value) for value in expected],
        "normalized_odd_invariants": expected_invariants,
        "contains_empirical_result": False,
        "parent_issue": "remain open",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    print(json.dumps(validate_contract(contract), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
