#!/usr/bin/env python3
"""Regenerate the exact `(a,d)=(3,2)` Pell polynomial regression pair."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any, Mapping

from exact_matching_polynomial import bernstein_counts, bernstein_to_power
from matched_torus_reference import axis_geometry, diamond_geometry


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "pell_3_2_exact_regression.json"
EXPECTED_SCHEMA = "matching-one/pell-3-2-exact-regression/v1"
FORBIDDEN_KEYS = {
    "effective_exponent",
    "estimate",
    "monte_carlo",
    "physical_root",
    "root_gap",
    "samples",
    "standard_error",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _walk_forbidden(value: Any, path: str = "$") -> None:
    if isinstance(value, dict):
        bad = sorted(set(value) & FORBIDDEN_KEYS)
        _require(not bad, f"{path} contains forbidden result keys: {','.join(bad)}")
        for key, child in value.items():
            _walk_forbidden(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _walk_forbidden(child, f"{path}[{index}]")


def _primitive(coefficients: list[int]) -> bool:
    content = 0
    for coefficient in coefficients:
        content = math.gcd(content, abs(coefficient))
    return content == 1


def _regenerate(kind: str, length: int) -> dict[str, Any]:
    if kind == "axis":
        geometry = axis_geometry(length)
    elif kind == "diamond":
        geometry = diamond_geometry(length)
    else:
        raise ValueError(f"unknown geometry kind: {kind}")
    bernstein = bernstein_counts(geometry)
    power = bernstein_to_power(bernstein)
    return {
        "kind": kind,
        "L": length,
        "N": geometry.n,
        "bernstein_counts": bernstein,
        "power_coefficients_ascending": power,
    }


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _walk_forbidden(contract)
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 19, "issue must be 19")
    _require(contract.get("status") == "exact_regression_only", "status drifted")
    _require(
        contract.get("source")
        == {
            "enumerator": "scripts/exact_matching_polynomial.py",
            "geometry": "scripts/matched_torus_reference.py",
        },
        "source declaration drifted",
    )

    pair = contract.get("pell_pair", {})
    a, d = pair.get("axis_a"), pair.get("diamond_d")
    _require((a, d) == (3, 2), "Pell pair must remain (3,2)")
    residual = a * a - 2 * d * d
    _require(residual == 1 and pair.get("residual") == residual, "Pell residual drifted")
    axis_sites, diamond_sites = a * a, 2 * d * d
    _require(pair.get("axis_sites") == axis_sites, "axis site count drifted")
    _require(pair.get("diamond_sites") == diamond_sites, "diamond site count drifted")
    _require(axis_sites - diamond_sites == 1, "Pell site counts must differ by one")
    squared_ratio = Fraction(axis_sites, diamond_sites)
    _require(pair.get("squared_length_ratio") == str(squared_ratio), "length ratio drifted")

    expected_layout = [("axis_L3", "axis", 3, 9), ("diamond_L2", "diamond", 2, 8)]
    geometries = contract.get("geometries")
    _require(isinstance(geometries, list) and len(geometries) == 2, "two geometry regressions are required")
    rows = []
    for declared, layout in zip(geometries, expected_layout):
        identifier, kind, length, site_count = layout
        _require(
            (declared.get("id"), declared.get("kind"), declared.get("L"), declared.get("N"))
            == layout,
            f"{identifier} declaration drifted",
        )
        regenerated = _regenerate(kind, length)
        _require(regenerated["N"] == site_count, f"{identifier} regenerated site count drifted")
        _require(
            declared.get("bernstein_counts") == regenerated["bernstein_counts"],
            f"{identifier} Bernstein coefficients drifted",
        )
        coefficients = regenerated["power_coefficients_ascending"]
        _require(
            declared.get("power_coefficients_ascending") == coefficients,
            f"{identifier} power coefficients drifted",
        )
        _require(len(coefficients) - 1 == site_count, f"{identifier} degree is not N")
        _require(_primitive(coefficients), f"{identifier} polynomial is not primitive")
        rows.append(
            {
                "id": identifier,
                "N": site_count,
                "degree": len(coefficients) - 1,
                "primitive": True,
                "bernstein_coefficients": len(regenerated["bernstein_counts"]),
            }
        )

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_exact_regression",
        "pell_residual": residual,
        "site_difference": axis_sites - diamond_sites,
        "squared_length_ratio": str(squared_ratio),
        "geometries": rows,
        "numerical_roots_evaluated": False,
        "contains_monte_carlo_result": False,
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
