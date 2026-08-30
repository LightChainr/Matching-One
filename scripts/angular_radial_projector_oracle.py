#!/usr/bin/env python3
"""Exact tensor-product oracle for angular and radial H4 filters.

The bounded model has two orientation samples and two sizes.  Angular
averaging kills the declared spin-4 orientation vector ``(1, -1)``.  The
radial row ``(-r, 1)`` kills a declared size vector ``(1/r, 1)``.  Because
the two maps act on different tensor factors, their compositions commute.

This is a synthetic linear-algebra gate.  It does not assert that a measured
observable is a pure separable H4 sector or determine a scaling exponent.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "angular_radial_projector_contract.json"
EXPECTED_SCHEMA = "matching-one/angular-radial-projector/v1"
Vector = tuple[Fraction, Fraction]
Grid = tuple[Vector, Vector]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def parse_fraction(value: Any, label: str) -> Fraction:
    _require(isinstance(value, (str, int)), f"{label} must be an exact rational")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{label} must be an exact rational") from exc


def parse_vector(raw: Any, label: str) -> Vector:
    _require(isinstance(raw, list) and len(raw) == 2, f"{label} must have length two")
    return parse_fraction(raw[0], f"{label}[0]"), parse_fraction(raw[1], f"{label}[1]")


def outer(left: Vector, right: Vector) -> Grid:
    return tuple(
        tuple(left_value * right_value for right_value in right)
        for left_value in left
    )  # type: ignore[return-value]


def add_grids(left: Grid, right: Grid) -> Grid:
    return tuple(
        tuple(left[i][j] + right[i][j] for j in range(2))
        for i in range(2)
    )  # type: ignore[return-value]


def dot(left: Vector, right: Vector) -> Fraction:
    return left[0] * right[0] + left[1] * right[1]


def angular_filter(grid: Grid, row: Vector) -> Vector:
    """Apply the orientation row and retain the two-size vector."""

    return tuple(
        row[0] * grid[0][size] + row[1] * grid[1][size]
        for size in range(2)
    )  # type: ignore[return-value]


def radial_filter(grid: Grid, row: Vector) -> Vector:
    """Apply the size row and retain the two-orientation vector."""

    return tuple(dot(grid[orientation], row) for orientation in range(2))  # type: ignore[return-value]


def compose_angular_then_radial(grid: Grid, angular: Vector, radial: Vector) -> Fraction:
    return dot(angular_filter(grid, angular), radial)


def compose_radial_then_angular(grid: Grid, angular: Vector, radial: Vector) -> Fraction:
    return dot(radial_filter(grid, radial), angular)


def payload(value: Fraction | Sequence[Fraction]) -> str | list[str]:
    if isinstance(value, Fraction):
        return str(value)
    return [str(item) for item in value]


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 8, "issue must be 8")
    _require(contract.get("status") == "synthetic_exact_tensor_gate", "status drifted")

    angular = parse_vector(contract.get("angular_scalar_row"), "angular scalar row")
    radial = parse_vector(contract.get("radial_h4_row"), "radial H4 row")
    scalar_orientation = parse_vector(contract.get("scalar_orientation"), "scalar orientation")
    h4_orientation = parse_vector(contract.get("h4_orientation"), "H4 orientation")
    scalar_sizes = parse_vector(contract.get("scalar_sizes"), "scalar sizes")
    h4_sizes = parse_vector(contract.get("h4_sizes"), "H4 sizes")

    _require(dot(angular, scalar_orientation) == 1, "angular row must retain the scalar sector")
    _require(dot(angular, h4_orientation) == 0, "angular row did not kill H4")
    _require(dot(radial, h4_sizes) == 0, "radial row did not kill the declared H4 size law")
    _require(dot(radial, scalar_sizes) != 0, "radial row must distinguish the scalar size law")

    pure_h4 = outer(h4_orientation, h4_sizes)
    pure_scalar = outer(scalar_orientation, scalar_sizes)
    mixed = add_grids(pure_scalar, pure_h4)
    sectors = {"pure_h4": pure_h4, "pure_scalar": pure_scalar, "mixed": mixed}
    rows: dict[str, Any] = {}
    for name, grid in sectors.items():
        ar = compose_angular_then_radial(grid, angular, radial)
        ra = compose_radial_then_angular(grid, angular, radial)
        _require(ar == ra, f"filters did not commute on {name}")
        rows[name] = {
            "angular_then_radial": str(ar),
            "radial_then_angular": str(ra),
            "commutes": True,
            "angular_output": payload(angular_filter(grid, angular)),
            "radial_output": payload(radial_filter(grid, radial)),
        }

    _require(rows["pure_h4"]["angular_then_radial"] == "0", "pure H4 composition survived")
    _require(angular_filter(pure_h4, angular) == (Fraction(0), Fraction(0)), "angular H4 kernel drifted")
    _require(radial_filter(pure_h4, radial) == (Fraction(0), Fraction(0)), "radial H4 kernel drifted")
    _require(rows["mixed"]["angular_then_radial"] != "0", "mixed-sector control lost scalar residue")

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_synthetic_angular_radial_gate",
        "angular_and_radial_filters_commute": True,
        "pure_h4_killed_by_each_filter": True,
        "pure_h4_double_filter_is_redundant": True,
        "mixed_sector_retains_scalar_residue": True,
        "sectors": rows,
        "contains_production_data": False,
        "identifies_h4_exponent": False,
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
