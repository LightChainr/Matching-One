#!/usr/bin/env python3
"""Exact conjugation and norm relations for the formal Gaussian Jordan transfer.

This extends the existing multiplier-semigroup oracle without changing it.  For
Gaussian conjugation z -> zbar, the quartic character and every formal matrix
coefficient conjugate entrywise.  Multiplying z by zbar gives the positive
integer norm, so multiplicativity requires

    T(zbar) = conjugate(T(z)),
    T(z) T(zbar) = T(N(z)).

All checks use exact Gaussian integers, Fractions, and formal logarithm symbols.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

from formal_jordan_multiplier_semigroup import (
    ComplexFraction,
    LinearFormalLog,
    Matrix,
    gaussian_multiply,
    gaussian_norm,
    gaussian_pair,
    jordan_transfer,
    matrix_multiply,
    matrix_payload,
    parse_expected_character,
    quartic_character,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "formal_jordan_conjugate_semigroup_contract.json"
EXPECTED_SCHEMA = "matching-one/formal-jordan-conjugate-semigroup/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def gaussian_conjugate(value: tuple[int, int]) -> tuple[int, int]:
    return value[0], -value[1]


def complex_conjugate(value: ComplexFraction) -> ComplexFraction:
    return ComplexFraction(value.real, -value.imag)


def formal_entry_conjugate(value: LinearFormalLog) -> LinearFormalLog:
    return LinearFormalLog(
        constant=complex_conjugate(value.constant),
        log2=complex_conjugate(value.log2),
        log5=complex_conjugate(value.log5),
    )


def matrix_conjugate(matrix: Matrix) -> Matrix:
    return tuple(
        tuple(formal_entry_conjugate(entry) for entry in row)
        for row in matrix
    )  # type: ignore[return-value]


def certify_multiplier(value: tuple[int, int]) -> dict[str, Any]:
    conjugate = gaussian_conjugate(value)
    norm = gaussian_norm(value)
    norm_scalar = (norm, 0)
    _require(gaussian_conjugate(conjugate) == value, "Gaussian conjugation is not involutive")
    _require(
        gaussian_multiply(value, conjugate) == norm_scalar,
        "Gaussian norm factorization failed",
    )

    character = quartic_character(value)
    conjugate_character = quartic_character(conjugate)
    _require(
        conjugate_character == complex_conjugate(character),
        "quartic character did not intertwine conjugation",
    )

    transfer = jordan_transfer(value)
    conjugate_transfer = jordan_transfer(conjugate)
    _require(
        conjugate_transfer == matrix_conjugate(transfer),
        "formal Jordan transfer did not intertwine conjugation",
    )
    norm_transfer = jordan_transfer(norm_scalar)
    forward = matrix_multiply(transfer, conjugate_transfer)
    reverse = matrix_multiply(conjugate_transfer, transfer)
    _require(
        forward == reverse == norm_transfer,
        "formal Jordan norm-factorization paths did not close",
    )
    return {
        "multiplier": list(value),
        "conjugate": list(conjugate),
        "norm": norm,
        "character": character.payload(),
        "conjugate_character": conjugate_character.payload(),
        "conjugation_intertwines_character": True,
        "conjugation_intertwines_transfer": True,
        "norm_factorization_closes_both_orders": True,
        "norm_transfer": matrix_payload(norm_transfer),
    }


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 145, "issue must be 145")
    _require(contract.get("status") == "exact_formal_conjugation_only", "status drifted")
    raw_multipliers = contract.get("multipliers")
    _require(isinstance(raw_multipliers, Mapping), "multipliers must be a mapping")
    multipliers = {
        name: gaussian_pair(raw_multipliers.get(name), name)
        for name in ("norm2", "norm5", "norm10")
    }
    _require(gaussian_norm(multipliers["norm2"]) == 2, "norm2 multiplier drifted")
    _require(gaussian_norm(multipliers["norm5"]) == 5, "norm5 multiplier drifted")
    _require(gaussian_norm(multipliers["norm10"]) == 10, "norm10 multiplier drifted")
    _require(
        gaussian_multiply(multipliers["norm2"], multipliers["norm5"])
        == multipliers["norm10"],
        "declared multiplier product did not close",
    )
    _require(
        gaussian_multiply(
            gaussian_conjugate(multipliers["norm2"]),
            gaussian_conjugate(multipliers["norm5"]),
        )
        == gaussian_conjugate(multipliers["norm10"]),
        "conjugated multiplier product did not close",
    )

    expected = contract.get("expected_quartic_characters")
    _require(isinstance(expected, Mapping), "expected characters must be a mapping")
    certificates: dict[str, Any] = {}
    for name, multiplier in multipliers.items():
        declared = parse_expected_character(expected.get(name), f"{name} character")
        _require(quartic_character(multiplier) == declared, f"{name} quartic character drifted")
        certificates[name] = certify_multiplier(multiplier)

    combined_transfer = jordan_transfer(multipliers["norm10"])
    product_transfer = matrix_multiply(
        jordan_transfer(multipliers["norm2"]),
        jordan_transfer(multipliers["norm5"]),
    )
    _require(combined_transfer == product_transfer, "direct and composite transfers drifted")
    _require(
        matrix_conjugate(product_transfer)
        == matrix_multiply(
            jordan_transfer(gaussian_conjugate(multipliers["norm2"])),
            jordan_transfer(gaussian_conjugate(multipliers["norm5"])),
        ),
        "conjugation did not preserve the declared composite path",
    )

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_formal_jordan_conjugate_semigroup",
        "certificates": certificates,
        "conjugation_involutive": True,
        "quartic_character_star_compatible": True,
        "jordan_transfer_star_compatible": True,
        "norm_factorization_exact": True,
        "composite_path_star_compatible": True,
        "constructs_group_inverse": False,
        "contains_physical_jordan_claim": False,
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
