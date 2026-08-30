#!/usr/bin/env python3
"""Exact formal rank-2 Jordan representation of a Gaussian multiplier semigroup."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping

from exact_jet_algebra import canonical_fraction


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "formal_jordan_multiplier_semigroup_contract.json"
EXPECTED_SCHEMA = "matching-one/formal-jordan-multiplier-semigroup/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


@dataclass(frozen=True)
class ComplexFraction:
    real: Fraction = Fraction(0)
    imag: Fraction = Fraction(0)

    def __add__(self, other: "ComplexFraction") -> "ComplexFraction":
        return ComplexFraction(self.real + other.real, self.imag + other.imag)

    def __neg__(self) -> "ComplexFraction":
        return ComplexFraction(-self.real, -self.imag)

    def __sub__(self, other: "ComplexFraction") -> "ComplexFraction":
        return self + (-other)

    def __mul__(self, other: "ComplexFraction") -> "ComplexFraction":
        return ComplexFraction(
            self.real * other.real - self.imag * other.imag,
            self.real * other.imag + self.imag * other.real,
        )

    def inverse(self) -> "ComplexFraction":
        denominator = self.real**2 + self.imag**2
        _require(denominator != 0, "zero Gaussian multiplier has no inverse")
        return ComplexFraction(self.real / denominator, -self.imag / denominator)

    def power(self, exponent: int) -> "ComplexFraction":
        _require(exponent >= 0, "complex exponent must be nonnegative")
        result = ComplexFraction(Fraction(1), Fraction(0))
        for _ in range(exponent):
            result = result * self
        return result

    def payload(self) -> list[str]:
        return [str(self.real), str(self.imag)]


ZERO_COMPLEX = ComplexFraction()
ONE_COMPLEX = ComplexFraction(Fraction(1), Fraction(0))


@dataclass(frozen=True)
class LinearFormalLog:
    """Complex coefficient linear in the independent symbols log(2), log(5)."""

    constant: ComplexFraction = ZERO_COMPLEX
    log2: ComplexFraction = ZERO_COMPLEX
    log5: ComplexFraction = ZERO_COMPLEX

    def __add__(self, other: "LinearFormalLog") -> "LinearFormalLog":
        return LinearFormalLog(
            self.constant + other.constant,
            self.log2 + other.log2,
            self.log5 + other.log5,
        )

    def __mul__(self, other: "LinearFormalLog") -> "LinearFormalLog":
        left_has_log = self.log2 != ZERO_COMPLEX or self.log5 != ZERO_COMPLEX
        right_has_log = other.log2 != ZERO_COMPLEX or other.log5 != ZERO_COMPLEX
        _require(not (left_has_log and right_has_log), "formal product exceeded log degree one")
        return LinearFormalLog(
            self.constant * other.constant,
            self.constant * other.log2 + self.log2 * other.constant,
            self.constant * other.log5 + self.log5 * other.constant,
        )

    def payload(self) -> dict[str, list[str]]:
        return {
            "constant": self.constant.payload(),
            "log2": self.log2.payload(),
            "log5": self.log5.payload(),
        }


ZERO_ENTRY = LinearFormalLog()
Matrix = tuple[
    tuple[LinearFormalLog, LinearFormalLog],
    tuple[LinearFormalLog, LinearFormalLog],
]


def constant_entry(value: ComplexFraction) -> LinearFormalLog:
    return LinearFormalLog(constant=value)


def matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            left[i][0] * right[0][j] + left[i][1] * right[1][j]
            for j in range(2)
        )
        for i in range(2)
    )  # type: ignore[return-value]


def matrix_payload(matrix: Matrix) -> list[list[dict[str, list[str]]]]:
    return [[entry.payload() for entry in row] for row in matrix]


def gaussian_pair(raw: Any, label: str) -> tuple[int, int]:
    _require(
        isinstance(raw, list)
        and len(raw) == 2
        and all(isinstance(value, int) for value in raw),
        f"{label} must be an integer pair",
    )
    pair = raw[0], raw[1]
    _require(pair != (0, 0), f"{label} must be nonzero")
    return pair


def gaussian_multiply(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    a, b = left
    c, d = right
    return a * c - b * d, a * d + b * c


def gaussian_norm(value: tuple[int, int]) -> int:
    return value[0] ** 2 + value[1] ** 2


def quartic_character(value: tuple[int, int]) -> ComplexFraction:
    gaussian = ComplexFraction(Fraction(value[0]), Fraction(value[1]))
    return gaussian.inverse().power(4)


def formal_log_vector(norm: int) -> tuple[int, int]:
    _require(norm > 0, "norm must be positive")
    exponents = []
    remainder = norm
    for prime in (2, 5):
        exponent = 0
        while remainder % prime == 0:
            exponent += 1
            remainder //= prime
        exponents.append(exponent)
    _require(remainder == 1, "norm contains an undeclared formal logarithm prime")
    return exponents[0], exponents[1]


def jordan_transfer(value: tuple[int, int]) -> Matrix:
    character = quartic_character(value)
    log2_power, log5_power = formal_log_vector(gaussian_norm(value))
    upper = LinearFormalLog(
        log2=ComplexFraction(character.real * log2_power, character.imag * log2_power),
        log5=ComplexFraction(character.real * log5_power, character.imag * log5_power),
    )
    diagonal = constant_entry(character)
    return ((diagonal, upper), (ZERO_ENTRY, diagonal))


def parse_expected_character(raw: Any, label: str) -> ComplexFraction:
    _require(isinstance(raw, list) and len(raw) == 2, f"{label} must be a pair")
    return ComplexFraction(
        canonical_fraction(raw[0], f"{label} real"),
        canonical_fraction(raw[1], f"{label} imaginary"),
    )


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 145, "issue must be 145")
    _require(contract.get("status") == "exact_formal_rank2_only", "status drifted")
    raw_multipliers = contract.get("multipliers")
    _require(isinstance(raw_multipliers, Mapping), "multipliers must be a mapping")
    multipliers = {
        name: gaussian_pair(raw_multipliers.get(name), name)
        for name in ("norm2", "norm5", "norm10")
    }
    _require(gaussian_norm(multipliers["norm2"]) == 2, "norm2 multiplier drifted")
    _require(gaussian_norm(multipliers["norm5"]) == 5, "norm5 multiplier drifted")
    _require(gaussian_norm(multipliers["norm10"]) == 10, "norm10 multiplier drifted")
    product_25 = gaussian_multiply(multipliers["norm2"], multipliers["norm5"])
    product_52 = gaussian_multiply(multipliers["norm5"], multipliers["norm2"])
    _require(product_25 == product_52 == multipliers["norm10"], "Gaussian multiplier paths do not close")

    raw_characters = contract.get("expected_quartic_characters")
    raw_logs = contract.get("expected_formal_log_vectors")
    _require(isinstance(raw_characters, Mapping), "expected characters must be a mapping")
    _require(isinstance(raw_logs, Mapping), "expected log vectors must be a mapping")
    characters: dict[str, ComplexFraction] = {}
    log_vectors: dict[str, tuple[int, int]] = {}
    for name, multiplier in multipliers.items():
        character = quartic_character(multiplier)
        expected = parse_expected_character(raw_characters.get(name), f"{name} character")
        _require(character == expected, f"{name} quartic character drifted")
        characters[name] = character
        vector = formal_log_vector(gaussian_norm(multiplier))
        raw_vector = raw_logs.get(name)
        _require(isinstance(raw_vector, list) and tuple(raw_vector) == vector, f"{name} formal log vector drifted")
        log_vectors[name] = vector

    _require(
        characters["norm2"] * characters["norm5"] == characters["norm10"],
        "quartic character did not compose",
    )
    direct = jordan_transfer(multipliers["norm10"])
    path_25 = matrix_multiply(
        jordan_transfer(multipliers["norm2"]),
        jordan_transfer(multipliers["norm5"]),
    )
    path_52 = matrix_multiply(
        jordan_transfer(multipliers["norm5"]),
        jordan_transfer(multipliers["norm2"]),
    )
    _require(direct == path_25 == path_52, "formal Jordan transfer paths do not close")
    _require(direct[0][1].log2 != ZERO_COMPLEX, "norm10 log2 cocycle vanished")
    _require(direct[0][1].log5 != ZERO_COMPLEX, "norm10 log5 cocycle vanished")

    nilpotent: Matrix = (
        (ZERO_ENTRY, constant_entry(ONE_COMPLEX)),
        (ZERO_ENTRY, ZERO_ENTRY),
    )
    zero_matrix: Matrix = ((ZERO_ENTRY, ZERO_ENTRY), (ZERO_ENTRY, ZERO_ENTRY))
    _require(matrix_multiply(nilpotent, nilpotent) == zero_matrix, "declared generator is not square-zero")

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_exact_formal_jordan_multiplier_semigroup",
        "gaussian_products_close_both_orders": True,
        "quartic_characters": {name: value.payload() for name, value in characters.items()},
        "quartic_character_multiplicative": True,
        "formal_log_vectors": {name: list(value) for name, value in log_vectors.items()},
        "common_nilpotent_squares_to_zero": True,
        "direct_norm10_transfer": matrix_payload(direct),
        "norm2_norm5_matrix_paths_close": True,
        "contains_physical_jordan_claim": False,
        "contains_target_data": False,
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
