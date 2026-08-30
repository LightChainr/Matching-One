#!/usr/bin/env python3
"""Structural square-bond dual transport and centered-parity certificate.

For every checked periodic length, geometric dual transport is complement
followed by a permutation P of the 2 L^2 primal bonds.  P swaps horizontal and
vertical bonds, while P^2 is the common torus translation (-1, -1).  Hence the
transport is bijective, swaps primal and dual wrapping up to a homology-
preserving translation, and sends occupation k to B-k.

Consequently, any transport-even Bernstein aggregate has coefficients
a_k=a_(B-k), while a transport-odd aggregate has a_k=-a_(B-k).  Substitution
p=1/2+t makes the former an even polynomial and the latter an odd polynomial.
No configuration enumeration is needed for this structural certificate.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from square_bond_centered_parity import centered_from_bernstein_sums
from square_bond_duality_exact import square_bond_pairs


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "square_bond_transport_parity_contract.json"
EXPECTED_SCHEMA = "matching-one/square-bond-transport-parity/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def dual_bond_permutation(length: int) -> tuple[int, ...]:
    """Return P where the dual edge crossing bond i is primal bond P(i)."""

    pairs = square_bond_pairs(length)
    primal_index = {pair.primal: index for index, pair in enumerate(pairs)}
    _require(len(primal_index) == len(pairs), "primal bond registry is not injective")
    try:
        permutation = tuple(primal_index[pair.dual] for pair in pairs)
    except KeyError as exc:
        raise ValueError("a crossing dual edge is absent from the primal registry") from exc
    _require(len(set(permutation)) == len(pairs), "dual bond map is not a permutation")
    return permutation


def translated_bond_index(length: int, index: int, dx: int, dy: int) -> int:
    cell, direction = divmod(index, 2)
    y, x = divmod(cell, length)
    translated_cell = ((y + dy) % length) * length + ((x + dx) % length)
    return 2 * translated_cell + direction


def certify_length(length: int) -> dict[str, Any]:
    _require(length >= 2, "length must be at least two")
    permutation = dual_bond_permutation(length)
    bond_count = 2 * length * length
    squared = tuple(permutation[permutation[index]] for index in range(bond_count))
    translation = tuple(
        translated_bond_index(length, index, -1, -1)
        for index in range(bond_count)
    )
    swaps_orientation = all(index % 2 != target % 2 for index, target in enumerate(permutation))
    square_is_translation = squared == translation
    _require(swaps_orientation, f"L={length} dual map did not swap bond orientations")
    _require(square_is_translation, f"L={length} P^2 translation certificate failed")
    return {
        "length": length,
        "bond_count": bond_count,
        "permutation_is_bijective": True,
        "permutation_swaps_orientations": swaps_orientation,
        "permutation_square_translation": [-1, -1],
        "permutation_square_preserves_wrapping_channels": square_is_translation,
        "transport_flips_occupation": True,
    }


def reversal_fixture(degree: int, sign: int) -> tuple[Fraction, ...]:
    _require(degree >= 1, "degree must be positive")
    _require(sign in (-1, 1), "sign must be plus or minus one")
    values = [Fraction(0)] * (degree + 1)
    for occupied in range((degree // 2) + 1):
        reflected = degree - occupied
        value = Fraction((occupied + 1) ** 2, occupied + 2)
        if occupied == reflected and sign == -1:
            value = Fraction(0)
        values[occupied] = value
        values[reflected] = sign * value
    return tuple(values)


def certify_centered_parity(degree: int, sign: int) -> dict[str, Any]:
    coefficients = reversal_fixture(degree, sign)
    reversal = all(
        coefficients[degree - occupied] == sign * coefficients[occupied]
        for occupied in range(degree + 1)
    )
    centered = centered_from_bernstein_sums(coefficients)
    forbidden_parity = 1 if sign == 1 else 0
    forbidden_zero = all(
        value == 0
        for power, value in enumerate(centered)
        if power % 2 == forbidden_parity
    )
    _require(reversal, "coefficient reversal fixture drifted")
    _require(forbidden_zero, "coefficient reversal did not imply centered parity")
    return {
        "degree": degree,
        "transport_sign": sign,
        "coefficient_reversal": True,
        "centered_polynomial_parity": "even" if sign == 1 else "odd",
        "forbidden_coefficients_zero": True,
    }


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 42, "issue must be 42")
    _require(contract.get("status") == "exact_transport_structure_only", "status drifted")
    lengths = contract.get("lengths")
    _require(
        isinstance(lengths, list)
        and lengths
        and all(isinstance(length, int) for length in lengths),
        "lengths must be a nonempty integer list",
    )
    _require(len(set(lengths)) == len(lengths), "lengths must be distinct")
    _require(lengths == sorted(lengths), "lengths must be sorted")
    _require(contract.get("permutation_square_translation") == [-1, -1], "translation drifted")

    structures = [certify_length(length) for length in lengths]
    parity_rows = [
        certify_centered_parity(row["bond_count"], sign)
        for row in structures
        for sign in (1, -1)
    ]
    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_square_bond_transport_parity_certificate",
        "lengths_checked": lengths,
        "structures": structures,
        "transport_is_complement_plus_permutation": True,
        "occupation_coefficient_reversal_certified": True,
        "transport_even_implies_centered_even": True,
        "transport_odd_implies_centered_odd": True,
        "parity_checks": parity_rows,
        "enumerates_configurations": False,
        "contains_continuum_amplitude_claim": False,
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
