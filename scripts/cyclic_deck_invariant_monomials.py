#!/usr/bin/env python3
"""Exact bounded monomial invariants of nontrivial cyclic deck characters."""

from __future__ import annotations

import argparse
from itertools import product
import json
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "cyclic_deck_invariant_monomial_contract.json"
Exponent = Tuple[int, ...]


def total_charge(order: int, exponent: Sequence[int]) -> int:
    if order < 2 or len(exponent) != order - 1:
        raise ValueError("exponent vector must contain every nontrivial cyclic charge")
    if any(value < 0 for value in exponent):
        raise ValueError("monomial exponents must be nonnegative")
    return sum(charge * exponent[charge - 1] for charge in range(1, order)) % order


def exponent_vectors(variable_count: int, degree: int) -> List[Exponent]:
    if variable_count < 1 or degree < 0:
        raise ValueError("variable count must be positive and degree nonnegative")
    return [
        tuple(values)
        for values in product(range(degree + 1), repeat=variable_count)
        if sum(values) == degree
    ]


def neutral_monomials(order: int, degree: int) -> List[Exponent]:
    return [
        exponent
        for exponent in exponent_vectors(order - 1, degree)
        if total_charge(order, exponent) == 0
    ]


def is_primitive_neutral(order: int, exponent: Exponent) -> bool:
    degree = sum(exponent)
    if degree == 0 or total_charge(order, exponent) != 0:
        return False
    for left_degree in range(1, degree):
        for left in neutral_monomials(order, left_degree):
            if all(left[index] <= exponent[index] for index in range(order - 1)):
                right = tuple(exponent[index] - left[index] for index in range(order - 1))
                if total_charge(order, right) == 0:
                    return False
    return True


def primitive_neutral_generators(order: int, max_degree: int) -> List[Exponent]:
    return [
        exponent
        for degree in range(1, max_degree + 1)
        for exponent in neutral_monomials(order, degree)
        if is_primitive_neutral(order, exponent)
    ]


def invariant_census(order: int, max_degree: int) -> Dict[str, object]:
    by_degree = {
        str(degree): [list(exponent) for exponent in neutral_monomials(order, degree)]
        for degree in range(max_degree + 1)
    }
    return {
        "order": order,
        "variables": [f"z{charge}" for charge in range(1, order)],
        "variable_charges": list(range(1, order)),
        "max_degree": max_degree,
        "hilbert_counts": [len(by_degree[str(degree)]) for degree in range(max_degree + 1)],
        "neutral_exponents_by_degree": by_degree,
        "primitive_neutral_generators": [
            list(exponent) for exponent in primitive_neutral_generators(order, max_degree)
        ],
    }


def build_contract() -> Dict[str, object]:
    return {
        "schema": "matching-one/cyclic-deck-invariant-monomials/v1",
        "status": "valid_exact_bounded_invariant_monomial_census",
        "parent_issue": "remain open",
        "neutrality_rule": "sum_(charge=1)^(Q-1) charge*exponent_charge is zero modulo Q",
        "degree_bound": 5,
        "c2": invariant_census(2, 5),
        "c5": invariant_census(5, 5),
        "claim_boundary": {
            "proves_nonzero_lattice_overlap": False,
            "computes_transfer_amplitude": False,
            "reads_measured_response": False,
            "identifies_operator": False,
            "makes_production_recommendation": False,
        },
    }


def validate_contract(path: Path = CONTRACT) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        frozen = json.load(handle)
    actual = build_contract()
    if frozen != actual:
        raise AssertionError("checked-in cyclic invariant monomial contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
