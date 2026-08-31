#!/usr/bin/env python3
"""Exact cyclic deck-charge neutrality rules for invariant observables."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Sequence


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "cyclic_deck_charge_selection_contract.json"


def total_charge(order: int, charges: Sequence[int]) -> int:
    if order < 2:
        raise ValueError("cyclic order must be at least two")
    return sum(charges) % order


def invariant_allowed(order: int, charges: Sequence[int]) -> bool:
    return total_charge(order, charges) == 0


def translation_exponents(order: int, charges: Sequence[int]) -> List[int]:
    """Exponents acquired by a product character over every deck translate."""

    charge = total_charge(order, charges)
    return [(charge * shift) % order for shift in range(order)]


def residue_counts(order: int, charges: Sequence[int]) -> List[int]:
    exponents = translation_exponents(order, charges)
    return [exponents.count(residue) for residue in range(order)]


def orbit_cancels(order: int, charges: Sequence[int]) -> bool:
    """Return whether the exact root-of-unity orbit sum vanishes."""

    charge = total_charge(order, charges)
    if charge == 0:
        return False
    cycle_length = order // math.gcd(order, charge)
    return cycle_length > 1


def minimal_neutral_tensor_power(order: int, charge: int) -> int:
    charge %= order
    if charge == 0:
        return 1
    return order // math.gcd(order, charge)


def selection_matrix(order: int) -> List[List[bool]]:
    return [
        [invariant_allowed(order, (left, right)) for right in range(order)]
        for left in range(order)
    ]


def record(order: int, charges: Sequence[int]) -> Dict[str, object]:
    return {
        "charges": list(charges),
        "total_charge": total_charge(order, charges),
        "translation_exponents": translation_exponents(order, charges),
        "residue_counts": residue_counts(order, charges),
        "invariant_allowed": invariant_allowed(order, charges),
        "orbit_sum_zero": orbit_cancels(order, charges),
    }


def build_contract() -> Dict[str, object]:
    return {
        "schema": "matching-one/cyclic-deck-charge-selection/v1",
        "status": "valid_exact_cyclic_charge_neutrality_certificate",
        "parent_issue": "remain open",
        "rule": "an invariant scalar can contain a character product only if sum(charges)=0 mod q",
        "C2": {
            "selection_matrix": selection_matrix(2),
            "linear_detail": record(2, (1,)),
            "quadratic_detail": record(2, (1, 1)),
            "minimal_tensor_power_charge1": minimal_neutral_tensor_power(2, 1),
        },
        "C5": {
            "selection_matrix": selection_matrix(5),
            "linear_characters": [record(5, (charge,)) for charge in range(1, 5)],
            "conjugate_pairs": [record(5, pair) for pair in ((1, 4), (2, 3))],
            "same_charge_controls": [record(5, pair) for pair in ((1, 1), (2, 2))],
            "marked_observable_charge4_with_score_charge1": record(5, (4, 1)),
            "minimal_tensor_power": {
                str(charge): minimal_neutral_tensor_power(5, charge)
                for charge in range(1, 5)
            },
        },
        "asserts_nonzero_lattice_overlap": False,
        "uses_measured_responses": False,
    }


def validate_contract(path: Path = CONTRACT) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        frozen = json.load(handle)
    actual = build_contract()
    if frozen != actual:
        raise AssertionError("checked-in cyclic deck charge contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
