#!/usr/bin/env python3
"""Exact C6 orbits and quotient characters for oriented hexagonal homology."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

try:
    from scripts.primitive_homology_c3_orbits import (
        ROTATION,
        apply,
        determinant,
        primitive_lines,
    )
except ModuleNotFoundError:  # Direct `python3 scripts/...` execution.
    from primitive_homology_c3_orbits import (  # type: ignore[no-redef]
        ROTATION,
        apply,
        determinant,
        primitive_lines,
    )


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "oriented_homology_c6_orbit_contract.json"
Vector = Tuple[int, int]


def negate(vector: Vector) -> Vector:
    return -vector[0], -vector[1]


def oriented_primitive_vectors(norm_bound: int) -> List[Vector]:
    vectors = set()
    for line in primitive_lines(norm_bound):
        vectors.add(line)
        vectors.add(negate(line))
    return sorted(vectors)


def rotate(vector: Vector) -> Vector:
    return apply(ROTATION, vector)


def c6_orbits(norm_bound: int) -> List[Tuple[Vector, ...]]:
    remaining = set(oriented_primitive_vectors(norm_bound))
    orbits: List[Tuple[Vector, ...]] = []
    while remaining:
        seed = min(remaining)
        cycle = [seed]
        for _ in range(5):
            cycle.append(rotate(cycle[-1]))
        if rotate(cycle[-1]) != seed or len(set(cycle)) != 6:
            raise AssertionError("oriented primitive vector did not have a C6 orbit")
        if not set(cycle).issubset(remaining):
            raise AssertionError("norm cutoff was not C6 invariant")
        start = min(range(6), key=lambda index: cycle[index])
        ordered = tuple(cycle[(start + offset) % 6] for offset in range(6))
        orbits.append(ordered)
        remaining.difference_update(cycle)
    return sorted(orbits)


def c6_character_inner(left_charge: int, right_charge: int) -> int:
    return 6 if (left_charge - right_charge) % 6 == 0 else 0


def oriented_spin_charge(spin: int) -> int:
    return spin % 6


def descends_to_unoriented_lines(charge: int) -> bool:
    """The central half-turn R^3 acts trivially exactly for even charge."""

    return charge % 2 == 0


def quotient_c3_charge(charge: int) -> int:
    if not descends_to_unoriented_lines(charge):
        raise ValueError("odd C6 charge does not descend through v~-v")
    return (charge // 2) % 3


def vector_lists(values: Sequence[Vector]) -> List[List[int]]:
    return [[a, b] for a, b in values]


def build_contract(norm_bound: int = 13) -> Dict[str, object]:
    vectors = oriented_primitive_vectors(norm_bound)
    orbits = c6_orbits(norm_bound)
    probe = (2, 1)
    rotated = probe
    for _ in range(3):
        rotated = rotate(rotated)
    if rotated != negate(probe):
        raise AssertionError("rotation cube must be the central sign involution")
    return {
        "schema": "matching-one/oriented-homology-c6-orbits/v1",
        "status": "valid_exact_oriented_c6_character_certificate",
        "parent_issue": "remain open",
        "basis": "(1, omega), omega=exp(2*pi*i/3)",
        "rotation_matrix": [list(row) for row in ROTATION],
        "rotation_determinant": determinant(ROTATION),
        "rotation_cube": "-I",
        "rotation_sixth_power": "I",
        "norm": "a^2-a*b+b^2",
        "norm_bound": norm_bound,
        "oriented_primitive_vector_count": len(vectors),
        "orbit_count": len(orbits),
        "orbits": [vector_lists(orbit) for orbit in orbits],
        "character_gram": [
            [c6_character_inner(left, right) for right in range(6)]
            for left in range(6)
        ],
        "quotient_rule": "C6 charge descends through v~-v iff it is even; quotient charge is charge/2 modulo 3",
        "spin_characters": {
            f"H{spin}": {
                "c6_charge": oriented_spin_charge(spin),
                "descends": descends_to_unoriented_lines(oriented_spin_charge(spin)),
                "c3_quotient_charge": quotient_c3_charge(oriented_spin_charge(spin)),
            }
            for spin in (4, 8, 12)
        },
        "odd_spin_controls": {
            f"H{spin}": {
                "c6_charge": oriented_spin_charge(spin),
                "descends": descends_to_unoriented_lines(oriented_spin_charge(spin)),
            }
            for spin in (1, 3, 5)
        },
        "uses_continuum_baseline_or_production_data": False,
    }


def validate_contract(path: Path = CONTRACT) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        frozen = json.load(handle)
    actual = build_contract()
    if frozen != actual:
        raise AssertionError("checked-in oriented C6 orbit contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
