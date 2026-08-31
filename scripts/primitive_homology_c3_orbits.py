#!/usr/bin/env python3
"""Exact C3 orbits and character projectors for primitive hexagonal lines."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "primitive_homology_c3_orbit_contract.json"
Vector = Tuple[int, int]
ROTATION = ((1, -1), (1, 0))


def determinant(matrix: Sequence[Sequence[int]]) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def apply(matrix: Sequence[Sequence[int]], vector: Vector) -> Vector:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )


def canonical_line(vector: Vector) -> Vector:
    a, b = vector
    if (a, b) == (0, 0):
        raise ValueError("the zero vector does not define a primitive line")
    divisor = math.gcd(abs(a), abs(b))
    if divisor != 1:
        raise ValueError("line representative must be primitive")
    if a < 0 or (a == 0 and b < 0):
        return -a, -b
    return a, b


def hexagonal_norm(vector: Vector) -> int:
    """Norm in the basis (1, omega), omega=exp(2*pi*i/3)."""

    a, b = vector
    return a * a - a * b + b * b


def rotate_line(vector: Vector) -> Vector:
    return canonical_line(apply(ROTATION, canonical_line(vector)))


def primitive_lines(norm_bound: int) -> List[Vector]:
    if norm_bound < 1:
        raise ValueError("norm bound must be positive")
    extent = 2 * math.isqrt(norm_bound) + 2
    lines = set()
    for a in range(-extent, extent + 1):
        for b in range(-extent, extent + 1):
            if (a, b) == (0, 0) or math.gcd(abs(a), abs(b)) != 1:
                continue
            if hexagonal_norm((a, b)) <= norm_bound:
                lines.add(canonical_line((a, b)))
    return sorted(lines)


def c3_orbits(norm_bound: int) -> List[Tuple[Vector, Vector, Vector]]:
    lines = set(primitive_lines(norm_bound))
    orbits: List[Tuple[Vector, Vector, Vector]] = []
    while lines:
        seed = min(lines)
        cycle = [seed, rotate_line(seed), rotate_line(rotate_line(seed))]
        if rotate_line(cycle[-1]) != seed or len(set(cycle)) != 3:
            raise AssertionError("primitive unoriented line did not have a C3 orbit")
        if not set(cycle).issubset(lines):
            raise AssertionError("hexagonal norm cutoff was not rotation invariant")
        start = min(range(3), key=lambda index: cycle[index])
        ordered = tuple(cycle[(start + offset) % 3] for offset in range(3))
        orbits.append(ordered)  # type: ignore[arg-type]
        lines.difference_update(cycle)
    return sorted(orbits)


def spin_charge(spin: int) -> int:
    """Return the C3 charge of an even spin under a 60-degree rotation."""

    if spin % 2:
        raise ValueError("unoriented-line characters require even spin")
    return (spin // 2) % 3


def projector_exponents(charge: int) -> List[int]:
    """Exponents of zeta3 in sum_j zeta3^(-charge*j) P_j."""

    charge %= 3
    return [(-charge * index) % 3 for index in range(3)]


def character_inner(left_charge: int, right_charge: int) -> int:
    """Exact C3 character inner product, avoiding floating complex roots."""

    return 3 if (left_charge - right_charge) % 3 == 0 else 0


def vectors(values: Sequence[Vector]) -> List[List[int]]:
    return [[a, b] for a, b in values]


def build_contract(norm_bound: int = 13) -> Dict[str, object]:
    lines = primitive_lines(norm_bound)
    orbits = c3_orbits(norm_bound)
    if determinant(ROTATION) != 1:
        raise AssertionError("rotation must be unimodular")
    sample = (2, 1)
    if apply(ROTATION, apply(ROTATION, apply(ROTATION, sample))) != (-2, -1):
        raise AssertionError("rotation cube must be minus the identity")
    return {
        "schema": "matching-one/primitive-homology-c3-orbits/v1",
        "status": "valid_exact_primitive_line_character_certificate",
        "parent_issue": "remain open",
        "basis": "(1, omega), omega=exp(2*pi*i/3)",
        "rotation_matrix": [list(row) for row in ROTATION],
        "rotation_determinant": determinant(ROTATION),
        "rotation_cube_on_oriented_vectors": "-I",
        "rotation_order_on_unoriented_lines": 3,
        "norm": "a^2-a*b+b^2",
        "norm_bound": norm_bound,
        "primitive_line_count": len(lines),
        "orbit_count": len(orbits),
        "orbits": [vectors(orbit) for orbit in orbits],
        "character_gram": [
            [character_inner(left, right) for right in range(3)]
            for left in range(3)
        ],
        "spin_characters": {
            f"H{spin}": {
                "charge": spin_charge(spin),
                "projector_exponents": projector_exponents(spin_charge(spin)),
            }
            for spin in (4, 8, 12)
        },
        "uses_continuum_baseline_or_production_data": False,
    }


def validate_contract(path: Path = CONTRACT) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        frozen = json.load(handle)
    actual = build_contract()
    if frozen != actual:
        raise AssertionError("checked-in primitive C3 orbit contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
