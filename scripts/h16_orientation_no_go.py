#!/usr/bin/env python3
"""Exact orientation-count obstruction to adding H16 at fixed N=1105."""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

try:
    from scripts.minimal_four_orientation_gaussian_torus import (
        primitive_first_octant_representations,
        primitive_orientation_count_formula,
    )
except ModuleNotFoundError:  # Direct execution from scripts/.
    from minimal_four_orientation_gaussian_torus import (  # type: ignore
        primitive_first_octant_representations,
        primitive_orientation_count_formula,
    )


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "h16_orientation_no_go_contract.json"
N1105 = 5 * 13 * 17
FIRST_FIVE_ORIENTATION_N = 5 * 13 * 17 * 29


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor = 3 if divisor == 2 else divisor + 2
    return True


def first_split_primes(count: int) -> List[int]:
    if count < 0:
        raise ValueError("count must be nonnegative")
    out: List[int] = []
    candidate = 5
    while len(out) < count:
        if candidate % 4 == 1 and is_prime(candidate):
            out.append(candidate)
        candidate += 4
    return out


def minimum_primitive_layer(required_orbits: int) -> Tuple[int, int, List[int]]:
    """Return the least N that can have at least the requested orbit count."""

    if required_orbits < 1:
        raise ValueError("required_orbits must be positive")
    split_count = 0
    capacity = 1
    while capacity < required_orbits:
        split_count += 1
        capacity = 1 << max(0, split_count - 1)
    if split_count == 0:
        return 1, 1, []
    primes = first_split_primes(split_count)
    n = 1
    for prime in primes:
        n *= prime
    return n, capacity, primes


def gaussian_mul(left: Tuple[int, int], right: Tuple[int, int]) -> Tuple[int, int]:
    a, b = left
    c, d = right
    return a * c - b * d, a * d + b * c


def gaussian_pow(value: Tuple[int, int], exponent: int) -> Tuple[int, int]:
    result = (1, 0)
    base = value
    while exponent:
        if exponent & 1:
            result = gaussian_mul(result, base)
        base = gaussian_mul(base, base)
        exponent >>= 1
    return result


def harmonic_value(n: int, orientation: Tuple[int, int], index: int) -> Fraction:
    if index == 0:
        return Fraction(1)
    real, _imaginary = gaussian_pow(orientation, 4 * index)
    return Fraction(real, n ** (2 * index))


def matrix_rank(matrix: Sequence[Sequence[Fraction]]) -> int:
    rows = [list(map(Fraction, row)) for row in matrix]
    if not rows:
        return 0
    row = 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(row, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[row], rows[pivot] = rows[pivot], rows[row]
        scale = rows[row][column]
        rows[row] = [value / scale for value in rows[row]]
        for i in range(len(rows)):
            if i != row and rows[i][column]:
                scale = rows[i][column]
                rows[i] = [x - scale * y for x, y in zip(rows[i], rows[row])]
        row += 1
        if row == len(rows):
            break
    return row


def harmonic_rank(n: int, orientations: Sequence[Tuple[int, int]]) -> int:
    matrix = [
        [harmonic_value(n, orientation, index) for index in range(5)]
        for orientation in orientations
    ]
    return matrix_rank(matrix)


def first_full_rank_subset(n: int) -> Tuple[Tuple[int, int], ...]:
    orientations = primitive_first_octant_representations(n)
    for subset in itertools.combinations(orientations, 5):
        if harmonic_rank(n, subset) == 5:
            return subset
    raise AssertionError("no five-orientation subset resolves H0 through H16")


def pairs(values: Sequence[Tuple[int, int]]) -> List[List[int]]:
    return [[a, b] for a, b in values]


def build_contract() -> Dict[str, object]:
    n4, capacity4, primes4 = minimum_primitive_layer(4)
    n5, capacity5, primes5 = minimum_primitive_layer(5)
    reps1105 = primitive_first_octant_representations(N1105)
    reps_first = primitive_first_octant_representations(n5)
    subset = first_full_rank_subset(n5)
    return {
        "schema": "matching-one/h16-orientation-no-go/v1",
        "status": "valid_exact_orientation_count_obstruction",
        "parent_issue": "remain open",
        "n1105": {
            "n": N1105,
            "factorization": primes4,
            "primitive_d4_orbit_count": len(reps1105),
            "formula_count": primitive_orientation_count_formula(N1105),
            "orientations": pairs(reps1105),
            "h0_through_h16_rank_upper_bound": len(reps1105),
            "can_resolve_five_generic_harmonics": False,
        },
        "first_layer_with_at_least_five_orbits": {
            "n": n5,
            "factorization": primes5,
            "available_orbits": capacity5,
            "orientations": pairs(reps_first),
            "first_full_rank_subset": pairs(subset),
            "h0_through_h16_rank": harmonic_rank(n5, subset),
        },
        "minimum_four_orbit_layer": n4,
        "minimum_four_orbit_capacity": capacity4,
        "recommends_a_production_size": False,
    }


def validate_contract(path: Path = CONTRACT) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        frozen = json.load(handle)
    actual = build_contract()
    if frozen != actual:
        raise AssertionError("checked-in H16 orientation no-go contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
