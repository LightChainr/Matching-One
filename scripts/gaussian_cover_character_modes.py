#!/usr/bin/env python3
"""Exact deck-group and character oracle for selected Gaussian covers.

This is deliberately a small rank-two integer calculation.  Character values
are stored as exact exponents of roots of unity; no floating point or Monte
Carlo input enters the oracle.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import reduce
import json
from math import gcd
from pathlib import Path
from typing import Iterable


Vector = tuple[int, int]
Matrix = tuple[tuple[int, int], tuple[int, int]]

MULTIPLIERS: tuple[tuple[str, Vector, bool], ...] = (
    ("1+i", (1, 1), True),
    ("2+i", (2, 1), True),
    ("2-i", (2, -1), True),
    ("2i", (0, 2), True),
    ("3+i", (3, 1), True),
    ("3-i", (3, -1), False),  # conjugation target required by D4
)

IDENTITY: Matrix = ((1, 0), (0, 1))
ROTATION_90: Matrix = ((0, -1), (1, 0))
CONJUGATION: Matrix = ((1, 0), (0, -1))


def matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def matrix_vector(matrix: Matrix, vector: Vector) -> Vector:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )


def transpose(matrix: Matrix) -> Matrix:
    return ((matrix[0][0], matrix[1][0]), (matrix[0][1], matrix[1][1]))


def determinant(matrix: Matrix) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def multiplication_matrix(multiplier: Vector) -> Matrix:
    a, b = multiplier
    return ((a, -b), (b, a))


def smith_invariants(matrix: Matrix) -> tuple[int, int]:
    """Smith invariants of a nonsingular 2x2 integer matrix."""

    det = abs(determinant(matrix))
    if det == 0:
        raise ValueError("cover multiplier must be nonzero")
    d1 = reduce(gcd, (abs(value) for row in matrix for value in row))
    d2 = det // d1
    if d2 % d1:
        raise ArithmeticError("invalid Smith divisibility")
    return d1, d2


def adjugate(matrix: Matrix) -> Matrix:
    return (
        (matrix[1][1], -matrix[0][1]),
        (-matrix[1][0], matrix[0][0]),
    )


def coset_key(matrix: Matrix, vector: Vector) -> Vector:
    """Exact key for vector modulo the column lattice of matrix."""

    det = abs(determinant(matrix))
    lifted = matrix_vector(adjugate(matrix), vector)
    return lifted[0] % det, lifted[1] % det


def _candidate_key(vector: Vector) -> tuple[int, int, int, int, int]:
    x, y = vector
    return (
        x * x + y * y,
        abs(x) + abs(y),
        int(x < 0 or y < 0),
        -x,
        -y,
    )


def canonical_representatives(matrix: Matrix) -> list[Vector]:
    """Choose deterministic short representatives for every lattice coset."""

    order = abs(determinant(matrix))
    representatives: dict[Vector, Vector] = {}
    radius = 0
    while len(representatives) < order:
        shell = [
            (x, y)
            for x in range(-radius, radius + 1)
            for y in range(-radius, radius + 1)
            if max(abs(x), abs(y)) == radius
        ]
        for vector in sorted(shell, key=_candidate_key):
            key = coset_key(matrix, vector)
            representatives.setdefault(key, vector)
        radius += 1
        if radius > order + 1:
            raise RuntimeError("failed to enumerate quotient representatives")
    return sorted(representatives.values(), key=_candidate_key)


def phase_mod_one(character: Vector, element: Vector, matrix: Matrix) -> Fraction:
    """Return q^T M^-1 x mod 1 for a dual label q."""

    det = determinant(matrix)
    lifted = matrix_vector(adjugate(matrix), element)
    numerator = character[0] * lifted[0] + character[1] * lifted[1]
    value = Fraction(numerator, det)
    return value - value.numerator // value.denominator


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _group_structure(invariants: tuple[int, int]) -> str:
    factors = [f"Z/{value}" for value in invariants if value > 1]
    return " x ".join(factors) if factors else "trivial"


def quotient_data(label: str, multiplier: Vector, required: bool) -> dict:
    matrix = multiplication_matrix(multiplier)
    invariants = smith_invariants(matrix)
    exponent = invariants[1]
    elements = canonical_representatives(matrix)
    characters = canonical_representatives(transpose(matrix))
    element_index = {coset_key(matrix, value): i for i, value in enumerate(elements)}

    addition = [
        [
            element_index[
                coset_key(matrix, (left[0] + right[0], left[1] + right[1]))
            ]
            for right in elements
        ]
        for left in elements
    ]
    phase_table = [
        [phase_mod_one(character, element, matrix) for element in elements]
        for character in characters
    ]
    exponent_table = []
    for row in phase_table:
        integer_row = []
        for phase in row:
            scaled = phase * exponent
            if scaled.denominator != 1:
                raise ArithmeticError("character phase exceeds Smith exponent")
            integer_row.append(scaled.numerator % exponent)
        exponent_table.append(integer_row)

    # Exact homomorphism and character-orthogonality checks.
    for row in exponent_table:
        for first in range(len(elements)):
            for second in range(len(elements)):
                if row[addition[first][second]] != (
                    row[first] + row[second]
                ) % exponent:
                    raise ArithmeticError("character row is not a homomorphism")
    for first in range(len(characters)):
        for second in range(len(characters)):
            differences = [
                (a - b) % exponent
                for a, b in zip(exponent_table[first], exponent_table[second])
            ]
            if first == second:
                if any(differences):
                    raise ArithmeticError("diagonal character inner product failed")
                continue
            step = reduce(gcd, differences + [exponent])
            character_order = exponent // step
            histogram = Counter(differences)
            expected = len(elements) // character_order
            if histogram != Counter(
                {index * step: expected for index in range(character_order)}
            ):
                raise ArithmeticError("off-diagonal character orthogonality failed")

    return {
        "label": label,
        "required_multiplier": required,
        "multiplier": list(multiplier),
        "multiplication_matrix": [list(row) for row in matrix],
        "norm_degree": abs(determinant(matrix)),
        "smith_invariants": list(invariants),
        "additive_group": _group_structure(invariants),
        "group_exponent": exponent,
        "element_representatives": [list(value) for value in elements],
        "character_representatives_dual_quotient": [
            list(value) for value in characters
        ],
        "addition_table_indices": addition,
        "character_root_order": exponent,
        "character_exponent_table": exponent_table,
        "character_phase_table_mod_1": [
            [fraction_text(value) for value in row] for row in phase_table
        ],
        "exact_checks": {
            "character_homomorphism": True,
            "character_orthogonality": True,
            "number_of_elements": len(elements),
            "number_of_characters": len(characters),
        },
        "_matrix": matrix,
        "_elements": elements,
        "_characters": characters,
        "_phase_table": phase_table,
    }


def _associate(first: Vector, second: Vector) -> bool:
    rotations = [
        first,
        (-first[1], first[0]),
        (-first[0], -first[1]),
        (first[1], -first[0]),
    ]
    return second in rotations


def d4_operations() -> tuple[tuple[str, Matrix, bool], ...]:
    rotations = [IDENTITY]
    for _ in range(3):
        rotations.append(matrix_multiply(ROTATION_90, rotations[-1]))
    return tuple(
        [(f"rotation_{90 * index}", matrix, False) for index, matrix in enumerate(rotations)]
        + [
            (f"reflection_after_rotation_{90 * index}", matrix_multiply(matrix, CONJUGATION), True)
            for index, matrix in enumerate(rotations)
        ]
    )


def _public(group: dict) -> dict:
    return {key: value for key, value in group.items() if not key.startswith("_")}


def _target_for_reflection(source: dict, groups: dict[str, dict]) -> str:
    a, b = source["multiplier"]
    conjugate = (a, -b)
    exact = [label for label, group in groups.items() if tuple(group["multiplier"]) == conjugate]
    if exact:
        return exact[0]
    associates = [
        label
        for label, group in groups.items()
        if _associate(tuple(group["multiplier"]), conjugate)
    ]
    if not associates:
        raise ArithmeticError("missing conjugate quotient")
    if source["label"] in associates:
        return source["label"]
    return sorted(associates)[0]


def action_payload(source: dict, target: dict, transform: Matrix) -> dict:
    target_element_index = {
        coset_key(target["_matrix"], value): index
        for index, value in enumerate(target["_elements"])
    }
    images = [
        target_element_index[
            coset_key(target["_matrix"], matrix_vector(transform, element))
        ]
        for element in source["_elements"]
    ]
    if len(set(images)) != len(images):
        raise ArithmeticError("D4 element action is not bijective")

    character_images = []
    for source_row in source["_phase_table"]:
        matches = []
        for target_index, target_character in enumerate(target["_characters"]):
            pulled_row = [
                phase_mod_one(
                    target_character,
                    matrix_vector(transform, element),
                    target["_matrix"],
                )
                for element in source["_elements"]
            ]
            if pulled_row == source_row:
                matches.append(target_index)
        if len(matches) != 1:
            raise ArithmeticError("D4 character pushforward is not unique")
        character_images.append(matches[0])
    return {
        "element_image_indices": images,
        "character_pushforward_indices": character_images,
    }


def composition_oracle(groups: dict[str, dict]) -> dict:
    degree2 = groups["1+i"]
    degree4 = groups["2i"]
    bit_rows: list[tuple[int, int, int]] = []
    for index, (x, y) in enumerate(degree4["_elements"]):
        coarse = (x + y) % 2
        # z-coarse is divisible by 1+i.  Division gives
        # ((x-coarse)+y + i*(y-(x-coarse)))/2.
        real = (x - coarse + y) // 2
        imag = (y - (x - coarse)) // 2
        detail = (real + imag) % 2
        reconstructed = (coarse + detail, detail)
        if coset_key(degree4["_matrix"], reconstructed) != coset_key(
            degree4["_matrix"], (x, y)
        ):
            raise ArithmeticError("two-stage bit reconstruction failed")
        bit_rows.append((index, coarse, detail))
    if len({(coarse, detail) for _, coarse, detail in bit_rows}) != 4:
        raise ArithmeticError("two-stage bits are not a bijection")

    bit_order = [(0, 0), (0, 1), (1, 0), (1, 1)]
    element_by_bits = {
        (coarse, detail): index for index, coarse, detail in bit_rows
    }
    character_composition = []
    for coarse_frequency, detail_frequency in bit_order:
        desired = [
            (coarse_frequency * coarse + detail_frequency * detail) % 2
            for coarse, detail in bit_order
        ]
        matches = []
        for character_index, row in enumerate(degree4["character_exponent_table"]):
            reordered = [row[element_by_bits[bits]] for bits in bit_order]
            if reordered == desired:
                matches.append(character_index)
        if len(matches) != 1:
            raise ArithmeticError("two-stage character factorization failed")
        index = matches[0]
        character_composition.append(
            {
                "stage_frequencies": [coarse_frequency, detail_frequency],
                "meaning": (
                    "trivial"
                    if (coarse_frequency, detail_frequency) == (0, 0)
                    else "coarse_pullback"
                    if (coarse_frequency, detail_frequency) == (1, 0)
                    else "new_detail"
                    if (coarse_frequency, detail_frequency) == (0, 1)
                    else "coarse_times_detail"
                ),
                "degree4_character_index": index,
                "dual_representative": degree4[
                    "character_representatives_dual_quotient"
                ][index],
                "exponents_in_stage_bit_order": desired,
            }
        )

    reduction_images = []
    degree2_index = {
        coset_key(degree2["_matrix"], value): index
        for index, value in enumerate(degree2["_elements"])
    }
    for element in degree4["_elements"]:
        reduction_images.append(degree2_index[coset_key(degree2["_matrix"], element)])
    kernel = [index for index, image in enumerate(reduction_images) if image == 0]

    return {
        "identity": "(1+i)^2=2i",
        "exact_sequence": "0 -> Z/2(detail) -> K_2i -> K_(1+i) -> 0",
        "degree2_group": degree2["additive_group"],
        "degree4_group": degree4["additive_group"],
        "reduction_K2i_to_K1plusi_indices": reduction_images,
        "kernel_element_indices": kernel,
        "basis_specific_split": "z = coarse + (1+i)*detail modulo 2i",
        "stage_bit_order": [list(bits) for bits in bit_order],
        "degree4_element_indices_in_stage_bit_order": [
            element_by_bits[bits] for bits in bit_order
        ],
        "character_composition": character_composition,
        "hadamard_exponent_table_mod_2": [
            [
                (frequency[0] * bits[0] + frequency[1] * bits[1]) % 2
                for bits in bit_order
            ]
            for frequency in bit_order
        ],
        "conclusion": "the one-step degree-4 character table is exactly the tensor product of two degree-2 character tables in the declared split basis",
    }


def render() -> dict:
    groups = {
        label: quotient_data(label, multiplier, required)
        for label, multiplier, required in MULTIPLIERS
    }
    actions = {}
    for label, source in groups.items():
        entries = []
        for operation, transform, reflected in d4_operations():
            target_label = (
                _target_for_reflection(source, groups) if reflected else label
            )
            entry = {
                "operation": operation,
                "matrix": [list(row) for row in transform],
                "target_multiplier": target_label,
            }
            entry.update(action_payload(source, groups[target_label], transform))
            entries.append(entry)
        actions[label] = entries

    public_groups = {label: _public(group) for label, group in groups.items()}
    return {
        "schema": "matching-one.gaussian-cover-character-oracle.v1",
        "issue": 226,
        "claim_level": "exact_integer_character_algebra",
        "required_multipliers": [
            label for label, _, required in MULTIPLIERS if required
        ],
        "groups": public_groups,
        "d4_and_conjugation_actions": actions,
        "degree2_to_degree4_composition": composition_oracle(groups),
        "norm4_conclusion": {
            "derived_smith_invariants_for_2i": [2, 2],
            "derived_group": "Z/2 x Z/2",
            "group_exponent": 2,
            "cyclic_Z4": False,
            "gaussian_scalar_scope": "all Gaussian integers of norm 4 are associates of 2 or 2i and have Smith invariants (2,2); no cyclic norm-4 scalar Gaussian cover exists",
        },
        "interpretation_boundary": [
            "These are exact deck-group characters, not measured RG eigenoperators.",
            "D4 reflections map non-self-conjugate multipliers to their conjugate quotient.",
            "The 3-i entry is included only to close conjugation of the required 3+i oracle.",
            "The two-stage split is exact in the declared Gaussian basis; abstract extension splittings should not be called canonical without that basis.",
            "No Monte Carlo response or continuum survival claim is made.",
        ],
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(list(argv) if argv is not None else None)
    text = json.dumps(render(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
