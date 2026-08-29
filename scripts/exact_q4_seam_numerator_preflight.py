#!/usr/bin/env python3
"""Exact Q=4 seam/numerator preflight linking Issues 257 and 258.

The tiny row transfer is an actual critical Q=4 Potts transfer matrix.  Its
central S4 projectors make the distinction between an unnormalized sector
numerator and a seam-normalized expectation executable before an H4-specific
production kernel exists.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path
from typing import Sequence


HERE = Path(__file__).resolve().parent
SELECTION_PATH = HERE / "global_matching_spin4_selection.py"
SPEC = importlib.util.spec_from_file_location("global_matching_selection", SELECTION_PATH)
SELECTION = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(SELECTION)

Matrix = list[list[Fraction]]
Permutation = tuple[int, ...]
Q = 4
V_CRITICAL = 2
GROUP_ORDER = 24
IRREPS = ("singlet", "two_row_2")


def fraction_record(value: Fraction) -> dict[str, int | str]:
    return {"numerator": value.numerator, "denominator": value.denominator,
            "text": str(value)}


def multiply(left: Matrix, right: Matrix) -> Matrix:
    return [[sum(left[i][k] * right[k][j] for k in range(len(right)))
             for j in range(len(right[0]))] for i in range(len(left))]


def add(first: Matrix, second: Matrix) -> Matrix:
    return [[a + b for a, b in zip(row_a, row_b)]
            for row_a, row_b in zip(first, second)]


def scale(matrix: Matrix, scalar: Fraction) -> Matrix:
    return [[scalar * value for value in row] for row in matrix]


def identity(size: int) -> Matrix:
    return [[Fraction(i == j) for j in range(size)] for i in range(size)]


def power(matrix: Matrix, exponent: int) -> Matrix:
    result = identity(len(matrix))
    factor = matrix
    while exponent:
        if exponent & 1:
            result = multiply(result, factor)
        factor = multiply(factor, factor)
        exponent //= 2
    return result


def trace(matrix: Matrix) -> Fraction:
    return sum(matrix[i][i] for i in range(len(matrix)))


def cycle_type(permutation: Permutation) -> tuple[int, ...]:
    seen: set[int] = set()
    cycles = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        length = 0
        value = start
        while value not in seen:
            seen.add(value)
            length += 1
            value = permutation[value]
        cycles.append(length)
    return tuple(sorted(cycles, reverse=True))


def character(irrep: str, permutation: Permutation) -> int:
    table = {
        "singlet": {(1, 1, 1, 1): 1, (2, 1, 1): 1, (2, 2): 1,
                    (3, 1): 1, (4,): 1},
        "two_row_2": {(1, 1, 1, 1): 2, (2, 1, 1): 0, (2, 2): 2,
                      (3, 1): -1, (4,): 0},
    }
    return table[irrep][cycle_type(permutation)]


def row_states(width: int) -> list[tuple[int, ...]]:
    return list(product(range(Q), repeat=width))


def colour_action(width: int, permutation: Permutation) -> Matrix:
    states = row_states(width)
    index = {state: i for i, state in enumerate(states)}
    matrix = [[Fraction(0) for _ in states] for _ in states]
    for column, state in enumerate(states):
        image = tuple(permutation[colour] for colour in state)
        matrix[index[image]][column] = Fraction(1)
    return matrix


def transfer_and_logv_derivative(width: int) -> tuple[Matrix, Matrix]:
    """Unsymmetrized row transfer and its total log(v) derivative at v=2."""
    states = row_states(width)
    transfer: Matrix = []
    derivative: Matrix = []
    for old in states:
        transfer_row = []
        derivative_row = []
        for new in states:
            factors = []
            for x in range(width):
                factors.append(old[x] == new[x])
                factors.append(new[x] == new[(x + 1) % width])
            weight = Fraction(1)
            logarithmic_derivative = Fraction(0)
            for equal in factors:
                factor = 1 + V_CRITICAL * int(equal)
                weight *= factor
                if equal:
                    logarithmic_derivative += Fraction(V_CRITICAL, factor)
            transfer_row.append(weight)
            derivative_row.append(weight * logarithmic_derivative)
        transfer.append(transfer_row)
        derivative.append(derivative_row)
    return transfer, derivative


def inserted_power_derivative(transfer: Matrix, derivative: Matrix, height: int) -> Matrix:
    result = [[Fraction(0) for _ in transfer] for _ in transfer]
    for location in range(height):
        term = multiply(power(transfer, location),
                        multiply(derivative, power(transfer, height - 1 - location)))
        result = add(result, term)
    return result


def central_projector(width: int, irrep: str) -> Matrix:
    dimension = character(irrep, tuple(range(Q)))
    result = [[Fraction(0) for _ in row_states(width)] for _ in row_states(width)]
    for permutation in permutations(range(Q)):
        result = add(result, scale(colour_action(width, permutation),
                                   Fraction(dimension * character(irrep, permutation), GROUP_ORDER)))
    return result


def commutator_zero(first: Matrix, second: Matrix) -> bool:
    return multiply(first, second) == multiply(second, first)


def sector_numerator(projector: Matrix, seam: Matrix, insertion: Matrix) -> Fraction:
    return trace(multiply(projector, multiply(seam, insertion)))


def direct_spin_partition_sum(width: int, height: int,
                              seam_permutation: Permutation) -> Fraction:
    """Independent spin-configuration check of the time-cycle seam trace."""
    total = Fraction(0)
    for flat in product(range(Q), repeat=width * height):
        weight = Fraction(1)
        for y in range(height):
            for x in range(width):
                colour = flat[y * width + x]
                horizontal = flat[y * width + (x + 1) % width]
                next_y = (y + 1) % height
                vertical = flat[next_y * width + x]
                if y == height - 1:
                    vertical = seam_permutation[vertical]
                weight *= 1 + V_CRITICAL * int(colour == horizontal)
                weight *= 1 + V_CRITICAL * int(colour == vertical)
        total += weight
    return total


def numerator_record(projector: Matrix, insertion: Matrix, seams: dict[str, Matrix],
                     partition_sums: dict[str, Fraction]) -> dict:
    values = {name: sector_numerator(projector, seam, insertion)
              for name, seam in seams.items()}
    ratio = values["transposition"] / values["identity"]
    normalized_ratio = ((values["transposition"] / partition_sums["transposition"]) /
                        (values["identity"] / partition_sums["identity"]))
    restored = normalized_ratio * partition_sums["transposition"] / partition_sums["identity"]
    return {
        "unnormalized_numerator": {name: fraction_record(value) for name, value in values.items()},
        "unnormalized_twist_to_identity_ratio": fraction_record(ratio),
        "normalized_expectation_twist_to_identity_ratio": fraction_record(normalized_ratio),
        "partition_restored_ratio": fraction_record(restored),
    }


def render(width: int = 2, height: int = 2) -> dict:
    if (width, height) != (2, 2):
        raise ValueError("the frozen exact preflight is the 2x2 torus")
    transfer, derivative = transfer_and_logv_derivative(width)
    transfer_power = power(transfer, height)
    inserted = inserted_power_derivative(transfer, derivative, height)
    seam_permutations = {
        "identity": tuple(range(Q)),
        "transposition": (1, 0, 2, 3),
    }
    seams = {name: colour_action(width, permutation)
             for name, permutation in seam_permutations.items()}
    partition_sums = {name: trace(multiply(seam, transfer_power))
                      for name, seam in seams.items()}
    direct_partition_sums = {
        name: direct_spin_partition_sum(width, height, permutation)
        for name, permutation in seam_permutations.items()
    }
    selection_targets = SELECTION.s4_transposition_oracle()
    sectors = {}
    for irrep in IRREPS:
        projector = central_projector(width, irrep)
        sectors[irrep] = {
            "S4_irrep": {"singlet": "[4]", "two_row_2": "[2,2]"}[irrep],
            "projector_rank": fraction_record(trace(projector)),
            "projector_commutes_with_transfer": commutator_zero(projector, transfer),
            "seam_character_target": selection_targets[irrep]["twist_to_identity_trace_ratio"],
            "transfer_power_witness": numerator_record(
                projector, transfer_power, seams, partition_sums),
            "logv_inserted_witness": numerator_record(
                projector, inserted, seams, partition_sums),
        }
    return {
        "schema": "matching-one.exact-q4-seam-numerator-preflight.v1",
        "issues": [257, 258],
        "status": "exact_tiny_preflight_complete_h4_production_not_run",
        "geometry": {
            "model": "critical Q=4 square-lattice Potts row transfer",
            "Q": Q,
            "v": V_CRITICAL,
            "width": width,
            "height": height,
            "row_states": Q ** width,
            "spin_configurations": Q ** (width * height),
            "parallel_edge_note": "width=2 torus retains the two directed nearest-neighbour bonds per row",
        },
        "seams": {
            "identity": list(seam_permutations["identity"]),
            "single_colour_transposition": list(seam_permutations["transposition"]),
            "placement": "time-cycle trace Tr[U_g K]",
        },
        "partition_sums": {name: fraction_record(value) for name, value in partition_sums.items()},
        "direct_spin_enumeration_partition_sums": {
            name: fraction_record(value) for name, value in direct_partition_sums.items()
        },
        "partition_twist_to_identity_ratio": fraction_record(
            partition_sums["transposition"] / partition_sums["identity"]),
        "sector_numerator_oracles": sectors,
        "exact_checks": {
            "all_projectors_commute_with_transfer": all(
                row["projector_commutes_with_transfer"] for row in sectors.values()),
            "both_witnesses_match_frozen_character_targets": all(
                row[witness]["unnormalized_twist_to_identity_ratio"]["text"]
                == row["seam_character_target"]["text"]
                for row in sectors.values()
                for witness in ("transfer_power_witness", "logv_inserted_witness")),
            "normalization_restoration_recovers_numerator_ratio": all(
                row[witness]["partition_restored_ratio"]
                == row[witness]["unnormalized_twist_to_identity_ratio"]
                for row in sectors.values()
                for witness in ("transfer_power_witness", "logv_inserted_witness")),
            "partition_function_is_seam_dependent": partition_sums["identity"]
            != partition_sums["transposition"],
            "row_transfer_equals_direct_spin_enumeration": partition_sums == direct_partition_sums,
        },
        "scope_boundary": {
            "proved": (
                "For an S4-equivariant transfer insertion inside a projected isotypic sector, the unnormalized "
                "identity/transposition numerator ratio is chi_lambda((12))/dim(lambda): 1 for the singlet and "
                "0 for [2,2]. The seam partition factor must be restored if normalized expectations are stored."
            ),
            "not_proved": (
                "The log(v) witness is a nonzero equivariant insertion, not a pure continuum H4 operator. "
                "A production runner must implement the separately frozen orientation-resolved H4 kernel."
            ),
            "issue_258_cross_validation": (
                "The observable numerator and the partition normalization are separate derivatives/objects; "
                "a seam-dependent denominator cannot be absorbed into the representation character."
            ),
        },
        "frozen_minimal_acquisition_schema": {
            "status": "prereveal_no_production_started",
            "model": "critical Q=4 square-lattice Potts/FK transfer",
            "seams": ["identity", "single transposition (01) on one fixed torus cycle"],
            "seam_orientation": "repeat with the same cycle and H4 handedness; do not average conjugacy representatives",
            "required_raw_outputs": [
                "Z_identity", "Z_transposition",
                "N_H4_singlet_identity", "N_H4_singlet_transposition",
                "N_H4_[2,2]_identity", "N_H4_[2,2]_transposition",
            ],
            "primary_scores": {
                "singlet": "N_H4_singlet_transposition/N_H4_singlet_identity = 1",
                "two_row_2": "N_H4_[2,2]_transposition/N_H4_[2,2]_identity = 0",
            },
            "normalized_output_conversion": (
                "If E_g=N_g/Z_g is recorded, score (E_transposition/E_identity)*"
                "(Z_transposition/Z_identity), never E_transposition/E_identity alone."
            ),
            "acceptance_gates": [
                "tiny exact identity/transposition Z values reproduce this artifact",
                "sector projectors commute with the zero-source transfer",
                "identity-seam H4 numerators are nonzero in both target sectors",
                "the H4 kernel has a declared orientation/handedness and is not the log(v) witness",
            ],
            "reason_no_large_run": "the repository currently has no colour-seam H4 kernel satisfying the last gate",
        },
        "scientific_layers": {
            "exact": "finite Q4 transfer/projector/seam traces and normalization restoration",
            "mechanism_inference": "one seam separates singlet and [2,2] colour charge without an exponent fit",
            "exploratory": "a future nonzero [2,2] H4 numerator would be the controlled V22-like channel needed to compare with local pivotal H4",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--width", type=int, default=2)
    parser.add_argument("--height", type=int, default=2)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render(args.width, args.height)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
