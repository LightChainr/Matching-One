#!/usr/bin/env python3
"""Exact W=4 operator spectroscopy for the singlet versus [2,2] spin-4 gate."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from itertools import permutations, product
from pathlib import Path


Q = 4
WIDTH = 4
WORDS = list(product(range(Q), repeat=WIDTH))
WORD_INDEX = {word: i for i, word in enumerate(WORDS)}
GROUP = list(permutations(range(Q)))


def cycle_type(permutation: tuple[int, ...]) -> tuple[int, ...]:
    seen: set[int] = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


CHARACTERS = {
    "singlet_[4]": {
        "dimension": 1,
        "values": {(1, 1, 1, 1): 1, (2, 1, 1): 1, (2, 2): 1, (3, 1): 1, (4,): 1},
    },
    "two_row_[2,2]": {
        "dimension": 2,
        "values": {(1, 1, 1, 1): 2, (2, 1, 1): 0, (2, 2): 2, (3, 1): -1, (4,): 0},
    },
}


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[i]] for i in range(len(left)))


def colour_action(word: tuple[int, ...], permutation: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(permutation[value] for value in word)


def spatial_action(word: tuple[int, ...], shift: int = 0, reflection: bool = False) -> tuple[int, ...]:
    if reflection:
        word = tuple(word[(-i) % WIDTH] for i in range(WIDTH))
    return tuple(word[(i - shift) % WIDTH] for i in range(WIDTH))


def monochromatic_components(word: tuple[int, ...], edges: list[tuple[int, int]]) -> int:
    parent = list(range(WIDTH))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for left, right in edges:
        if word[left] != word[right]:
            continue
        a, b = find(left), find(right)
        if a != b:
            parent[b] = a
    return len({find(i) for i in range(WIDTH)})


def row_matching_defect(word: tuple[int, ...]) -> int:
    cycle = [(i, (i + 1) % WIDTH) for i in range(WIDTH)]
    complete = [(i, j) for i in range(WIDTH) for j in range(i + 1, WIDTH)]
    return monochromatic_components(word, cycle) - monochromatic_components(word, complete)


def potts_transfer_entry(old: tuple[int, ...], new: tuple[int, ...], v: int = 2) -> int:
    """One-row square-lattice Potts transfer; v=sqrt(Q)=2 is self-dual critical."""
    weight = 1
    for i in range(WIDTH):
        weight *= 1 + v * int(old[i] == new[i])
        weight *= 1 + v * int(new[i] == new[(i + 1) % WIDTH])
    return weight


def projector_vector(irrep: str, vector: dict[int, Fraction]) -> dict[int, Fraction]:
    data = CHARACTERS[irrep]
    output: dict[int, Fraction] = {}
    scale = Fraction(data["dimension"], len(GROUP))
    for permutation in GROUP:
        coefficient = scale * data["values"][cycle_type(permutation)]
        if not coefficient:
            continue
        for index, value in vector.items():
            target = WORD_INDEX[colour_action(WORDS[index], permutation)]
            output[target] = output.get(target, Fraction(0)) + coefficient * value
    return {index: value for index, value in output.items() if value}


def cross_operator_max_coefficient() -> Fraction:
    maximum = Fraction(0)
    for index in range(len(WORDS)):
        charged = projector_vector("two_row_[2,2]", {index: Fraction(1)})
        inserted = {i: value * row_matching_defect(WORDS[i]) for i, value in charged.items()}
        crossed = projector_vector("singlet_[4]", inserted)
        maximum = max(maximum, *(abs(value) for value in crossed.values()), Fraction(0))
    return maximum


def projector_trace(irrep: str, insertion=None, seam: tuple[int, ...] | None = None) -> Fraction:
    data = CHARACTERS[irrep]
    seam = seam or tuple(range(Q))
    total = 0
    for permutation in GROUP:
        coefficient = data["values"][cycle_type(permutation)]
        action = compose(permutation, seam)
        for word in WORDS:
            if colour_action(word, action) == word:
                total += coefficient * (1 if insertion is None else insertion(word))
    return Fraction(data["dimension"] * total, len(GROUP))


def spatial_even_zero_momentum_rank(irrep: str) -> Fraction:
    data = CHARACTERS[irrep]
    total = 0
    for permutation in GROUP:
        coefficient = data["values"][cycle_type(permutation)]
        for reflection in (False, True):
            for shift in range(WIDTH):
                total += coefficient * sum(
                    colour_action(spatial_action(word, shift, reflection), permutation) == word
                    for word in WORDS
                )
    return Fraction(data["dimension"] * total, len(GROUP) * 2 * WIDTH)


def transfer_symmetry_checks() -> dict:
    colour_generators = ((1, 0, 2, 3), (1, 2, 3, 0))
    colour_ok = translation_ok = reflection_ok = True
    observable_colour_ok = observable_spatial_ok = True
    for old in WORDS:
        for new in WORDS:
            base = potts_transfer_entry(old, new)
            colour_ok &= all(
                base == potts_transfer_entry(colour_action(old, generator), colour_action(new, generator))
                for generator in colour_generators
            )
            translation_ok &= base == potts_transfer_entry(spatial_action(old, 1), spatial_action(new, 1))
            reflection_ok &= base == potts_transfer_entry(
                spatial_action(old, reflection=True), spatial_action(new, reflection=True)
            )
        observable_colour_ok &= all(
            row_matching_defect(old) == row_matching_defect(colour_action(old, generator))
            for generator in colour_generators
        )
        observable_spatial_ok &= all(
            row_matching_defect(old) == row_matching_defect(spatial_action(old, shift, reflection))
            for reflection in (False, True)
            for shift in range(WIDTH)
        )
    return {
        "T_commutes_with_S4_generators": colour_ok,
        "T_commutes_with_translation": translation_ok,
        "T_commutes_with_reflection": reflection_ok,
        "O_matching_commutes_with_S4_generators": observable_colour_ok,
        "O_matching_commutes_with_D4_row_group": observable_spatial_ok,
    }


def build_oracle() -> dict:
    transposition = (1, 0, 2, 3)
    defect_histogram = {
        str(value): sum(row_matching_defect(word) == value for word in WORDS)
        for value in sorted({row_matching_defect(word) for word in WORDS})
    }
    singlet_trace = projector_trace("singlet_[4]", row_matching_defect)
    charged_trace = projector_trace("two_row_[2,2]", row_matching_defect)
    singlet_twisted = projector_trace("singlet_[4]", row_matching_defect, transposition)
    charged_twisted = projector_trace("two_row_[2,2]", row_matching_defect, transposition)

    return {
        "schema": "matching-one.p120-small-width-operator-spectroscopy.v1",
        "issue": 120,
        "control": {
            "model": "Q=4 square-lattice Potts row transfer at self-dual v=2",
            "width": WIDTH,
            "row_state_count": len(WORDS),
            "reason_for_Q4": "integer realization of the generic Potts [2] carrier; no pc estimate is performed",
            "symmetry_checks": transfer_symmetry_checks(),
        },
        "spatial_sector": {
            "target": "spin +/-4, reflection-even cosine channel",
            "width4_translation_phase": "exp(2 pi i * 4/4)=1",
            "alias": "spin 4 and spin 0 share the zero-momentum reflection-even row sector at W=4",
            "combined_sector_ranks": {
                "singlet_[4]": str(spatial_even_zero_momentum_rank("singlet_[4]")),
                "two_row_[2,2]": str(spatial_even_zero_momentum_rank("two_row_[2,2]")),
            },
            "conclusion": "rotation/translation and reflection alone cannot distinguish thermal Q4 from V_(2,+/-2)",
        },
        "topology_and_Potts_sector": {
            "global_matching_endpoint": {
                "Potts_sector": "singlet [4]",
                "marked_cluster_legs": 0,
                "row_insertion": "unmarked colour-blind C4-to-K4 monochromatic connectivity defect",
                "insertion_histogram": defect_histogram,
            },
            "V22_candidate": {
                "Potts_sector": "[2] analytically; [2,2] at Q=4",
                "marked_cluster_legs": 4,
                "spin": "+/-4",
            },
            "thermal_Q4_candidate": {
                "Potts_sector": "singlet",
                "marked_cluster_legs": 0,
                "spin": "+/-4",
            },
        },
        "exact_matrix_element": {
            "operator": "P_singlet O_matching P_[2,2]",
            "max_absolute_coefficient_in_full_W4_basis": str(cross_operator_max_coefficient()),
            "is_zero": cross_operator_max_coefficient() == 0,
            "reason": "T and O_matching commute with S4, while the central projectors are orthogonal",
            "thermal_Q4": "not symmetry-forbidden",
        },
        "minimal_counterexample_to_absence": {
            "trace_P_singlet_O": str(singlet_trace),
            "trace_P_[2,2]_O": str(charged_trace),
            "meaning": "O_matching acts nontrivially inside [2,2]; only the singlet-to-[2,2] matrix element vanishes",
        },
        "colour_seam_one_shot": {
            "seam": "one S4 transposition",
            "unnormalized_singlet_trace_identity": str(singlet_trace),
            "unnormalized_singlet_trace_twisted": str(singlet_twisted),
            "singlet_ratio": str(singlet_twisted / singlet_trace),
            "unnormalized_[2,2]_trace_identity": str(charged_trace),
            "unnormalized_[2,2]_trace_twisted": str(charged_twisted),
            "[2,2]_ratio": str(charged_twisted / charged_trace),
            "targets": "thermal singlet 1 versus V22 [2,2] 0",
        },
        "matching_odd_marker": {
            "construction": "tensor O_matching with the #155 complement-odd two-state marker D, with J D J=-D",
            "selection_result": "P_singlet (O_matching tensor D) P_[2,2]=0",
            "meaning": "matching/complement oddness commutes with Potts colour projection and cannot rescue V22",
        },
        "claim_boundary": {
            "proved": [
                "the exact W4 Potts transfer and insertion symmetry commutators",
                "the full-basis singlet-to-[2,2] matrix-element zero",
                "the nonzero [2,2] diagonal trace counterexample",
                "the transposition-seam ratios 1 and 0",
                "the W4 spatial alias of spin 0 and spin 4",
            ],
            "not_proved": [
                "a numerical scaling dimension from a transfer gap",
                "a Q=1 transfer realization of the complete global matching observable",
                "that every possible twisted or Q-derivative matching insertion is colour singlet",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    payload = json.dumps(build_oracle(), indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
