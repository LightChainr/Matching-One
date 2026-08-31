#!/usr/bin/env python3
"""Exhaustive automorphism and anti-automorphism census of the serial monoid."""

from __future__ import annotations

import argparse
from collections import deque
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.terminal_partition_canonical import enumerate_rgs
    from scripts.terminal_partition_serial_category import serial_compose
    from scripts.terminal_partition_serial_reversal import reverse_ports
except ModuleNotFoundError:
    from terminal_partition_canonical import enumerate_rgs
    from terminal_partition_serial_category import serial_compose
    from terminal_partition_serial_reversal import reverse_ports

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "matching-one/terminal-partition-serial-automorphisms/v1"
IDENTITY = 6


def multiplication_table() -> tuple[tuple[int, ...], ...]:
    states = enumerate_rgs(4)
    index = {state: i for i, state in enumerate(states)}
    return tuple(tuple(index[serial_compose(a, b)] for b in states) for a in states)


def closure(seed: Sequence[int], table: Sequence[Sequence[int]]) -> frozenset[int]:
    result = {IDENTITY, *seed}
    changed = True
    while changed:
        changed = False
        for a in tuple(result):
            for b in tuple(result):
                if table[a][b] not in result:
                    result.add(table[a][b])
                    changed = True
    return frozenset(result)


def word_representatives(seed: Sequence[int], table: Sequence[Sequence[int]]) -> dict[int, tuple[int, ...]]:
    words = {IDENTITY: ()}
    queue = deque([IDENTITY])
    while queue:
        source = queue.popleft()
        for position, generator in enumerate(seed):
            target = table[source][generator]
            if target not in words:
                words[target] = words[source] + (position,)
                queue.append(target)
    if len(words) != len(table):
        raise ValueError("seed does not generate every state")
    return words


def evaluate_word(word: Sequence[int], images: Sequence[int], table: Sequence[Sequence[int]]) -> int:
    value = IDENTITY
    for position in word:
        value = table[value][images[position]]
    return value


def compose_maps(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    if len(left) != len(right):
        raise ValueError("maps have different domains")
    return tuple(left[right[i]] for i in range(len(left)))


def build_artifact() -> dict[str, Any]:
    states = enumerate_rgs(4)
    table = multiplication_table()
    generating_sets = [seed for seed in combinations(range(15), 3) if len(closure(seed, table)) == 15]
    base = generating_sets[0]
    words = word_representatives(base, table)
    candidates = []
    automorphisms = set()
    for target_set in generating_sets:
        for images in permutations(target_set):
            mapping = tuple(evaluate_word(words[i], images, table) for i in range(15))
            candidates.append(mapping)
            if len(set(mapping)) == 15 and all(
                mapping[table[a][b]] == table[mapping[a]][mapping[b]]
                for a in range(15) for b in range(15)
            ):
                automorphisms.add(mapping)
    automorphisms = sorted(automorphisms)
    index = {state: i for i, state in enumerate(states)}
    reversal = tuple(index[reverse_ports(state)] for state in states)
    anti = sorted(tuple(reversal[mapping[i]] for i in range(15)) for mapping in automorphisms)
    identity_map = tuple(range(15))
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "complete_serial_monoid_automorphism_census",
        "base_generators": list(base),
        "minimum_generating_set_count": len(generating_sets),
        "ordered_generator_image_candidates": len(candidates),
        "automorphisms": [list(value) for value in automorphisms],
        "automorphism_count": len(automorphisms),
        "anti_automorphisms": [list(value) for value in anti],
        "anti_automorphism_count": len(anti),
        "exact_checks": {
            "candidate_search_is_complete_from_generators": len(candidates) == 48,
            "automorphism_group_has_order_two": len(automorphisms) == 2,
            "nontrivial_automorphism_is_an_involution": all(compose_maps(value, value) == identity_map for value in automorphisms),
            "automorphisms_are_closed_under_composition": all(compose_maps(a, b) in automorphisms for a in automorphisms for b in automorphisms),
            "anti_automorphisms_form_reversal_coset": len(anti) == len(automorphisms) == 2,
            "all_anti_maps_reverse_all_products": all(value[table[a][b]] == table[value[b]][value[a]] for value in anti for a in range(15) for b in range(15)),
        },
        "claim_boundary": {
            "included": "complete automorphism group and complete anti-automorphism coset of the finite typed serial monoid",
            "excluded": "graph-gadget automorphisms, planar duality, periodic gluing, reliability, or thresholds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("automorphism artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {"schema": SCHEMA, "status": "valid", "automorphisms": 2, "anti_automorphisms": 2}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_artifact(json.loads(args.validate.read_text())), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
