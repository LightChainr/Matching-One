#!/usr/bin/env python3
"""Exact Green relations and ideal lattice of the typed serial monoid."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.terminal_partition_canonical import enumerate_rgs
    from scripts.terminal_partition_serial_category import serial_compose
except ModuleNotFoundError:
    from terminal_partition_canonical import enumerate_rgs
    from terminal_partition_serial_category import serial_compose

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "matching-one/terminal-partition-serial-green/v1"


def multiplication_table() -> tuple[tuple[int, ...], ...]:
    states = enumerate_rgs(4)
    index = {state: i for i, state in enumerate(states)}
    return tuple(tuple(index[serial_compose(a, b)] for b in states) for a in states)


def equivalence_classes(keys: Sequence[Any]) -> list[list[int]]:
    buckets: dict[Any, list[int]] = {}
    for i, key in enumerate(keys):
        buckets.setdefault(key, []).append(i)
    return sorted(buckets.values())


def d_classes(left_keys: Sequence[Any], right_keys: Sequence[Any]) -> list[list[int]]:
    adjacency = [set() for _ in left_keys]
    for i in range(len(left_keys)):
        for j in range(len(left_keys)):
            if left_keys[i] == left_keys[j] or right_keys[i] == right_keys[j]:
                adjacency[i].add(j)
    seen: set[int] = set()
    result = []
    for start in range(len(left_keys)):
        if start in seen:
            continue
        component = {start}
        stack = [start]
        seen.add(start)
        while stack:
            for target in adjacency[stack.pop()] - seen:
                seen.add(target)
                component.add(target)
                stack.append(target)
        result.append(sorted(component))
    return sorted(result)


def build_artifact() -> dict[str, Any]:
    states = enumerate_rgs(4)
    table = multiplication_table()
    n = len(states)
    left = [frozenset(table[s][a] for s in range(n)) for a in range(n)]
    right = [frozenset(table[a][s] for s in range(n)) for a in range(n)]
    two_sided = [
        frozenset(table[table[s][a]][t] for s in range(n) for t in range(n))
        for a in range(n)
    ]
    l_classes = equivalence_classes(left)
    r_classes = equivalence_classes(right)
    h_classes = equivalence_classes([(left[i], right[i]) for i in range(n)])
    j_classes = equivalence_classes(two_sided)
    d = d_classes(left, right)
    ideals = []
    for mask in range(1 << n):
        subset = {i for i in range(n) if mask & (1 << i)}
        if all(table[table[a][x]][b] in subset for x in subset for a in range(n) for b in range(n)):
            ideals.append(sorted(subset))
    regular = [
        a for a in range(n)
        if any(table[table[a][b]][a] == a for b in range(n))
    ]
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_green_relations_and_ideal_lattice",
        "partition_catalog": [list(state) for state in states],
        "green_classes": {"L": l_classes, "R": r_classes, "H": h_classes, "D": d, "J": j_classes},
        "class_size_profiles": {
            name: sorted(len(value) for value in classes)
            for name, classes in (("L", l_classes), ("R", r_classes), ("H", h_classes), ("D", d), ("J", j_classes))
        },
        "two_sided_ideals": ideals,
        "principal_two_sided_ideal_sizes": [len(value) for value in two_sided],
        "regular_element_indices": regular,
        "exact_checks": {
            "all_15_elements_are_regular": regular == list(range(n)),
            "finite_regular_D_equals_J": d == j_classes,
            "ideal_lattice_is_a_four_element_chain": len(ideals) == 4 and all(set(ideals[i]) < set(ideals[i + 1]) for i in range(3)),
            "unique_nontrivial_H_class_has_order_two": [c for c in h_classes if len(c) > 1] == [[6, 8]],
        },
        "claim_boundary": {
            "included": "Green L/R/H/D/J relations, regularity, and complete two-sided ideal lattice of the 15-state typed serial monoid",
            "excluded": "planarity, duality, reliability, thresholds, or infinite composition limits",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("Green-relations artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {"schema": SCHEMA, "status": "valid", "states": 15, "ideals": len(expected["two_sided_ideals"])}


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
