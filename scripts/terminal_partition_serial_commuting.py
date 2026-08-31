#!/usr/bin/env python3
"""Exact center and commutative-submonoid census for typed serial composition."""

from __future__ import annotations

import argparse
from collections import Counter
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
SCHEMA = "matching-one/terminal-partition-serial-commuting/v1"
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


def all_submonoids(table: Sequence[Sequence[int]]) -> set[frozenset[int]]:
    return {
        closure([i for i in range(15) if mask & (1 << i)], table)
        for mask in range(1 << 15)
    }


def build_artifact() -> dict[str, Any]:
    table = multiplication_table()
    commuting = [[table[a][b] == table[b][a] for b in range(15)] for a in range(15)]
    center = [a for a in range(15) if all(commuting[a])]
    submonoids = all_submonoids(table)
    commutative = [value for value in submonoids if all(commuting[a][b] for a in value for b in value)]
    maximal = sorted(
        (sorted(value) for value in commutative if not any(value < other for other in commutative)),
        key=lambda value: (len(value), value),
    )
    vertices = [i for i in range(15) if i not in center]
    degrees = [sum(commuting[v][w] for w in vertices if w != v) for v in vertices]
    seen: set[int] = set()
    components = []
    for start in vertices:
        if start in seen:
            continue
        component = {start}
        stack = [start]
        seen.add(start)
        while stack:
            source = stack.pop()
            for target in vertices:
                if target not in seen and target != source and commuting[source][target]:
                    seen.add(target)
                    component.add(target)
                    stack.append(target)
        components.append(sorted(component))
    clique_number = 0
    for mask in range(1 << len(vertices)):
        subset = [vertices[i] for i in range(len(vertices)) if mask & (1 << i)]
        if len(subset) > clique_number and all(commuting[a][b] for a in subset for b in subset):
            clique_number = len(subset)
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_center_and_commutative_submonoids",
        "center_indices": center,
        "ordered_commuting_pair_count": sum(sum(row) for row in commuting),
        "unordered_distinct_commuting_pair_count": sum(commuting[a][b] for a in range(15) for b in range(a + 1, 15)),
        "commuting_graph_without_center": {
            "vertices": vertices,
            "degree_sequence": sorted(degrees),
            "components": components,
            "clique_number": clique_number,
        },
        "commutative_submonoid_count": len(commutative),
        "commutative_submonoid_size_histogram": {str(k): v for k, v in sorted(Counter(map(len, commutative)).items())},
        "maximal_commutative_submonoids": maximal,
        "exact_checks": {
            "center_is_only_wire_identity": center == [IDENTITY],
            "exactly_29_commutative_submonoids": len(commutative) == 29,
            "exactly_11_are_inclusion_maximal": len(maximal) == 11,
            "center_removed_commuting_graph_is_connected": len(components) == 1,
            "center_removed_clique_number_is_three": clique_number == 3,
        },
        "claim_boundary": {
            "included": "center, commuting graph, all commutative submonoids, and inclusion-maximal commutative sectors",
            "excluded": "probabilistic commutation, physical independence, planar duality, reliability, or thresholds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("commuting artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {"schema": SCHEMA, "status": "valid", "center": 1, "commutative_submonoids": 29}


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
