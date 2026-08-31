#!/usr/bin/env python3
"""Exact monoid rank and minimal generators of typed serial composition."""

from __future__ import annotations

import argparse
from collections import Counter, deque
from itertools import combinations
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
SCHEMA = "matching-one/terminal-partition-serial-generators/v1"
IDENTITY = 6


def multiplication_table() -> tuple[tuple[int, ...], ...]:
    states = enumerate_rgs(4)
    index = {state: i for i, state in enumerate(states)}
    return tuple(tuple(index[serial_compose(a, b)] for b in states) for a in states)


def generated_closure(seed: Sequence[int], table: Sequence[Sequence[int]]) -> frozenset[int]:
    n = len(table)
    if any(type(value) is not int or not 0 <= value < n for value in seed):
        raise ValueError("generator index outside Cayley table")
    closure = {IDENTITY, *seed}
    changed = True
    while changed:
        changed = False
        for left in tuple(closure):
            for right in tuple(closure):
                target = table[left][right]
                if target not in closure:
                    closure.add(target)
                    changed = True
    return frozenset(closure)


def shortest_word_lengths(seed: Sequence[int], table: Sequence[Sequence[int]]) -> tuple[int, ...]:
    lengths = {IDENTITY: 0}
    queue = deque([IDENTITY])
    while queue:
        source = queue.popleft()
        for generator in seed:
            target = table[source][generator]
            if target not in lengths:
                lengths[target] = lengths[source] + 1
                queue.append(target)
    if len(lengths) != len(table):
        raise ValueError("seed does not generate the full monoid")
    return tuple(lengths[i] for i in range(len(table)))


def build_artifact() -> dict[str, Any]:
    states = enumerate_rgs(4)
    table = multiplication_table()
    full = frozenset(range(15))
    generating_sets = []
    rank = None
    rejected = Counter()
    for size in range(4):
        current = []
        for seed in combinations(range(15), size):
            closure_size = len(generated_closure(seed, table))
            if closure_size == 15:
                current.append(seed)
            else:
                rejected[size] += 1
        if current:
            rank = size
            generating_sets = current
            break
    if rank is None:
        raise ValueError("rank search bound exhausted")
    profiles = []
    for seed in generating_sets:
        lengths = shortest_word_lengths(seed, table)
        profiles.append({"generators": list(seed), "shortest_word_lengths": list(lengths), "diameter": max(lengths)})
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_serial_monoid_rank_and_generators",
        "partition_catalog": [list(state) for state in states],
        "monoid_rank": rank,
        "minimal_generating_sets": [list(seed) for seed in generating_sets],
        "minimal_generating_set_count": len(generating_sets),
        "rejected_subset_counts_below_rank": {str(k): rejected[k] for k in range(rank)},
        "word_metric_profiles": profiles,
        "exact_checks": {
            "no_set_of_size_at_most_two_generates": sum(rejected[k] for k in range(rank)) == sum(1 for k in range(rank) for _ in combinations(range(15), k)),
            "all_reported_sets_generate_all_states": all(generated_closure(seed, table) == full for seed in generating_sets),
            "exactly_eight_minimal_rank_sets": rank == 3 and len(generating_sets) == 8,
            "shortest_word_diameters_split_four_and_five": Counter(profile["diameter"] for profile in profiles) == {4: 4, 5: 4},
        },
        "claim_boundary": {
            "included": "exact finite-monoid rank, all minimum-cardinality generating sets, and right-Cayley shortest-word profiles",
            "excluded": "physical generator interpretation, planar realization, reliability, thresholds, or asymptotic word growth",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("generator artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {"schema": SCHEMA, "status": "valid", "rank": expected["monoid_rank"], "minimum_sets": expected["minimal_generating_set_count"]}


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
