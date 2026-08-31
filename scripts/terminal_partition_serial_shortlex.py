#!/usr/bin/env python3
"""Canonical shortest serial words for every minimum generating set."""

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

SCHEMA = "matching-one/terminal-partition-serial-shortlex/v1"
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


def shortlex_words(seed: Sequence[int], table: Sequence[Sequence[int]]) -> dict[int, tuple[int, ...]]:
    if list(seed) != sorted(seed) or len(set(seed)) != len(seed):
        raise ValueError("generators must be distinct and increasing")
    words = {IDENTITY: ()}
    queue = deque([IDENTITY])
    while queue:
        source = queue.popleft()
        for generator in seed:
            target = table[source][generator]
            if target not in words:
                words[target] = words[source] + (generator,)
                queue.append(target)
    if len(words) != len(table):
        raise ValueError("seed does not generate every state")
    return words


def evaluate_word(word: Sequence[int], table: Sequence[Sequence[int]]) -> tuple[int, list[int]]:
    value = IDENTITY
    trace = [value]
    for generator in word:
        if type(generator) is not int or not 0 <= generator < len(table):
            raise ValueError("word generator outside Cayley table")
        value = table[value][generator]
        trace.append(value)
    return value, trace


def build_artifact() -> dict[str, Any]:
    states = enumerate_rgs(4)
    table = multiplication_table()
    full = frozenset(range(15))
    generating_sets = [seed for seed in combinations(range(15), 3) if closure(seed, table) == full]
    profiles = []
    for seed in generating_sets:
        words = shortlex_words(seed, table)
        records = []
        for target in range(15):
            value, trace = evaluate_word(words[target], table)
            records.append({"target": target, "word": list(words[target]), "length": len(words[target]), "evaluation_trace": trace, "evaluates_to": value})
        profiles.append({"generators": list(seed), "diameter": max(record["length"] for record in records), "normal_forms": records})
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "canonical_shortlex_words_for_all_minimum_generators",
        "partition_catalog": [list(state) for state in states],
        "profiles": profiles,
        "profile_count": len(profiles),
        "diameter_histogram": {str(k): v for k, v in sorted(Counter(profile["diameter"] for profile in profiles).items())},
        "exact_checks": {
            "all_eight_minimum_generating_sets_covered": len(profiles) == 8,
            "every_normal_form_evaluates_to_target": all(record["evaluates_to"] == record["target"] for profile in profiles for record in profile["normal_forms"]),
            "all_traces_have_length_plus_one": all(len(record["evaluation_trace"]) == record["length"] + 1 for profile in profiles for record in profile["normal_forms"]),
            "diameters_split_four_and_five": Counter(profile["diameter"] for profile in profiles) == {4: 4, 5: 4},
            "bfs_words_reproduce_on_second_pass": all(shortlex_words(profile["generators"], table)[record["target"]] == tuple(record["word"]) for profile in profiles for record in profile["normal_forms"]),
        },
        "claim_boundary": {
            "included": "canonical right-Cayley shortlex words, evaluation traces, and finite diameters for all minimum generating sets",
            "excluded": "preferred physical generators, planar composition words, periodic gluing, reliability, or thresholds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("shortlex artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {"schema": SCHEMA, "status": "valid", "profiles": expected["profile_count"]}


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
