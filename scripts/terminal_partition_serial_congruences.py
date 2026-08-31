#!/usr/bin/env python3
"""Exact congruence lattice of the 15-state typed serial monoid."""

from __future__ import annotations

import argparse
from collections import Counter, deque
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.terminal_partition_canonical import enumerate_rgs
    from scripts.terminal_partition_serial_category import serial_compose
except ModuleNotFoundError:
    from terminal_partition_canonical import enumerate_rgs
    from terminal_partition_serial_category import serial_compose

SCHEMA = "matching-one/terminal-partition-serial-congruences/v1"


def multiplication_table() -> tuple[tuple[int, ...], ...]:
    states = enumerate_rgs(4)
    index = {state: i for i, state in enumerate(states)}
    return tuple(tuple(index[serial_compose(a, b)] for b in states) for a in states)


def _canonical(parent: Sequence[int]) -> tuple[int, ...]:
    labels: dict[int, int] = {}
    result = []
    for value in parent:
        labels.setdefault(value, len(labels))
        result.append(labels[value])
    return tuple(result)


def generated_congruence(
    pairs: Sequence[tuple[int, int]], table: Sequence[Sequence[int]]
) -> tuple[int, ...]:
    n = len(table)
    if any(len(row) != n for row in table):
        raise ValueError("Cayley table must be square")
    parent = list(range(n))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: int, right: int) -> bool:
        if type(left) is not int or type(right) is not int or not (0 <= left < n and 0 <= right < n):
            raise ValueError("congruence pair outside Cayley table")
        left, right = find(left), find(right)
        if left == right:
            return False
        if left > right:
            left, right = right, left
        parent[right] = left
        return True

    for pair in pairs:
        union(*pair)
    changed = True
    while changed:
        changed = False
        related = [(a, b) for a in range(n) for b in range(a + 1, n) if find(a) == find(b)]
        for a, b in related:
            for x in range(n):
                changed |= union(table[x][a], table[x][b])
                changed |= union(table[a][x], table[b][x])
    roots = [find(i) for i in range(n)]
    return _canonical(roots)


def join(left: Sequence[int], right: Sequence[int], table: Sequence[Sequence[int]]) -> tuple[int, ...]:
    if len(left) != len(right) or len(left) != len(table):
        raise ValueError("congruence width mismatch")
    pairs = [
        (a, b)
        for a in range(len(left))
        for b in range(a + 1, len(left))
        if left[a] == left[b] or right[a] == right[b]
    ]
    return generated_congruence(pairs, table)


def refines(left: Sequence[int], right: Sequence[int]) -> bool:
    return all(left[a] != left[b] or right[a] == right[b] for a in range(len(left)) for b in range(a + 1, len(left)))


def congruence_classes(labels: Sequence[int]) -> list[list[int]]:
    return [[i for i, value in enumerate(labels) if value == label] for label in sorted(set(labels))]


def all_congruences(table: Sequence[Sequence[int]]) -> list[tuple[int, ...]]:
    n = len(table)
    identity = tuple(range(n))
    principals = sorted({generated_congruence([(a, b)], table) for a in range(n) for b in range(a + 1, n)})
    found = {identity}
    queue = deque([identity])
    while queue:
        current = queue.popleft()
        for principal in principals:
            value = join(current, principal, table)
            if value not in found:
                found.add(value)
                queue.append(value)
    return sorted(found, key=lambda value: (-len(set(value)), value))


def congruence_failure(labels: Sequence[int], table: Sequence[Sequence[int]]) -> Optional[dict[str, Any]]:
    if len(labels) != len(table) or any(type(value) is not int or value < 0 for value in labels):
        raise ValueError("invalid equivalence labels")
    labels = _canonical(labels)
    for a in range(len(table)):
        for b in range(a + 1, len(table)):
            if labels[a] != labels[b]:
                continue
            for x in range(len(table)):
                if labels[table[x][a]] != labels[table[x][b]]:
                    return {"related_pair": [a, b], "multiplier": x, "side": "left", "products": [table[x][a], table[x][b]]}
                if labels[table[a][x]] != labels[table[b][x]]:
                    return {"related_pair": [a, b], "multiplier": x, "side": "right", "products": [table[a][x], table[b][x]]}
    return None


def build_artifact() -> dict[str, Any]:
    states = enumerate_rgs(4)
    table = multiplication_table()
    congruences = all_congruences(table)
    covers = []
    for i, lower in enumerate(congruences):
        for j, upper in enumerate(congruences):
            if i == j or not refines(lower, upper):
                continue
            if not any(k not in (i, j) and refines(lower, middle) and refines(middle, upper) for k, middle in enumerate(congruences)):
                covers.append([i, j])
    raw_failure = None
    for a in range(15):
        for b in range(a + 1, 15):
            labels = list(range(15))
            labels[b] = labels[a]
            failure = congruence_failure(labels, table)
            if failure is not None:
                raw_failure = {"unsupported_merge": [a, b], **failure}
                break
        if raw_failure is not None:
            break
    sizes = [len(set(value)) for value in congruences]
    nontrivial = [size for size in sizes if size < 15]
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "complete_serial_monoid_congruence_lattice",
        "partition_catalog": [list(state) for state in states],
        "congruences": [{"labels": list(value), "classes": congruence_classes(value), "quotient_size": len(set(value))} for value in congruences],
        "congruence_count": len(congruences),
        "quotient_size_histogram": {str(k): v for k, v in sorted(Counter(sizes).items())},
        "cover_relations": covers,
        "largest_proper_quotient_size": max(nontrivial),
        "raw_merge_failure_witness": raw_failure,
        "exact_checks": {
            "identity_and_universal_present": sizes[0] == 15 and sizes[-1] == 1,
            "every_reported_relation_is_compatible": all(congruence_failure(value, table) is None for value in congruences),
            "closed_under_binary_joins": all(join(a, b, table) in congruences for a in congruences for b in congruences),
            "cover_relations_are_strict": all(refines(congruences[i], congruences[j]) and i != j for i, j in covers),
            "unsupported_raw_merge_has_witness": raw_failure is not None,
        },
        "claim_boundary": {
            "included": "complete two-sided congruence lattice, quotient sizes, covers, and a fail-closed raw-merge witness",
            "excluded": "preferred physical quotient, planar realization, duality, periodic gluing, reliability, or thresholds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("congruence artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {"schema": SCHEMA, "status": "valid", "congruences": expected["congruence_count"]}


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
