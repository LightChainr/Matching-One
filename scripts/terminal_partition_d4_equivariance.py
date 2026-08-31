#!/usr/bin/env python3
"""Exact D4 covariance of four-terminal partition and declared-interface gluing."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.terminal_partition_canonical import (
        RGS,
        apply_permutation,
        enumerate_rgs,
        validate_group,
        validate_permutation,
        validate_rgs,
    )
    from scripts.terminal_partition_gluing_algebra import (
        DisjointSet,
        output_partition,
        partition_join,
        union_partition_blocks,
    )
except ModuleNotFoundError:
    from terminal_partition_canonical import RGS, apply_permutation, enumerate_rgs, validate_group, validate_permutation, validate_rgs
    from terminal_partition_gluing_algebra import DisjointSet, output_partition, partition_join, union_partition_blocks


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "terminal_partition_d4_equivariance_certificate.json"
SCHEMA = "matching-one/terminal-partition-d4-equivariance/v1"


def d4_group() -> tuple[tuple[int, ...], ...]:
    rotations = [tuple((index + shift) % 4 for index in range(4)) for shift in range(4)]
    reflections = [tuple((shift - index) % 4 for index in range(4)) for shift in range(4)]
    return validate_group(rotations + reflections, 4)


def declared_interface_glue(
    left: Sequence[int],
    right: Sequence[int],
    *,
    left_outer: Sequence[int],
    left_interface: Sequence[int],
    right_interface: Sequence[int],
    right_outer: Sequence[int],
) -> RGS:
    """Glue two labelled four-terminal states using explicit ordered port declarations."""

    left = validate_rgs(left, 4)
    right = validate_rgs(right, 4)
    declarations = tuple(tuple(value) for value in (left_outer, left_interface, right_interface, right_outer))
    if any(len(value) != 2 for value in declarations):
        raise ValueError("every port declaration must contain exactly two terminals")
    if set(left_outer) | set(left_interface) != set(range(4)) or set(left_outer) & set(left_interface):
        raise ValueError("left outer/interface declarations must partition the terminals")
    if set(right_outer) | set(right_interface) != set(range(4)) or set(right_outer) & set(right_interface):
        raise ValueError("right outer/interface declarations must partition the terminals")
    if any(isinstance(item, bool) or not isinstance(item, int) for declaration in declarations for item in declaration):
        raise TypeError("terminal declarations must be integers")

    dsu = DisjointSet(8)
    union_partition_blocks(dsu, left, 0)
    union_partition_blocks(dsu, right, 4)
    for left_terminal, right_terminal in zip(left_interface, right_interface):
        dsu.union(left_terminal, 4 + right_terminal)
    return output_partition(dsu, tuple(left_outer) + tuple(4 + item for item in right_outer))


def d4_orbits(partitions: Sequence[RGS], group: Sequence[Sequence[int]]) -> list[list[int]]:
    index = {partition: position for position, partition in enumerate(partitions)}
    unseen = set(range(len(partitions)))
    result = []
    while unseen:
        start = min(unseen)
        orbit = sorted({index[apply_permutation(partitions[start], permutation)] for permutation in group})
        unseen.difference_update(orbit)
        result.append(orbit)
    return result


def build_artifact() -> dict[str, Any]:
    partitions = enumerate_rgs(4)
    index = {partition: position for position, partition in enumerate(partitions)}
    group = d4_group()
    orbits = d4_orbits(partitions, group)
    join_failures = 0
    declared_interface_failures = 0
    standard = {
        "left_outer": (0, 1),
        "left_interface": (2, 3),
        "right_interface": (0, 1),
        "right_outer": (2, 3),
    }
    for permutation in group:
        moved = {
            key: tuple(permutation[item] for item in value)
            for key, value in standard.items()
        }
        for left in partitions:
            for right in partitions:
                expected_join = apply_permutation(partition_join(left, right), permutation)
                actual_join = partition_join(
                    apply_permutation(left, permutation), apply_permutation(right, permutation)
                )
                join_failures += expected_join != actual_join

                expected_interface = declared_interface_glue(left, right, **standard)
                actual_interface = declared_interface_glue(
                    apply_permutation(left, permutation),
                    apply_permutation(right, permutation),
                    **moved,
                )
                declared_interface_failures += expected_interface != actual_interface
    orbit_rows = []
    for orbit in orbits:
        stabilizers = sum(
            apply_permutation(partitions[orbit[0]], permutation) == partitions[orbit[0]]
            for permutation in group
        )
        orbit_rows.append(
            {
                "representative": list(partitions[orbit[0]]),
                "members": [list(partitions[item]) for item in orbit],
                "orbit_size": len(orbit),
                "stabilizer_size": stabilizers,
            }
        )
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_four_terminal_d4_equivariance",
        "terminal_cycle": [0, 1, 2, 3],
        "d4_group": [list(value) for value in group],
        "partition_catalog": [list(value) for value in partitions],
        "d4_orbits": orbit_rows,
        "counts": {
            "group_order": len(group),
            "partition_orbits": len(orbits),
            "join_covariance_cases": len(group) * len(partitions) ** 2,
            "declared_interface_covariance_cases": len(group) * len(partitions) ** 2,
            "orbit_size_histogram": {
                str(size): count for size, count in sorted(Counter(len(orbit) for orbit in orbits).items())
            },
        },
        "exact_checks": {
            "d4_has_order_eight": len(group) == 8,
            "join_covariance_failures_are_zero": join_failures == 0,
            "declared_interface_covariance_failures_are_zero": declared_interface_failures == 0,
            "orbit_stabilizer_products_are_eight": all(
                row["orbit_size"] * row["stabilizer_size"] == 8 for row in orbit_rows
            ),
            "orbits_partition_all_15_states": sum(len(orbit) for orbit in orbits) == 15,
        },
        "claim_boundary": {
            "included": "D4 relabeling action and covariance of join and explicitly declared two-port gluing",
            "excluded": "noncrossing selection, planar duality, self-duality, reliability, composition search, thresholds, or bounds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("D4-equivariance artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_four_terminal_d4_equivariance",
        "group_order": 8,
        "partition_orbits": expected["counts"]["partition_orbits"],
        "cases": expected["counts"]["join_covariance_cases"] + expected["counts"]["declared_interface_covariance_cases"],
    }


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
