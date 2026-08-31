#!/usr/bin/env python3
"""Exact submonoid census and symmetry orbits of typed serial composition."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.terminal_partition_canonical import apply_permutation, enumerate_rgs
    from scripts.terminal_partition_serial_category import serial_compose
    from scripts.terminal_partition_serial_reversal import reverse_ports
except ModuleNotFoundError:
    from terminal_partition_canonical import apply_permutation, enumerate_rgs
    from terminal_partition_serial_category import serial_compose
    from terminal_partition_serial_reversal import reverse_ports

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "matching-one/terminal-partition-serial-submonoids/v1"
IDENTITY = 6
LANE_SWAP = (1, 0, 3, 2)


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


def image(submonoid: frozenset[int], mapping: Sequence[int]) -> frozenset[int]:
    if len(mapping) != 15 or set(mapping) != set(range(15)):
        raise ValueError("state map must be a permutation of all 15 states")
    return frozenset(mapping[value] for value in submonoid)


def build_artifact() -> dict[str, Any]:
    states = enumerate_rgs(4)
    index = {state: i for i, state in enumerate(states)}
    table = multiplication_table()
    submonoids = {
        closure([i for i in range(15) if mask & (1 << i)], table)
        for mask in range(1 << 15)
    }
    ordered = sorted((sorted(value) for value in submonoids), key=lambda value: (len(value), value))
    proper = [value for value in submonoids if len(value) < 15]
    maximal = sorted(
        (sorted(value) for value in proper if not any(value < other and len(other) < 15 for other in submonoids)),
        key=lambda value: (len(value), value),
    )
    reversal = tuple(index[reverse_ports(state)] for state in states)
    lane_swap = tuple(index[apply_permutation(state, LANE_SWAP)] for state in states)
    reverse_lane = tuple(reversal[lane_swap[i]] for i in range(15))
    transformations = (tuple(range(15)), reversal, lane_swap, reverse_lane)
    seen: set[frozenset[int]] = set()
    orbits = []
    for value in sorted(submonoids, key=lambda item: (len(item), sorted(item))):
        if value in seen:
            continue
        orbit = {image(value, mapping) for mapping in transformations}
        if not orbit <= submonoids:
            raise ValueError("symmetry image escaped the submonoid census")
        seen.update(orbit)
        orbits.append(sorted((sorted(item) for item in orbit), key=lambda item: (len(item), item)))
    orbit_sizes = Counter(map(len, orbits))
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "complete_typed_serial_submonoid_census",
        "submonoid_count": len(submonoids),
        "submonoid_size_histogram": {str(k): v for k, v in sorted(Counter(map(len, submonoids)).items())},
        "submonoids": ordered,
        "maximal_proper_submonoids": maximal,
        "symmetry": {
            "lane_swap_index_map": list(lane_swap),
            "port_reversal_index_map": list(reversal),
            "reversal_stable_count": sum(image(value, reversal) == value for value in submonoids),
            "lane_swap_stable_count": sum(image(value, lane_swap) == value for value in submonoids),
            "orbit_count": len(orbits),
            "orbit_size_histogram": {str(k): v for k, v in sorted(orbit_sizes.items())},
            "orbits": orbits,
        },
        "exact_checks": {
            "all_32768_seeds_reduce_to_228_submonoids": len(submonoids) == 228,
            "exactly_five_maximal_proper_submonoids": len(maximal) == 5,
            "symmetry_orbits_partition_the_census": seen == submonoids,
            "symmetry_orbit_profile_is_22_43_30": orbit_sizes == {1: 22, 2: 43, 4: 30},
            "reversal_and_lane_swap_stability_counts_are_exact": sum(image(value, reversal) == value for value in submonoids) == 32 and sum(image(value, lane_swap) == value for value in submonoids) == 84,
        },
        "claim_boundary": {
            "included": "all finite submonoids, maximal proper submonoids, and lane-swap/port-reversal symmetry orbits",
            "excluded": "subcategories with external objects, planar realizability, reliability, thresholds, or asymptotic growth",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("submonoid artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {"schema": SCHEMA, "status": "valid", "submonoids": 228, "symmetry_orbits": expected["symmetry"]["orbit_count"]}


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
