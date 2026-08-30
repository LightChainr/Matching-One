#!/usr/bin/env python3
"""Independent Burnside certificate for bounded terminal-gadget graph orbits."""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.terminal_partition_canonical import full_symmetric_group
except ModuleNotFoundError:
    from terminal_partition_canonical import full_symmetric_group


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "gadget_burnside_orbit_certificate.json"
SCHEMA = "matching-one/gadget-burnside-orbit-certificate/v1"


def cycle_type(permutation: Sequence[int]) -> tuple[int, ...]:
    seen = set()
    lengths = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            current = permutation[current]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def induced_edge_cycle_lengths(terminal_permutation: Sequence[int]) -> tuple[int, ...]:
    terminal_count = len(terminal_permutation)
    vertex_count = terminal_count + 1
    internal = terminal_count
    mapping = tuple(terminal_permutation) + (internal,)
    slots = tuple(combinations(range(vertex_count), 2))
    slot_index = {edge: index for index, edge in enumerate(slots)}
    induced = []
    for u, v in slots:
        moved = tuple(sorted((mapping[u], mapping[v])))
        induced.append(slot_index[moved])
    return cycle_type(induced)


def multiply_by_cycle(polynomial: Sequence[int], cycle_length: int) -> list[int]:
    result = [0] * (len(polynomial) + cycle_length)
    for degree, coefficient in enumerate(polynomial):
        result[degree] += coefficient
        result[degree + cycle_length] += coefficient
    return result


def fixed_graph_edge_polynomial(edge_cycles: Sequence[int]) -> list[int]:
    polynomial = [1]
    for length in edge_cycles:
        polynomial = multiply_by_cycle(polynomial, length)
    return polynomial


def build_row(terminal_count: int) -> dict[str, Any]:
    if terminal_count not in (3, 4):
        raise ValueError("Burnside certificate supports three or four terminals")
    group = full_symmetric_group(terminal_count)
    classes: dict[tuple[int, ...], dict[str, Any]] = {}
    edge_count = (terminal_count + 1) * terminal_count // 2
    orbit_edge_numerators = [0] * (edge_count + 1)
    fixed_graph_sum = 0
    for permutation in group:
        terminal_type = cycle_type(permutation)
        edge_cycles = induced_edge_cycle_lengths(permutation)
        polynomial = fixed_graph_edge_polynomial(edge_cycles)
        fixed_graphs = sum(polynomial)
        fixed_graph_sum += fixed_graphs
        for degree, coefficient in enumerate(polynomial):
            orbit_edge_numerators[degree] += coefficient
        key = terminal_type
        if key not in classes:
            classes[key] = {
                "terminal_cycle_type": list(terminal_type),
                "class_size": 0,
                "induced_edge_cycle_lengths": list(edge_cycles),
                "fixed_graphs": fixed_graphs,
                "fixed_graph_edge_polynomial": polynomial,
            }
        row = classes[key]
        if row["induced_edge_cycle_lengths"] != list(edge_cycles):
            raise ArithmeticError("edge action varies within a conjugacy class")
        row["class_size"] += 1
    group_order = len(group)
    if fixed_graph_sum % group_order:
        raise ArithmeticError("Burnside fixed-point sum is not divisible by group order")
    if any(value % group_order for value in orbit_edge_numerators):
        raise ArithmeticError("edge-refined Burnside coefficients are not integral")
    orbit_histogram = [value // group_order for value in orbit_edge_numerators]
    return {
        "terminal_count": terminal_count,
        "internal_count": 1,
        "group_order": group_order,
        "conjugacy_classes": [classes[key] for key in sorted(classes, reverse=True)],
        "fixed_graph_sum": fixed_graph_sum,
        "canonical_orbits_by_burnside": fixed_graph_sum // group_order,
        "orbit_edge_count_histogram_by_burnside": orbit_histogram,
    }


def build_artifact() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_independent_burnside_certificate",
        "derivation": (
            "a graph fixed by a terminal permutation is constant on each induced edge-slot "
            "cycle; its edge enumerator is the product of (1+x^cycle_length)"
        ),
        "rows": [build_row(3), build_row(4)],
        "cross_check": {
            "canonical_census_orbits": [20, 90],
            "canonical_census_edge_histograms": [
                [1, 2, 4, 6, 4, 2, 1],
                [1, 2, 5, 11, 17, 18, 17, 11, 5, 2, 1],
            ],
            "exact_agreement": True,
        },
        "claim_boundary": {
            "included": "total and edge-refined graph-orbit counts by Burnside for t in {3,4}, one internal vertex",
            "excluded": "connectivity filters, probability, planarity, tiling, self-duality, ranking, thresholds, or bounds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("Burnside artifact does not exactly reproduce")
    for row, expected_orbits, expected_histogram in zip(
        artifact["rows"],
        artifact["cross_check"]["canonical_census_orbits"],
        artifact["cross_check"]["canonical_census_edge_histograms"],
    ):
        if row["canonical_orbits_by_burnside"] != expected_orbits:
            raise ValueError("Burnside orbit count disagrees with canonical census")
        if row["orbit_edge_count_histogram_by_burnside"] != expected_histogram:
            raise ValueError("edge-refined Burnside count disagrees with canonical census")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_independent_burnside_certificate",
        "orbit_counts": [row["canonical_orbits_by_burnside"] for row in artifact["rows"]],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        artifact = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_artifact(artifact), indent=2, sort_keys=True))
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
