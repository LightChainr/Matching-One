#!/usr/bin/env python3
"""Exact typed serial composition of two-input/two-output partition states."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.terminal_partition_canonical import RGS, blocks_to_rgs, enumerate_rgs, rgs_to_blocks, validate_rgs
    from scripts.terminal_partition_gluing_algebra import bilinear_compose, interface_glue, interface_table
except ModuleNotFoundError:
    from terminal_partition_canonical import RGS, blocks_to_rgs, enumerate_rgs, rgs_to_blocks, validate_rgs
    from terminal_partition_gluing_algebra import bilinear_compose, interface_glue, interface_table


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "terminal_partition_serial_category_certificate.json"
SCHEMA = "matching-one/terminal-partition-serial-category/v1"


def serial_compose(left: Sequence[int], right: Sequence[int]) -> RGS:
    """Compose `(L0,L1,R0,R1)` states by identifying left R with right L."""

    return interface_glue(validate_rgs(left, 4), validate_rgs(right, 4))


def triple_graph_output(left: Sequence[int], middle: Sequence[int], right: Sequence[int]) -> RGS:
    """Independent explicit graph oracle for a three-morphism serial word."""

    states = tuple(validate_rgs(state, 4) for state in (left, middle, right))
    edges = []
    next_node = 12
    for state_index, state in enumerate(states):
        for block in rgs_to_blocks(state):
            if len(block) < 2:
                continue
            connector = next_node
            next_node += 1
            edges.extend((4 * state_index + terminal, connector) for terminal in block)
    edges.extend(((2, 4), (3, 5), (6, 8), (7, 9)))
    adjacency = [set() for _ in range(next_node)]
    for left_node, right_node in edges:
        adjacency[left_node].add(right_node)
        adjacency[right_node].add(left_node)
    components = [-1] * next_node
    component = 0
    for start in range(next_node):
        if components[start] != -1:
            continue
        components[start] = component
        stack = [start]
        while stack:
            node = stack.pop()
            for neighbor in adjacency[node]:
                if components[neighbor] == -1:
                    components[neighbor] = component
                    stack.append(neighbor)
        component += 1
    outputs = (0, 1, 10, 11)
    labels = sorted({components[node] for node in outputs})
    return blocks_to_rgs(
        (tuple(index for index, node in enumerate(outputs) if components[node] == label) for label in labels),
        4,
    )


def compose_measures(
    left: Sequence[Fraction], right: Sequence[Fraction], table: Sequence[Sequence[int]]
) -> tuple[Fraction, ...]:
    return bilinear_compose(left, right, table)


def build_artifact() -> dict[str, Any]:
    partitions = enumerate_rgs(4)
    index = {partition: position for position, partition in enumerate(partitions)}
    table = interface_table(partitions)
    identity = (0, 1, 0, 1)
    associativity_failures = 0
    graph_failures = 0
    for left in partitions:
        for middle in partitions:
            for right in partitions:
                left_grouped = serial_compose(serial_compose(left, middle), right)
                right_grouped = serial_compose(left, serial_compose(middle, right))
                associativity_failures += left_grouped != right_grouped
                graph_failures += left_grouped != triple_graph_output(left, middle, right)

    denominator = sum(range(1, 16))
    measures = [
        tuple(Fraction((index + offset) % 15 + 1, denominator) for index in range(15))
        for offset in (0, 4, 9)
    ]
    measure_left = compose_measures(compose_measures(measures[0], measures[1], table), measures[2], table)
    measure_right = compose_measures(measures[0], compose_measures(measures[1], measures[2], table), table)
    idempotents = [index[state] for state in partitions if serial_compose(state, state) == state]
    commutativity_failures = sum(
        serial_compose(left, right) != serial_compose(right, left)
        for left in partitions
        for right in partitions
    )
    identities = [
        index[state]
        for state in partitions
        if all(serial_compose(state, other) == other == serial_compose(other, state) for other in partitions)
    ]
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_typed_two_port_serial_category",
        "terminal_order": ["L0", "L1", "R0", "R1"],
        "partition_catalog": [list(value) for value in partitions],
        "serial_cayley_table": table,
        "monoid": {
            "identity_partition": list(identity),
            "identity_index": index[identity],
            "identity_candidates": identities,
            "idempotent_indices": idempotents,
            "idempotent_count": len(idempotents),
            "ordered_commutativity_failures": commutativity_failures,
        },
        "exhaustive_triples": {
            "cases": len(partitions) ** 3,
            "associativity_failures": associativity_failures,
            "explicit_graph_oracle_failures": graph_failures,
        },
        "measure_associativity_control": {
            "all_input_masses": [str(sum(measure)) for measure in measures],
            "left_grouped_mass": str(sum(measure_left)),
            "right_grouped_mass": str(sum(measure_right)),
            "vectors_equal": measure_left == measure_right,
            "output": [str(value) for value in measure_left],
        },
        "exact_checks": {
            "all_3375_state_triples_are_associative": associativity_failures == 0,
            "all_3375_triples_match_explicit_graphs": graph_failures == 0,
            "wire_identity_is_unique": identities == [index[identity]],
            "serial_product_is_not_commutative": commutativity_failures > 0,
            "exact_probability_composition_is_associative": measure_left == measure_right,
            "probability_mass_is_preserved": sum(measure_left) == sum(measure_right) == 1,
        },
        "claim_boundary": {
            "included": "typed 2-in/2-out serial composition, monoid laws, triple graph oracle, and exact measure associativity",
            "excluded": "planar restriction, duality, composition-word search, reliability inputs, thresholds, or bounds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("serial-category artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_typed_two_port_serial_category",
        "states": len(expected["partition_catalog"]),
        "triples": expected["exhaustive_triples"]["cases"],
        "idempotents": expected["monoid"]["idempotent_count"],
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
