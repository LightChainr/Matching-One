#!/usr/bin/env python3
"""Exact finite terminal-partition gluing algebra for Issue 13."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Union

try:
    from scripts.terminal_partition_canonical import (
        RGS,
        blocks_to_rgs,
        enumerate_rgs,
        rgs_to_blocks,
        validate_rgs,
    )
except ModuleNotFoundError:
    from terminal_partition_canonical import RGS, blocks_to_rgs, enumerate_rgs, rgs_to_blocks, validate_rgs


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "terminal_partition_gluing_algebra_certificate.json"
SCHEMA = "matching-one/terminal-partition-gluing-algebra/v1"
ExactInput = Union[int, str, Fraction]


class DisjointSet:
    def __init__(self, size: int) -> None:
        if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
            raise ValueError("disjoint-set size must be a positive integer")
        self.parent = list(range(size))

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def exact_fraction(value: ExactInput, *, field: str) -> Fraction:
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError(f"{field} must be exact; floats and booleans are forbidden")
    try:
        return value if isinstance(value, Fraction) else Fraction(value)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"invalid exact value for {field}") from exc


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}" if value.denominator != 1 else str(value.numerator)


def union_partition_blocks(dsu: DisjointSet, partition: Sequence[int], node_offset: int = 0) -> None:
    valid = validate_rgs(partition)
    for block in rgs_to_blocks(valid):
        anchor = node_offset + block[0]
        for terminal in block[1:]:
            dsu.union(anchor, node_offset + terminal)


def output_partition(dsu: DisjointSet, output_nodes: Sequence[int]) -> RGS:
    if len(set(output_nodes)) != len(output_nodes) or not output_nodes:
        raise ValueError("output nodes must be distinct and nonempty")
    roots = [dsu.find(node) for node in output_nodes]
    labels = {}
    encoded = []
    for root in roots:
        if root not in labels:
            labels[root] = len(labels)
        encoded.append(labels[root])
    return validate_rgs(encoded, len(output_nodes))


def partition_join(left: Sequence[int], right: Sequence[int]) -> RGS:
    left_valid = validate_rgs(left)
    right_valid = validate_rgs(right, len(left_valid))
    dsu = DisjointSet(len(left_valid))
    union_partition_blocks(dsu, left_valid)
    union_partition_blocks(dsu, right_valid)
    return output_partition(dsu, tuple(range(len(left_valid))))


def interface_glue(left: Sequence[int], right: Sequence[int]) -> RGS:
    """Glue (A,B,x,y) to (x',y',C,D) and retain (A,B,C,D)."""

    left_valid = validate_rgs(left, 4)
    right_valid = validate_rgs(right, 4)
    dsu = DisjointSet(8)
    union_partition_blocks(dsu, left_valid, 0)
    union_partition_blocks(dsu, right_valid, 4)
    dsu.union(2, 4)
    dsu.union(3, 5)
    return output_partition(dsu, (0, 1, 6, 7))


def representative_graph(partition: Sequence[int], terminal_offset: int, first_internal: int) -> tuple[list[tuple[int, int]], int]:
    valid = validate_rgs(partition, 4)
    edges = []
    next_internal = first_internal
    for block in rgs_to_blocks(valid):
        if len(block) < 2:
            continue
        connector = next_internal
        next_internal += 1
        edges.extend((terminal_offset + terminal, connector) for terminal in block)
    return edges, next_internal


def graph_composition_output(left: Sequence[int], right: Sequence[int]) -> RGS:
    """Independently realize partitions as tiny star forests and glue the interface."""

    left_valid = validate_rgs(left, 4)
    right_valid = validate_rgs(right, 4)
    left_edges, next_node = representative_graph(left_valid, 0, 8)
    right_edges, next_node = representative_graph(right_valid, 4, next_node)
    edges = left_edges + right_edges + [(2, 4), (3, 5)]
    adjacency = [set() for _ in range(next_node)]
    for u, v in edges:
        adjacency[u].add(v)
        adjacency[v].add(u)
    component = [-1] * next_node
    component_id = 0
    for start in range(next_node):
        if component[start] != -1:
            continue
        component[start] = component_id
        stack = [start]
        while stack:
            vertex = stack.pop()
            for neighbor in adjacency[vertex]:
                if component[neighbor] == -1:
                    component[neighbor] = component_id
                    stack.append(neighbor)
        component_id += 1
    return blocks_to_rgs(
        (
            tuple(index for index, node in enumerate((0, 1, 6, 7)) if component[node] == label)
            for label in sorted({component[node] for node in (0, 1, 6, 7)})
        ),
        4,
    )


def interface_table(partitions: Sequence[RGS]) -> list[list[int]]:
    index = {partition: position for position, partition in enumerate(partitions)}
    if len(index) != len(partitions):
        raise ValueError("partition catalog contains duplicates")
    return [[index[interface_glue(left, right)] for right in partitions] for left in partitions]


def join_table(partitions: Sequence[RGS]) -> list[list[int]]:
    index = {partition: position for position, partition in enumerate(partitions)}
    if len(index) != len(partitions):
        raise ValueError("partition catalog contains duplicates")
    return [[index[partition_join(left, right)] for right in partitions] for left in partitions]


def bilinear_compose(
    left_weights: Sequence[ExactInput], right_weights: Sequence[ExactInput], table: Sequence[Sequence[int]]
) -> tuple[Fraction, ...]:
    if not table or any(len(row) != len(table) for row in table):
        raise ValueError("gluing table must be nonempty and square")
    size = len(table)
    if len(left_weights) != size or len(right_weights) != size:
        raise ValueError("weight vectors must match the gluing table")
    left = tuple(exact_fraction(value, field=f"left_weights[{index}]") for index, value in enumerate(left_weights))
    right = tuple(exact_fraction(value, field=f"right_weights[{index}]") for index, value in enumerate(right_weights))
    result = [Fraction(0) for _ in range(size)]
    for left_index, left_weight in enumerate(left):
        for right_index, right_weight in enumerate(right):
            output_index = table[left_index][right_index]
            if isinstance(output_index, bool) or not isinstance(output_index, int) or not 0 <= output_index < size:
                raise ValueError("gluing table output index out of range")
            result[output_index] += left_weight * right_weight
    return tuple(result)


def build_artifact() -> dict[str, Any]:
    partitions = enumerate_rgs(4)
    index = {partition: position for position, partition in enumerate(partitions)}
    joins = join_table(partitions)
    interface = interface_table(partitions)
    discrete = index[(0, 1, 2, 3)]
    connected = index[(0, 0, 0, 0)]

    associativity_failures = 0
    for left in partitions:
        for middle in partitions:
            for right in partitions:
                if partition_join(partition_join(left, middle), right) != partition_join(left, partition_join(middle, right)):
                    associativity_failures += 1

    graph_control_failures = 0
    for left in partitions:
        for right in partitions:
            if interface_glue(left, right) != graph_composition_output(left, right):
                graph_control_failures += 1

    denominator = sum(range(1, len(partitions) + 1))
    left_weights = tuple(Fraction(index_value + 1, denominator) for index_value in range(len(partitions)))
    right_weights = tuple(reversed(left_weights))
    composed = bilinear_compose(left_weights, right_weights, interface)
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_terminal_partition_gluing_algebra",
        "partition_catalog": [list(partition) for partition in partitions],
        "join_cayley_table": joins,
        "join_laws": {
            "triples_checked_for_associativity": len(partitions) ** 3,
            "associativity_failures": associativity_failures,
            "commutative": all(joins[i][j] == joins[j][i] for i in range(len(partitions)) for j in range(len(partitions))),
            "idempotent": all(joins[i][i] == i for i in range(len(partitions))),
            "discrete_identity_index": discrete,
            "discrete_is_identity": all(joins[discrete][i] == i == joins[i][discrete] for i in range(len(partitions))),
            "connected_absorber_index": connected,
            "connected_is_absorbing": all(joins[connected][i] == connected == joins[i][connected] for i in range(len(partitions))),
        },
        "two_port_interface": {
            "left_terminals": ["A", "B", "x", "y"],
            "right_terminals": ["x_prime", "y_prime", "C", "D"],
            "identifications": [["x", "x_prime"], ["y", "y_prime"]],
            "surviving_output": ["A", "B", "C", "D"],
            "deterministic_output_index_table": interface,
            "input_pairs_checked": len(partitions) ** 2,
            "tiny_graph_control_failures": graph_control_failures,
        },
        "bilinear_probability_control": {
            "left_total": fraction_text(sum(left_weights)),
            "right_total": fraction_text(sum(right_weights)),
            "output_total": fraction_text(sum(composed)),
            "output_nonnegative": all(value >= 0 for value in composed),
            "output_weights": [fraction_text(value) for value in composed],
        },
        "exact_checks": {
            "bell_number_four_is_15": len(partitions) == 15,
            "join_associative": associativity_failures == 0,
            "join_commutative": all(joins[i][j] == joins[j][i] for i in range(len(partitions)) for j in range(len(partitions))),
            "join_idempotent": all(joins[i][i] == i for i in range(len(partitions))),
            "interface_matches_all_tiny_graph_compositions": graph_control_failures == 0,
            "bilinear_mass_preserved": sum(composed) == sum(left_weights) * sum(right_weights) == 1,
        },
        "claim_boundary": {
            "included": "finite RGS join algebra and one declared two-port gluing tensor",
            "excluded": "W5 periodic gluing, planar duality, self-duality, composition-word search, critical manifolds, thresholds, or bounds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("partition-gluing artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_terminal_partition_gluing_algebra",
        "partition_count": len(expected["partition_catalog"]),
        "interface_input_pairs": expected["two_port_interface"]["input_pairs_checked"],
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
