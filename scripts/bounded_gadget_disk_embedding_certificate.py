#!/usr/bin/env python3
"""Exact terminal-order disk-embedding certificate for bounded gadgets."""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

try:
    from scripts.bounded_gadget_census import is_connected, vertex_degrees
    from scripts.bounded_gadget_planarity_certificate import cyclic_orders, validate_simple_graph
    from scripts.gadget_graph_canonical import Graph, decode_graph, graph_orbit_catalog
    from scripts.terminal_partition_canonical import full_symmetric_group
except ModuleNotFoundError:
    from bounded_gadget_census import is_connected, vertex_degrees
    from bounded_gadget_planarity_certificate import cyclic_orders, validate_simple_graph
    from gadget_graph_canonical import Graph, decode_graph, graph_orbit_catalog
    from terminal_partition_canonical import full_symmetric_group


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "bounded_gadget_disk_embedding_certificate.json"
SCHEMA = "matching-one/bounded-gadget-disk-embedding/v1"


def normalize_terminal_order(vertex_count: int, order: Sequence[int]) -> tuple[int, ...]:
    if isinstance(order, (str, bytes)) or not isinstance(order, Sequence) or len(order) < 3:
        raise ValueError("terminal order must contain at least three vertices")
    normalized = tuple(order)
    if any(isinstance(vertex, bool) or not isinstance(vertex, int) for vertex in normalized):
        raise TypeError("terminal order entries must be integers")
    if len(set(normalized)) != len(normalized):
        raise ValueError("terminal order must not repeat vertices")
    if any(vertex < 0 or vertex >= vertex_count for vertex in normalized):
        raise ValueError("terminal order vertex out of range")
    return normalized


def rotation_faces(graph: Graph, rotations: Mapping[int, Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    darts = tuple((u, v) for edge in graph for u, v in (edge, tuple(reversed(edge))))
    if not darts:
        return (tuple(),)
    successor = {}
    for vertex, order_source in rotations.items():
        order = tuple(order_source)
        for index, neighbor in enumerate(order):
            successor[(vertex, neighbor)] = (vertex, order[(index + 1) % len(order)])
    face_permutation = {(u, v): successor[(v, u)] for u, v in darts}
    unseen = set(darts)
    faces = []
    while unseen:
        start = min(unseen)
        current = start
        face = []
        while True:
            unseen.remove(current)
            face.append(current[0])
            current = face_permutation[current]
            if current == start:
                break
        faces.append(tuple(face))
    return tuple(sorted(faces))


def _ordered_subsequence_in_one_turn(face: Sequence[int], order: Sequence[int]) -> bool:
    if not face:
        return False
    doubled = tuple(face) + tuple(face)
    length = len(face)
    for start in range(length):
        if doubled[start] != order[0]:
            continue
        cursor = start + 1
        valid = True
        for terminal in order[1:]:
            while cursor < start + length and doubled[cursor] != terminal:
                cursor += 1
            if cursor >= start + length:
                valid = False
                break
            cursor += 1
        if valid:
            return True
    return False


def face_has_terminal_order(face: Sequence[int], order: Sequence[int]) -> bool:
    normalized = tuple(order)
    return _ordered_subsequence_in_one_turn(face, normalized) or _ordered_subsequence_in_one_turn(
        face, (normalized[0],) + tuple(reversed(normalized[1:]))
    )


def disk_embedding(
    vertex_count: int, edges: Iterable[Sequence[int]], terminal_order: Sequence[int]
) -> dict[str, Any]:
    graph = validate_simple_graph(vertex_count, edges)
    order = normalize_terminal_order(vertex_count, terminal_order)
    if not is_connected(vertex_count, graph):
        raise ValueError("disk-embedding search requires a connected graph")
    adjacency = {vertex: [] for vertex in range(vertex_count)}
    for u, v in graph:
        adjacency[u].append(v)
        adjacency[v].append(u)
    choices = [cyclic_orders(adjacency[vertex]) for vertex in range(vertex_count)]
    checked = 0
    planar_rotation_systems = 0
    for orders in product(*choices):
        checked += 1
        rotations = dict(zip(range(vertex_count), orders))
        faces = rotation_faces(graph, rotations)
        euler_characteristic = vertex_count - len(graph) + len(faces)
        numerator = 2 - euler_characteristic
        if numerator < 0 or numerator % 2:
            raise ArithmeticError("rotation system produced invalid orientable genus")
        if numerator:
            continue
        planar_rotation_systems += 1
        for face in faces:
            if face_has_terminal_order(face, order):
                return {
                    "disk_planar": True,
                    "rotation_systems_checked": checked,
                    "planar_rotation_systems_checked": planar_rotation_systems,
                    "terminal_order": list(order),
                    "witness_face": list(face),
                    "witness_rotation": {str(vertex): list(rotations[vertex]) for vertex in range(vertex_count)},
                }
    return {
        "disk_planar": False,
        "rotation_systems_checked": checked,
        "planar_rotation_systems_checked": planar_rotation_systems,
        "terminal_order": list(order),
        "witness_face": None,
        "witness_rotation": None,
    }


def canonical_cyclic_orders(terminal_count: int) -> tuple[tuple[int, ...], ...]:
    if terminal_count == 3:
        return ((0, 1, 2),)
    if terminal_count == 4:
        return ((0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 1, 3))
    raise ValueError("only three- and four-terminal cyclic orders are supported")


def terminal_partition(vertex_count: int, graph: Graph, terminal_count: int) -> tuple[tuple[int, ...], ...]:
    adjacency = [set() for _ in range(vertex_count)]
    for u, v in graph:
        adjacency[u].add(v)
        adjacency[v].add(u)
    unseen = set(range(vertex_count))
    blocks = []
    while unseen:
        start = unseen.pop()
        component = {start}
        stack = [start]
        while stack:
            vertex = stack.pop()
            for neighbor in adjacency[vertex] & unseen:
                unseen.remove(neighbor)
                component.add(neighbor)
                stack.append(neighbor)
        terminal_block = tuple(sorted(vertex for vertex in component if vertex < terminal_count))
        if terminal_block:
            blocks.append(terminal_block)
    return tuple(sorted(blocks))


def crossing_partition(order: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    if len(order) != 4:
        raise ValueError("crossing partition requires four terminals")
    return tuple(sorted((tuple(sorted((order[0], order[2]))), tuple(sorted((order[1], order[3]))))))


def crossing_subgraph_count(
    vertex_count: int, graph: Graph, terminal_order: Sequence[int]
) -> tuple[int, int]:
    order = normalize_terminal_order(vertex_count, terminal_order)
    forbidden = crossing_partition(order)
    count = 0
    for mask in range(1 << len(graph)):
        subgraph = tuple(edge for index, edge in enumerate(graph) if mask & (1 << index))
        if terminal_partition(vertex_count, subgraph, len(order)) == forbidden:
            count += 1
    return count, 1 << len(graph)


def build_row(terminal_count: int) -> dict[str, Any]:
    vertex_count = terminal_count + 1
    catalog = graph_orbit_catalog(vertex_count, terminal_count, full_symmetric_group(terminal_count))
    representatives = [(encoding, decode_graph(encoding)[2]) for encoding in catalog]
    connected = [(encoding, graph) for encoding, graph in representatives if is_connected(vertex_count, graph)]
    order_classes = canonical_cyclic_orders(terminal_count)
    records = []
    crossing_states_checked = 0
    for encoding, graph in connected:
        order_results = []
        for order in order_classes:
            result = disk_embedding(vertex_count, graph, order)
            if terminal_count == 4 and result["disk_planar"]:
                crossings, checked = crossing_subgraph_count(vertex_count, graph, order)
                crossing_states_checked += checked
                result["crossing_partition"] = [list(block) for block in crossing_partition(order)]
                result["crossing_subgraphs"] = crossings
                result["edge_subgraphs_checked"] = checked
            order_results.append(result)
        records.append({
            "encoding": encoding,
            "edge_count": len(graph),
            "internal_degree": vertex_degrees(vertex_count, graph)[terminal_count],
            "disk_planar_order_classes": sum(result["disk_planar"] for result in order_results),
            "orders": order_results,
        })
    class_histogram = Counter(record["disk_planar_order_classes"] for record in records)
    degree3 = [record for record in records if record["internal_degree"] >= 3]
    return {
        "terminal_count": terminal_count,
        "internal_count": 1,
        "connected_orbits": len(records),
        "cyclic_order_classes_per_orbit": len(order_classes),
        "disk_planar_order_class_histogram": {str(key): value for key, value in sorted(class_histogram.items())},
        "disk_planar_orbit_order_pairs": sum(record["disk_planar_order_classes"] for record in records),
        "non_disk_orbit_order_pairs": len(records) * len(order_classes) - sum(record["disk_planar_order_classes"] for record in records),
        "degree_at_least_3": {
            "orbits": len(degree3),
            "disk_planar_orbit_order_pairs": sum(record["disk_planar_order_classes"] for record in degree3),
        },
        "crossing_edge_subgraphs_checked": crossing_states_checked,
        "records": records,
    }


def build_artifact() -> dict[str, Any]:
    rows = [build_row(3), build_row(4)]
    crossing_failures = sum(
        result.get("crossing_subgraphs", 0)
        for row in rows
        for record in row["records"]
        for result in record["orders"]
    )
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_terminal_order_disk_embedding_certificate",
        "method": "exhaust orientable rotation systems and require one face to contain terminals in the declared cyclic order",
        "rows": rows,
        "exact_checks": {
            "all_disk_witnesses_have_zero_crossing_subgraphs": crossing_failures == 0,
            "three_terminal_connected_orbits_complete": rows[0]["connected_orbits"] == 11,
            "four_terminal_connected_orbits_complete": rows[1]["connected_orbits"] == 58,
            "four_terminal_has_three_cyclic_order_classes": rows[1]["cyclic_order_classes_per_orbit"] == 3,
        },
        "claim_boundary": {
            "included": "bounded connected one-internal-vertex disk embeddings with declared terminal cyclic order",
            "excluded": "planar duality, self-duality, critical manifolds, periodic tilings, candidate ranking, thresholds, or bounds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("disk-embedding artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_terminal_order_disk_embedding_certificate",
        "rows": [
            {
                "terminal_count": row["terminal_count"],
                "connected_orbits": row["connected_orbits"],
                "disk_planar_orbit_order_pairs": row["disk_planar_orbit_order_pairs"],
            }
            for row in expected["rows"]
        ],
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
