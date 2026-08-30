#!/usr/bin/env python3
"""Exact bounded planarity certificate by orientable rotation systems."""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import permutations, product
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

try:
    from scripts.bounded_gadget_census import is_connected, vertex_degrees
    from scripts.gadget_graph_canonical import Graph, decode_graph, graph_orbit_catalog
    from scripts.terminal_partition_canonical import full_symmetric_group
except ModuleNotFoundError:
    from bounded_gadget_census import is_connected, vertex_degrees
    from gadget_graph_canonical import Graph, decode_graph, graph_orbit_catalog
    from terminal_partition_canonical import full_symmetric_group


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "bounded_gadget_planarity_certificate.json"
SCHEMA = "matching-one/bounded-gadget-planarity-certificate/v1"


def validate_simple_graph(vertex_count: int, edges: Iterable[Sequence[int]]) -> Graph:
    normalized = []
    for edge in edges:
        if len(edge) != 2:
            raise ValueError("edge must have two endpoints")
        u, v = edge
        if not (0 <= u < vertex_count and 0 <= v < vertex_count):
            raise ValueError("edge endpoint out of range")
        if u == v:
            raise ValueError("self-loop is not a simple edge")
        normalized.append(tuple(sorted((u, v))))
    if len(normalized) != len(set(normalized)):
        raise ValueError("duplicate edge is not simple")
    return tuple(sorted(normalized))


def connected_components(vertex_count: int, graph: Graph) -> tuple[tuple[int, ...], ...]:
    adjacency = [set() for _ in range(vertex_count)]
    for u, v in graph:
        adjacency[u].add(v)
        adjacency[v].add(u)
    unseen = set(range(vertex_count))
    components = []
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
        components.append(tuple(sorted(component)))
    return tuple(sorted(components))


def cyclic_orders(neighbors: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    ordered = tuple(sorted(neighbors))
    if len(ordered) <= 1:
        return (ordered,)
    first = ordered[0]
    return tuple((first,) + tail for tail in permutations(ordered[1:]))


def face_count(graph: Graph, rotations: Mapping[int, Sequence[int]]) -> int:
    darts = tuple((u, v) for edge in graph for u, v in (edge, tuple(reversed(edge))))
    if not darts:
        return 1
    successor = {}
    for vertex, order_source in rotations.items():
        order = tuple(order_source)
        if order:
            for index, neighbor in enumerate(order):
                successor[(vertex, neighbor)] = (vertex, order[(index + 1) % len(order)])
    face_permutation = {}
    for u, v in darts:
        face_permutation[(u, v)] = successor[(v, u)]
    unseen = set(darts)
    faces = 0
    while unseen:
        faces += 1
        dart = unseen.pop()
        current = face_permutation[dart]
        while current != dart:
            unseen.remove(current)
            current = face_permutation[current]
    return faces


def connected_minimum_orientable_genus(
    vertices: Sequence[int], graph: Graph
) -> tuple[int, int, dict[str, list[int]]]:
    vertices = tuple(vertices)
    if len(vertices) == 1:
        return 0, 1, {str(vertices[0]): []}
    adjacency = {vertex: [] for vertex in vertices}
    for u, v in graph:
        adjacency[u].append(v)
        adjacency[v].append(u)
    choices = [cyclic_orders(adjacency[vertex]) for vertex in vertices]
    checked = 0
    best_genus = None
    best_rotation = None
    for orders in product(*choices):
        checked += 1
        rotations = dict(zip(vertices, orders))
        faces = face_count(graph, rotations)
        euler_characteristic = len(vertices) - len(graph) + faces
        numerator = 2 - euler_characteristic
        if numerator < 0 or numerator % 2:
            raise ArithmeticError("rotation system produced invalid orientable genus")
        genus = numerator // 2
        if best_genus is None or genus < best_genus:
            best_genus = genus
            best_rotation = rotations
        if genus == 0:
            break
    assert best_genus is not None and best_rotation is not None
    return best_genus, checked, {
        str(vertex): list(best_rotation[vertex]) for vertex in vertices
    }


def minimum_orientable_genus(
    vertex_count: int, edges: Iterable[Sequence[int]]
) -> tuple[int, int, list[dict[str, list[int]]]]:
    graph = validate_simple_graph(vertex_count, edges)
    total_genus = 0
    total_checked = 0
    witnesses = []
    for component in connected_components(vertex_count, graph):
        members = set(component)
        component_graph = tuple(edge for edge in graph if edge[0] in members)
        genus, checked, witness = connected_minimum_orientable_genus(
            component, component_graph
        )
        total_genus += genus
        total_checked += checked
        witnesses.append(witness)
    return total_genus, total_checked, witnesses


def summarize(graphs: Iterable[tuple[str, Graph]], terminal_count: int) -> dict[str, Any]:
    graph_list = tuple(graphs)
    genus_histogram = Counter()
    nonplanar = []
    exhaustive_nonplanar_rotation_systems = 0
    for encoding, graph in graph_list:
        genus, checked, _ = minimum_orientable_genus(terminal_count + 1, graph)
        genus_histogram[genus] += 1
        if genus:
            nonplanar.append({"encoding": encoding, "minimum_orientable_genus": genus, "rotation_systems_checked": checked})
            exhaustive_nonplanar_rotation_systems += checked
    return {
        "orbits": len(graph_list),
        "minimum_orientable_genus_histogram": {
            str(genus): count for genus, count in sorted(genus_histogram.items())
        },
        "planar_orbits": genus_histogram[0],
        "nonplanar_orbits": len(graph_list) - genus_histogram[0],
        "nonplanar_witnesses": nonplanar,
        "exhaustive_nonplanar_rotation_systems": exhaustive_nonplanar_rotation_systems,
    }


def build_row(terminal_count: int) -> dict[str, Any]:
    if terminal_count not in (3, 4):
        raise ValueError("bounded planarity certificate supports three or four terminals")
    vertex_count = terminal_count + 1
    catalog = graph_orbit_catalog(vertex_count, terminal_count, full_symmetric_group(terminal_count))
    representatives = [(encoding, decode_graph(encoding)[2]) for encoding in catalog]
    connected = [(encoding, graph) for encoding, graph in representatives if is_connected(vertex_count, graph)]
    degree3 = [
        (encoding, graph)
        for encoding, graph in connected
        if vertex_degrees(vertex_count, graph)[terminal_count] >= 3
    ]
    return {
        "terminal_count": terminal_count,
        "internal_count": 1,
        "all_canonical_orbits": summarize(representatives, terminal_count),
        "connected_carrier": summarize(connected, terminal_count),
        "connected_internal_degree_at_least_3": summarize(degree3, terminal_count),
    }


def build_artifact() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_bounded_planarity_certificate",
        "method": "exhaust orientable cyclic neighbor orders; compute faces and genus from V-E+F=2-2g",
        "rows": [build_row(3), build_row(4)],
        "claim_boundary": {
            "included": "minimum orientable genus for every bounded one-internal-vertex canonical graph orbit",
            "excluded": "periodic tiling, reliability, planar duality, self-duality, critical manifolds, ranking, thresholds, or bounds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    if artifact != build_artifact():
        raise ValueError("planarity artifact does not exactly reproduce")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_bounded_planarity_certificate",
        "rows": [
            {
                "terminal_count": row["terminal_count"],
                "planar_orbits": row["all_canonical_orbits"]["planar_orbits"],
                "nonplanar_orbits": row["all_canonical_orbits"]["nonplanar_orbits"],
            }
            for row in artifact["rows"]
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
