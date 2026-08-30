#!/usr/bin/env python3
"""Exact bounded census of canonical one-internal-vertex terminal gadgets."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from math import comb
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

try:
    from scripts.gadget_graph_canonical import (
        Graph,
        canonical_graph,
        decode_graph,
        enumerate_graphs,
    )
    from scripts.terminal_partition_canonical import full_symmetric_group
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from gadget_graph_canonical import Graph, canonical_graph, decode_graph, enumerate_graphs
    from terminal_partition_canonical import full_symmetric_group


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "bounded_gadget_census.json"
SCHEMA = "matching-one/bounded-gadget-census/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def vertex_degrees(vertex_count: int, edges: Iterable[Sequence[int]]) -> tuple[int, ...]:
    degrees = [0] * vertex_count
    for u, v in edges:
        degrees[u] += 1
        degrees[v] += 1
    return tuple(degrees)


def is_connected(vertex_count: int, edges: Iterable[Sequence[int]]) -> bool:
    adjacency = [set() for _ in range(vertex_count)]
    for u, v in edges:
        adjacency[u].add(v)
        adjacency[v].add(u)
    seen = {0}
    stack = [0]
    while stack:
        vertex = stack.pop()
        for neighbor in adjacency[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return len(seen) == vertex_count


def _histogram(graphs: Iterable[Graph]) -> dict[str, int]:
    return {str(edge_count): count for edge_count, count in sorted(Counter(map(len, graphs)).items())}


def build_census_row(terminal_count: int) -> dict[str, Any]:
    _require(terminal_count in (3, 4), "bounded census supports three or four terminals")
    internal_count = 1
    vertex_count = terminal_count + internal_count
    internal_vertex = terminal_count
    group = full_symmetric_group(terminal_count)
    labeled = enumerate_graphs(vertex_count, terminal_count)
    buckets: dict[str, list[Graph]] = {}
    for graph in labeled:
        key = canonical_graph(vertex_count, terminal_count, graph, group)
        buckets.setdefault(key, []).append(graph)

    representatives = [decode_graph(key)[2] for key in sorted(buckets)]
    connected = [graph for graph in representatives if is_connected(vertex_count, graph)]
    degree3 = [
        graph
        for graph in connected
        if vertex_degrees(vertex_count, graph)[internal_vertex] >= 3
    ]
    orbit_sizes = Counter(len(members) for members in buckets.values())
    return {
        "terminal_count": terminal_count,
        "internal_count": internal_count,
        "vertex_count": vertex_count,
        "edge_slots": comb(vertex_count, 2),
        "labeled_simple_graphs": len(labeled),
        "canonical_orbits": len(buckets),
        "orbit_multiplicity_histogram": {
            str(size): count for size, count in sorted(orbit_sizes.items())
        },
        "orbit_multiplicity_sum": sum(len(members) for members in buckets.values()),
        "all_orbit_edge_count_histogram": _histogram(representatives),
        "connected_carrier_orbits": len(connected),
        "connected_carrier_edge_count_histogram": _histogram(connected),
        "connected_internal_degree_at_least_3_orbits": len(degree3),
        "connected_internal_degree_at_least_3_edge_count_histogram": _histogram(degree3),
    }


def build_artifact() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_bounded_candidate_space_census",
        "symmetry": "full terminal symmetric group; one internal vertex is fixed",
        "rows": [build_census_row(3), build_census_row(4)],
        "filters": {
            "connected_carrier": "all terminal and internal vertices lie in one connected component",
            "connected_internal_degree_at_least_3": (
                "connected carrier and the sole nonterminal vertex has degree at least three; "
                "this removes isolated, leaf, and degree-two cases but is not a general "
                "series-parallel irreducibility certificate"
            ),
        },
        "claim_boundary": {
            "included": "exact census for t in {3,4} with exactly one internal vertex",
            "excluded": (
                "probability polynomials, planarity or periodic tiling, general series-parallel "
                "reduction, self-duality, critical manifolds, ranking, optimization, or a bound"
            ),
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    _require(artifact == build_artifact(), "census artifact does not exactly reproduce")
    _require(artifact.get("schema") == SCHEMA, "unknown census schema")
    _require(artifact.get("issue") == 13, "wrong issue")
    _require(
        artifact.get("claim_boundary", {}).get("parent_issue") == "remain open",
        "parent issue boundary drift",
    )
    return {
        "schema": SCHEMA,
        "status": "valid_exact_bounded_candidate_space_census",
        "rows": [
            {
                "terminal_count": row["terminal_count"],
                "labeled_simple_graphs": row["labeled_simple_graphs"],
                "canonical_orbits": row["canonical_orbits"],
                "connected_carrier_orbits": row["connected_carrier_orbits"],
                "connected_internal_degree_at_least_3_orbits": row[
                    "connected_internal_degree_at_least_3_orbits"
                ],
            }
            for row in artifact["rows"]
        ],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate is not None:
        artifact = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_artifact(artifact), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
