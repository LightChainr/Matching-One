#!/usr/bin/env python3
"""Exact complement involution on bounded canonical terminal-gadget orbits."""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

try:
    from scripts.gadget_graph_canonical import (
        Graph,
        canonical_graph,
        decode_graph,
        graph_orbit_catalog,
        relabel_graph,
    )
    from scripts.terminal_partition_canonical import full_symmetric_group
except ModuleNotFoundError:
    from gadget_graph_canonical import Graph, canonical_graph, decode_graph, graph_orbit_catalog, relabel_graph
    from terminal_partition_canonical import full_symmetric_group


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "gadget_complement_involution_certificate.json"
SCHEMA = "matching-one/gadget-complement-involution-certificate/v1"


def complement_graph(vertex_count: int, edges: Iterable[Sequence[int]]) -> Graph:
    graph = {tuple(sorted(edge)) for edge in edges}
    slots = set(combinations(range(vertex_count), 2))
    if not graph <= slots:
        raise ValueError("graph contains an invalid edge")
    return tuple(sorted(slots - graph))


def complement_commutes_with_relabeling(
    vertex_count: int,
    terminal_count: int,
    edges: Iterable[Sequence[int]],
    vertex_map: Sequence[int],
) -> bool:
    left = complement_graph(
        vertex_count,
        relabel_graph(vertex_count, terminal_count, edges, vertex_map),
    )
    right = relabel_graph(
        vertex_count,
        terminal_count,
        complement_graph(vertex_count, edges),
        vertex_map,
    )
    return left == right


def build_row(terminal_count: int) -> dict[str, Any]:
    if terminal_count not in (3, 4):
        raise ValueError("complement certificate supports three or four terminals")
    vertex_count = terminal_count + 1
    edge_slots = vertex_count * (vertex_count - 1) // 2
    group = full_symmetric_group(terminal_count)
    catalog = graph_orbit_catalog(vertex_count, terminal_count, group)
    mapping = {}
    for encoding in catalog:
        _, _, graph = decode_graph(encoding)
        target = canonical_graph(
            vertex_count,
            terminal_count,
            complement_graph(vertex_count, graph),
            group,
        )
        mapping[encoding] = target
    if set(mapping.values()) != set(mapping):
        raise ArithmeticError("complement map is not a permutation of canonical orbits")
    if any(mapping[mapping[source]] != source for source in mapping):
        raise ArithmeticError("complement map is not an involution")
    if any(
        len(decode_graph(source)[2]) + len(decode_graph(target)[2]) != edge_slots
        for source, target in mapping.items()
    ):
        raise ArithmeticError("complement map does not reverse edge count")

    self_complementary = sorted(source for source, target in mapping.items() if source == target)
    paired = sorted(
        (source, target)
        for source, target in mapping.items()
        if source < target
    )
    edge_histogram = Counter(len(decode_graph(encoding)[2]) for encoding in mapping)
    histogram = [edge_histogram[index] for index in range(edge_slots + 1)]
    if histogram != list(reversed(histogram)):
        raise ArithmeticError("orbit edge-count histogram is not palindromic")
    return {
        "terminal_count": terminal_count,
        "internal_count": 1,
        "edge_slots": edge_slots,
        "canonical_orbits": len(mapping),
        "self_complementary_orbits": len(self_complementary),
        "complement_pairs": len(paired),
        "orbit_accounting": len(self_complementary) + 2 * len(paired),
        "orbit_edge_count_histogram": histogram,
        "self_complementary_encodings": self_complementary,
        "involution_certified": True,
        "edge_reversal_certified": True,
        "palindromic_histogram_certified": True,
    }


def build_artifact() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_complement_involution_certificate",
        "derivation": (
            "complementation commutes with every vertex relabeling, hence descends to the "
            "terminal-symmetry quotient and reverses k edges to M-k"
        ),
        "rows": [build_row(3), build_row(4)],
        "claim_boundary": {
            "included": "complement involution on all one-internal-vertex canonical graph orbits for t in {3,4}",
            "excluded": "connected-filter preservation, probability or planar duality, tilings, critical manifolds, ranking, thresholds, or bounds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    if artifact != build_artifact():
        raise ValueError("complement artifact does not exactly reproduce")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_complement_involution_certificate",
        "rows": [
            {
                "terminal_count": row["terminal_count"],
                "self_complementary_orbits": row["self_complementary_orbits"],
                "complement_pairs": row["complement_pairs"],
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
