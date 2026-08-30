#!/usr/bin/env python3
"""Canonicalize finite simple terminal gadgets under explicit symmetries."""

from __future__ import annotations

import argparse
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence, Tuple

try:
    from scripts.terminal_partition_canonical import full_symmetric_group, validate_group
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from terminal_partition_canonical import full_symmetric_group, validate_group


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "analysis" / "gadget_graph_canonical_manifest.json"
Edge = Tuple[int, int]
Graph = Tuple[Edge, ...]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate_graph(vertex_count: int, terminal_count: int, edges: Iterable[Sequence[int]]) -> Graph:
    """Validate an undirected simple graph with terminals 0..terminal_count-1."""

    _require(type(vertex_count) is int and vertex_count >= 2, "vertex count must be at least two")
    _require(type(terminal_count) is int and 2 <= terminal_count <= vertex_count, "invalid terminal count")
    normalized = []
    for edge in edges:
        _require(len(edge) == 2, "each edge must have two endpoints")
        u, v = edge
        _require(type(u) is int and type(v) is int, "edge endpoints must be integers")
        _require(0 <= u < vertex_count and 0 <= v < vertex_count, "edge endpoint out of range")
        _require(u != v, "self-loops are not simple edges")
        normalized.append((min(u, v), max(u, v)))
    _require(len(normalized) == len(set(normalized)), "parallel/duplicate edges are forbidden")
    return tuple(sorted(normalized))


def _edge_slots(vertex_count: int) -> Tuple[Edge, ...]:
    return tuple(combinations(range(vertex_count), 2))


def encode_graph(vertex_count: int, terminal_count: int, edges: Iterable[Sequence[int]]) -> str:
    graph = set(validate_graph(vertex_count, terminal_count, edges))
    bits = "".join("1" if edge in graph else "0" for edge in _edge_slots(vertex_count))
    return "%d:%d:%s" % (terminal_count, vertex_count, bits)


def decode_graph(encoding: str) -> tuple[int, int, Graph]:
    _require(isinstance(encoding, str), "encoding must be a string")
    fields = encoding.split(":")
    _require(len(fields) == 3, "encoding must contain terminal, vertex, and bit fields")
    try:
        terminal_count, vertex_count = int(fields[0]), int(fields[1])
    except ValueError as exc:
        raise ValueError("encoding counts must be integers") from exc
    slots = _edge_slots(vertex_count)
    bits = fields[2]
    _require(len(bits) == len(slots) and set(bits) <= {"0", "1"}, "invalid adjacency bit field")
    graph = tuple(edge for edge, bit in zip(slots, bits) if bit == "1")
    return terminal_count, vertex_count, validate_graph(vertex_count, terminal_count, graph)


def _vertex_maps(
    vertex_count: int,
    terminal_count: int,
    terminal_group: Iterable[Sequence[int]],
) -> Tuple[Tuple[int, ...], ...]:
    group = validate_group(terminal_group, terminal_count)
    internal = tuple(range(terminal_count, vertex_count))
    maps = []
    for terminal_map in group:
        for moved_internal in permutations(internal):
            mapping = list(terminal_map) + list(moved_internal)
            maps.append(tuple(mapping))
    return tuple(sorted(maps))


def relabel_graph(
    vertex_count: int,
    terminal_count: int,
    edges: Iterable[Sequence[int]],
    vertex_map: Sequence[int],
) -> Graph:
    graph = validate_graph(vertex_count, terminal_count, edges)
    mapping = tuple(vertex_map)
    _require(len(mapping) == vertex_count and set(mapping) == set(range(vertex_count)), "vertex map must be a bijection")
    _require(set(mapping[:terminal_count]) == set(range(terminal_count)), "terminal vertices must map to terminals")
    _require(set(mapping[terminal_count:]) == set(range(terminal_count, vertex_count)), "internal vertices must map internally")
    return validate_graph(vertex_count, terminal_count, ((mapping[u], mapping[v]) for u, v in graph))


def canonical_graph(
    vertex_count: int,
    terminal_count: int,
    edges: Iterable[Sequence[int]],
    terminal_group: Iterable[Sequence[int]],
) -> str:
    graph = validate_graph(vertex_count, terminal_count, edges)
    return min(
        encode_graph(vertex_count, terminal_count, relabel_graph(vertex_count, terminal_count, graph, mapping))
        for mapping in _vertex_maps(vertex_count, terminal_count, terminal_group)
    )


def enumerate_graphs(vertex_count: int, terminal_count: int) -> Tuple[Graph, ...]:
    validate_graph(vertex_count, terminal_count, ())
    slots = _edge_slots(vertex_count)
    return tuple(
        tuple(slots[index] for index in range(len(slots)) if mask & (1 << index))
        for mask in range(1 << len(slots))
    )


def graph_orbit_catalog(
    vertex_count: int,
    terminal_count: int,
    terminal_group: Iterable[Sequence[int]],
) -> Mapping[str, int]:
    counts: dict[str, int] = {}
    group = validate_group(terminal_group, terminal_count)
    for graph in enumerate_graphs(vertex_count, terminal_count):
        key = canonical_graph(vertex_count, terminal_count, graph, group)
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def validate_manifest(manifest: Mapping[str, Any]) -> dict[str, Any]:
    _require(manifest.get("schema") == "matching-one/gadget-graph-canonical/v1", "unknown schema")
    _require(manifest.get("issue") == 13, "wrong issue")
    _require(manifest.get("status") == "graph_encoding_only", "scope status drift")
    contract = manifest.get("encoding_contract", {})
    _require(contract.get("graph_type") == "finite_undirected_simple", "graph type drift")
    _require(contract.get("terminals") == "vertices_0_through_t_minus_1", "terminal convention drift")
    _require(contract.get("edge_order") == "lexicographic_upper_triangle", "edge order drift")
    _require(contract.get("canonicalization") == "lexicographic_minimum_over_explicit_terminal_group_and_all_internal_relabelings", "canonicalization drift")

    audited = {}
    for row in manifest.get("exhaustive_checks", []):
        terminal_count = row.get("terminal_count")
        internal_count = row.get("internal_count")
        _require(terminal_count in (3, 4) and internal_count == 1, "unexpected exhaustive scope")
        vertex_count = terminal_count + internal_count
        catalog = graph_orbit_catalog(vertex_count, terminal_count, full_symmetric_group(terminal_count))
        graph_count = 1 << (vertex_count * (vertex_count - 1) // 2)
        _require(row.get("labeled_simple_graphs") == graph_count, "labeled graph count drift")
        _require(row.get("canonical_orbits") == len(catalog), "canonical orbit count drift")
        _require(sum(catalog.values()) == graph_count, "orbit accounting drift")
        audited[str(terminal_count)] = {"labeled_graphs": graph_count, "canonical_orbits": len(catalog)}
    _require(set(audited) == {"3", "4"}, "both terminal audits are required")
    boundary = manifest.get("claim_boundary", {})
    _require(boundary.get("parent_issue") == "remain open", "parent boundary drift")
    for phrase in ("probability", "planarity", "self-duality", "search ranking"):
        _require(phrase in boundary.get("excluded", ""), "missing excluded boundary: %s" % phrase)
    return {"schema": manifest["schema"], "status": "valid_graph_encoding_only", "audited": audited}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args(argv)
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    print(json.dumps(validate_manifest(manifest), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
