#!/usr/bin/env python3
"""Exact enriched W5 relative-dual state and checkerboard periodic gluing."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from itertools import product
import json
from pathlib import Path
from typing import Iterable, Sequence

try:  # script execution
    from p14_w5_terminal_duality import (
        Q,
        canonical_partition,
        partition_key,
        primal_partition,
        relative_dual_partition,
        spherical_transform,
    )
except ModuleNotFoundError:  # package import in tests
    from .p14_w5_terminal_duality import (
        Q,
        canonical_partition,
        partition_key,
        primal_partition,
        relative_dual_partition,
        spherical_transform,
    )


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "analysis" / "p14_w5_joint_relative_dual_manifest.json"
DEFAULT_OUTPUT = (
    ROOT / "results" / "terminal-reliability" / "p14-w5-joint-relative-dual.json"
)

IndexPartition = tuple[tuple[int, ...], ...]
RelativeState = tuple[tuple[int, ...], IndexPartition]
Vertex = tuple[object, ...]
Edge = tuple[Vertex, Vertex]


def union_partition(vertices: Sequence[int], edges: Iterable[tuple[int, int]]) -> IndexPartition:
    parent = {vertex: vertex for vertex in vertices}

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    def union(left: int, right: int) -> None:
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for left, right in edges:
        union(left, right)
    groups: dict[int, list[int]] = defaultdict(list)
    for vertex in vertices:
        groups[find(vertex)].append(vertex)
    return canonical_partition(groups.values())


def face_partition(bits: Sequence[int]) -> IndexPartition:
    """Connectivity of F_0,...,F_3 through closed primal spokes."""

    return union_partition(
        tuple(range(4)),
        (
            ((index - 1) % 4, index)
            for index in range(4)
            if not bits[4 + index]
        ),
    )


def relative_state(bits: Sequence[int]) -> RelativeState:
    """Connectivity sufficient statistics before the outer arcs are glued."""

    attachments = tuple(index for index in range(4) if not bits[index])
    return attachments, face_partition(bits)


def full_relative_partition(state: RelativeState) -> IndexPartition:
    """Partition of F_0..F_3,B_0..B_3 encoded as integer labels 0..7."""

    attachments, faces = state
    edges = []
    for group in faces:
        edges.extend((group[0], vertex) for vertex in group[1:])
    edges.extend((index, 4 + index) for index in attachments)
    return union_partition(tuple(range(8)), edges)


def labelled_partition_key(partition: IndexPartition) -> str:
    def label(value: int) -> str:
        return f"F{value}" if value < 4 else f"B{value - 4}"

    return "|".join(
        ",".join(label(value) for value in group) for group in partition
    )


def spherical_partition_from_state(state: RelativeState) -> IndexPartition:
    """Glue B_0,...,B_3 to the outer face O and read F_i -> q_-i."""

    attachments, faces = state
    outer = 4
    edges = []
    for group in faces:
        edges.extend((group[0], vertex) for vertex in group[1:])
    edges.extend((outer, index) for index in attachments)
    glued = union_partition(tuple(range(5)), edges)

    parent_group: dict[int, tuple[int, ...]] = {}
    for group in glued:
        for vertex in group:
            parent_group[vertex] = group
    terminal_groups: dict[tuple[int, ...], list[int]] = defaultdict(list)
    for face_index in range(4):
        terminal_groups[parent_group[face_index]].append(Q[(-face_index) % 4])
    return canonical_partition(terminal_groups.values())


def transform_state(state: RelativeState, sign: int, shift: int) -> RelativeState:
    transform = lambda index: (sign * index + shift) % 4
    attachments, faces = state
    return (
        tuple(sorted(transform(index) for index in attachments)),
        canonical_partition(
            tuple(transform(index) for index in group) for group in faces
        ),
    )


def state_key(state: RelativeState) -> str:
    attachments, faces = state
    mask = "".join("1" if index in attachments else "0" for index in range(4))
    return f"A={mask};F={partition_key(faces)}"


def d4_canonical_state(state: RelativeState) -> RelativeState:
    return min(
        (
            transform_state(state, sign, shift)
            for sign in (1, -1)
            for shift in range(4)
        ),
        key=state_key,
    )


def canonical_edge(left: Vertex, right: Vertex) -> Edge:
    return tuple(sorted((left, right), key=repr))  # type: ignore[return-value]


def degree_histogram(vertices: set[Vertex], edges: set[Edge]) -> dict[str, int]:
    degrees = Counter({vertex: 0 for vertex in vertices})
    for left, right in edges:
        degrees[left] += 1
        degrees[right] += 1
    return {str(degree): count for degree, count in sorted(Counter(degrees.values()).items())}


def checkerboard_periodic_embedding(length: int) -> dict[str, object]:
    if length < 4 or length % 2:
        raise ValueError("checkerboard torus requires an even length at least four")

    def mod(value: int) -> int:
        return value % length

    def q(x: int, y: int) -> Vertex:
        return ("q", mod(x), mod(y))

    def h(x: int, y: int) -> Vertex:
        return ("h", mod(x), mod(y))

    def f(x: int, y: int, index: int) -> Vertex:
        return ("f", mod(x), mod(y), index % 4)

    def w(x: int, y: int) -> Vertex:
        return ("w", mod(x), mod(y))

    black = [
        (x, y)
        for y in range(length)
        for x in range(length)
        if (x + y) % 2 == 0
    ]
    white = {
        (x, y)
        for y in range(length)
        for x in range(length)
        if (x + y) % 2 == 1
    }
    primal_vertices = {q(x, y) for y in range(length) for x in range(length)}
    primal_vertices.update(h(x, y) for x, y in black)
    dual_vertices = {w(x, y) for x, y in white}
    dual_vertices.update(f(x, y, index) for x, y in black for index in range(4))

    primal_to_dual: dict[Edge, Edge] = {}
    rim_owners: dict[Edge, tuple[int, int, int]] = {}
    for x, y in black:
        corners = (q(x, y), q(x + 1, y), q(x + 1, y + 1), q(x, y + 1))
        white_neighbors = (
            (mod(x), mod(y - 1)),
            (mod(x + 1), mod(y)),
            (mod(x), mod(y + 1)),
            (mod(x - 1), mod(y)),
        )
        for index in range(4):
            rim = canonical_edge(corners[index], corners[(index + 1) % 4])
            if rim in rim_owners:
                raise AssertionError("a square-grid rim edge has multiple W5 owners")
            rim_owners[rim] = (x, y, index)
            wx, wy = white_neighbors[index]
            primal_to_dual[rim] = canonical_edge(f(x, y, index), w(wx, wy))

            spoke = canonical_edge(h(x, y), corners[index])
            if spoke in primal_to_dual:
                raise AssertionError("duplicate W5 spoke")
            primal_to_dual[spoke] = canonical_edge(
                f(x, y, index - 1), f(x, y, index)
            )

    primal_edges = set(primal_to_dual)
    dual_edges = set(primal_to_dual.values())
    if len(primal_edges) != len(dual_edges):
        raise AssertionError("periodic primal/dual edge correspondence is not bijective")
    if any(vertex not in primal_vertices for edge in primal_edges for vertex in edge):
        raise AssertionError("primal edge endpoint missing from vertex set")
    if any(vertex not in dual_vertices for edge in dual_edges for vertex in edge):
        raise AssertionError("dual edge endpoint missing from vertex set")
    if len(rim_owners) != 2 * length * length:
        raise AssertionError("checkerboard W5 cells do not own every grid edge exactly once")

    primal_hist = degree_histogram(primal_vertices, primal_edges)
    dual_hist = degree_histogram(dual_vertices, dual_edges)
    expected_primal = {"4": length * length // 2, "6": length * length}
    expected_dual = {"3": 2 * length * length, "4": length * length // 2}
    if primal_hist != expected_primal or dual_hist != expected_dual:
        raise AssertionError("periodic checkerboard degree census changed")
    euler = len(primal_vertices) - len(primal_edges) + len(dual_vertices)
    if euler != 0:
        raise AssertionError("periodic checkerboard incidence violates torus Euler identity")

    return {
        "length": length,
        "black_w5_cells": len(black),
        "white_outer_faces": len(white),
        "primal": {
            "vertices": len(primal_vertices),
            "edges": len(primal_edges),
            "degree_histogram": primal_hist,
        },
        "disk_relative_dual": {
            "vertices": len(dual_vertices),
            "edges": len(dual_edges),
            "degree_histogram": dual_hist,
        },
        "edge_bijection_count": len(primal_to_dual),
        "every_grid_rim_owned_once": True,
        "torus_euler_characteristic": euler,
        "per_black_cell_density": {
            "primal_vertices": "3",
            "edges": "8",
            "dual_vertices_or_primal_faces": "5",
        },
        "graph_self_isomorphism_obstructed": (
            len(primal_vertices) != len(dual_vertices) or primal_hist != dual_hist
        ),
        "parameter_duality": {
            "primal_rim_r": "dual_leaf_1-r",
            "primal_spoke_s": "dual_internal_cycle_1-s",
        },
    }


def build_report(manifest: dict[str, object]) -> dict[str, object]:
    grouped: dict[RelativeState, list[tuple[int, ...]]] = defaultdict(list)
    primal_outputs: dict[RelativeState, set[str]] = defaultdict(set)
    spherical_outputs: dict[RelativeState, set[str]] = defaultdict(set)
    relative_outputs: dict[RelativeState, set[str]] = defaultdict(set)
    terminal_pair_outputs: dict[tuple[str, str], set[str]] = defaultdict(set)

    for bits in product((0, 1), repeat=8):
        state = relative_state(bits)
        primal_key = partition_key(primal_partition(bits))
        relative_key = partition_key(relative_dual_partition(bits))
        spherical_key = partition_key(primal_partition(spherical_transform(bits)))
        reconstructed = partition_key(spherical_partition_from_state(state))
        if reconstructed != spherical_key:
            raise AssertionError("relative-dual gluing did not reconstruct spherical output")
        grouped[state].append(bits)
        primal_outputs[state].add(primal_key)
        relative_outputs[state].add(relative_key)
        spherical_outputs[state].add(spherical_key)
        terminal_pair_outputs[(primal_key, relative_key)].add(spherical_key)

    if any(len(outputs) != 1 for outputs in primal_outputs.values()):
        raise AssertionError("enriched relative state does not determine primal partition")
    if any(len(outputs) != 1 for outputs in spherical_outputs.values()):
        raise AssertionError("enriched relative state does not determine spherical output")
    if any(len(outputs) != 1 for outputs in relative_outputs.values()):
        raise AssertionError("enriched relative state does not determine relative boundary output")

    rows = []
    for state in sorted(grouped, key=state_key):
        configurations = grouped[state]
        rows.append(
            {
                "state": state_key(state),
                "attachment_mask_closed_rim": state_key(state).split(";")[0].split("=")[1],
                "face_partition_closed_spokes": partition_key(state[1]),
                "full_relative_partition": labelled_partition_key(full_relative_partition(state)),
                "primal_partition": next(iter(primal_outputs[state])),
                "relative_boundary_partition": next(iter(relative_outputs[state])),
                "spherical_partition_after_outer_gluing": next(iter(spherical_outputs[state])),
                "configuration_multiplicity": len(configurations),
                "d4_orbit_representative": state_key(d4_canonical_state(state)),
            }
        )

    multiplicities = Counter(row["configuration_multiplicity"] for row in rows)
    ambiguous_terminal_pairs = {
        f"P={primal};D={relative}": sorted(outputs)
        for (primal, relative), outputs in sorted(terminal_pair_outputs.items())
        if len(outputs) > 1
    }
    periodic = checkerboard_periodic_embedding(int(manifest["torus_length"]))
    return {
        "schema": manifest["schema"],
        "parent_certificate": manifest["parent_certificate"],
        "finite_cell": {
            "configuration_count": 256,
            "labelled_enriched_state_count": len(rows),
            "d4_orbit_count": len({row["d4_orbit_representative"] for row in rows}),
            "configuration_multiplicity_histogram": {
                str(key): value for key, value in sorted(multiplicities.items())
            },
            "terminal_pair_state_count": len(terminal_pair_outputs),
            "terminal_pair_ambiguous_for_spherical_output": len(ambiguous_terminal_pairs),
            "terminal_pair_ambiguity_examples": dict(list(ambiguous_terminal_pairs.items())[:8]),
            "enriched_state_determines_primal_partition": True,
            "enriched_state_determines_relative_boundary_partition": True,
            "enriched_state_determines_spherical_partition_after_outer_gluing": True,
            "state_rows": rows,
        },
        "periodic_checkerboard_embedding": periodic,
        "decision": (
            "PASS finite-cell closure: the relative-dual face partition plus boundary-attachment mask "
            "is the missing connectivity state and reconstructs both primal and outer-glued spherical outputs. "
            "FAIL graph self-duality in the natural periodic checkerboard gluing: the primal has 3 vertices "
            "per W5 cell and degree types 4/6, whereas the disk-relative dual has 5 vertices and degree types 3/4."
        ),
        "next_object": (
            "use the explicit checkerboard primal/dual edge bijection to seek a stochastic comparison or "
            "local transformation between the nonisomorphic periodic graphs; do not return to a scalar W5 balance root"
        ),
        "claim_boundary": manifest["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    report = build_report(manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
