#!/usr/bin/env python3
"""Exact tiny oracle for the Riordan--Walters open-boundary block event.

The domain is two adjacent ``s by s`` site cells.  A block bond succeeds when
each half has a unique largest open cluster and those two clusters are joined
in the full ``2s by s`` rectangle.  Both square nearest-neighbour adjacency and
the square site-matching graph (nearest neighbours plus both diagonals of each
unit face) are supported.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Optional, Sequence, Set, Tuple


Vertex = Tuple[int, int]
Edge = Tuple[Vertex, Vertex]
GRAPHS = ("square", "matching")


@dataclass(frozen=True)
class EventDecision:
    success: bool
    reason: str
    left_largest: Tuple[Vertex, ...]
    right_largest: Tuple[Vertex, ...]


def ordered_vertices(s: int) -> list[Vertex]:
    if s < 1:
        raise ValueError("s must be positive")
    return [(x, y) for y in range(s) for x in range(2 * s)]


def graph_edges(s: int, graph: str) -> list[Edge]:
    if graph not in GRAPHS:
        raise ValueError("unknown graph: %s" % graph)
    width = 2 * s
    height = s
    steps = [(1, 0), (0, 1)]
    if graph == "matching":
        steps.extend([(1, 1), (1, -1)])
    edges = []
    for x in range(width):
        for y in range(height):
            for dx, dy in steps:
                other = (x + dx, y + dy)
                if 0 <= other[0] < width and 0 <= other[1] < height:
                    edges.append(((x, y), other))
    return edges


def open_vertices(mask: int, vertices: Sequence[Vertex]) -> Set[Vertex]:
    if mask < 0 or mask >= 1 << len(vertices):
        raise ValueError("mask is outside the rectangle")
    return {vertex for index, vertex in enumerate(vertices) if mask & (1 << index)}


def components(open_set: Set[Vertex], allowed: Set[Vertex], edges: Sequence[Edge]) -> list[Set[Vertex]]:
    active = open_set & allowed
    adjacency = {vertex: [] for vertex in active}
    for first, second in edges:
        if first in active and second in active:
            adjacency[first].append(second)
            adjacency[second].append(first)

    unseen = set(active)
    result = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        queue = deque([start])
        component = {start}
        while queue:
            current = queue.popleft()
            for neighbor in adjacency[current]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    component.add(neighbor)
                    queue.append(neighbor)
        result.append(component)
    result.sort(key=lambda part: (-len(part), tuple(sorted(part))))
    return result


def unique_largest(parts: Sequence[Set[Vertex]]) -> Tuple[Optional[Set[Vertex]], str]:
    if not parts:
        return None, "empty"
    largest_size = max(len(part) for part in parts)
    largest = [part for part in parts if len(part) == largest_size]
    if len(largest) != 1:
        return None, "tie"
    return largest[0], "unique"


def connected_in_union(
    left_cluster: Set[Vertex],
    right_cluster: Set[Vertex],
    union_parts: Sequence[Set[Vertex]],
) -> bool:
    return any(
        bool(part & left_cluster) and bool(part & right_cluster)
        for part in union_parts
    )


def decide_event(mask: int, s: int, graph: str) -> EventDecision:
    vertices = ordered_vertices(s)
    edges = graph_edges(s, graph)
    opened = open_vertices(mask, vertices)
    left = {vertex for vertex in vertices if vertex[0] < s}
    right = set(vertices) - left

    left_largest, left_status = unique_largest(components(opened, left, edges))
    if left_largest is None:
        return EventDecision(False, "left_%s" % left_status, (), ())
    right_largest, right_status = unique_largest(components(opened, right, edges))
    if right_largest is None:
        return EventDecision(
            False,
            "right_%s" % right_status,
            tuple(sorted(left_largest)),
            (),
        )

    union_parts = components(opened, set(vertices), edges)
    success = connected_in_union(left_largest, right_largest, union_parts)
    return EventDecision(
        success,
        "success" if success else "largest_clusters_disconnected",
        tuple(sorted(left_largest)),
        tuple(sorted(right_largest)),
    )


def enumerate_reliability(s: int, graph: str) -> dict[str, Any]:
    vertices = ordered_vertices(s)
    site_count = len(vertices)
    success_by_occupied = [0] * (site_count + 1)
    reasons: Counter[str] = Counter()
    successful_masks = []
    for mask in range(1 << site_count):
        decision = decide_event(mask, s, graph)
        reasons[decision.reason] += 1
        if decision.success:
            # ``int.bit_count`` starts in Python 3.10; CI still supports 3.9.
            occupied = bin(mask).count("1")
            success_by_occupied[occupied] += 1
            successful_masks.append(mask)

    total_successes = sum(success_by_occupied)
    probability_half = Fraction(total_successes, 1 << site_count)
    return {
        "s": s,
        "graph": graph,
        "site_count": site_count,
        "edge_count": len(graph_edges(s, graph)),
        "configuration_count": 1 << site_count,
        "success_count": total_successes,
        "success_by_occupied": success_by_occupied,
        "reliability_polynomial": (
            "sum_k c[k]*p^k*(1-p)^(%d-k)" % site_count
        ),
        "probability_at_half": "%d/%d" % (probability_half.numerator, probability_half.denominator),
        "reason_counts": dict(sorted(reasons.items())),
        "successful_masks": successful_masks,
    }


def single_site_addition_counterexamples(s: int, graph: str) -> list[dict[str, int]]:
    """Find success->failure transitions after opening one additional site."""

    vertices = ordered_vertices(s)
    examples = []
    for mask in range(1 << len(vertices)):
        if not decide_event(mask, s, graph).success:
            continue
        for index in range(len(vertices)):
            if mask & (1 << index):
                continue
            enlarged = mask | (1 << index)
            if not decide_event(enlarged, s, graph).success:
                examples.append({"mask": mask, "opened_index": index, "enlarged_mask": enlarged})
    return examples


def build_artifact() -> dict[str, Any]:
    cases = {
        "%s_s%d" % (graph, s): enumerate_reliability(s, graph)
        for graph in GRAPHS
        for s in (1, 2)
    }
    monotonicity = {
        graph: single_site_addition_counterexamples(2, graph)
        for graph in GRAPHS
    }

    assert cases["square_s1"]["success_by_occupied"] == [0, 0, 1]
    assert cases["matching_s1"]["success_by_occupied"] == [0, 0, 1]
    assert cases["square_s2"]["configuration_count"] == 256
    assert cases["matching_s2"]["configuration_count"] == 256

    return {
        "schema": "matching-one/pc-block-event-exact/v1",
        "issue": 112,
        "status": "tiny_exact_semantics_oracle",
        "vertex_order": "row-major (x,y), y outer, x inner, width=2s, height=s",
        "cells": {
            "left": "0<=x<s",
            "right": "s<=x<2s",
            "boundary": "open rectangle; no periodic identifications",
        },
        "graphs": {
            "square": "nearest-neighbour horizontal/vertical edges",
            "matching": "square edges plus both diagonals of every unit face",
        },
        "event": (
            "each half has a unique largest nonempty open cluster and those two selected "
            "clusters are connected in the full rectangle"
        ),
        "cases": cases,
        "single_site_addition_counterexamples_s2": monotonicity,
        "claim_boundary": {
            "proved": "exact event decisions and reliability coefficients for s=1,2",
            "not_proved": (
                "production-scale performance, event probabilities near pc, statistical independence "
                "of any future sampler, or a new critical-probability bound"
            ),
        },
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    lines = [
        "# Tiny exact open-boundary block-event oracle",
        "",
        "The rectangle consists of adjacent `s x s` cells with open boundaries. A success requires",
        "one unique largest nonempty open cluster in each half and connection of those selected",
        "clusters in the full `2s x s` rectangle.",
        "",
        "## Exact enumeration",
        "",
        "| graph | s | sites | edges | successes/configurations | P(E) at p=1/2 | coefficients c[k] |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for key in ("square_s1", "matching_s1", "square_s2", "matching_s2"):
        case = artifact["cases"][key]
        lines.append(
            "| %s | %d | %d | %d | %d/%d | `%s` | `%s` |"
            % (
                case["graph"],
                case["s"],
                case["site_count"],
                case["edge_count"],
                case["success_count"],
                case["configuration_count"],
                case["probability_at_half"],
                case["success_by_occupied"],
            )
        )
    lines.extend(
        [
            "",
            "The coefficient vector defines the exact reliability polynomial",
            "`sum_k c[k] p^k (1-p)^(2s^2-k)`.",
            "",
            "## Frozen semantics",
            "",
            "- an empty half fails because it has no largest open cluster;",
            "- a tie for largest cluster fails;",
            "- largest clusters are chosen inside each half, then connectivity is tested in the union;",
            "- the matching graph adds both diagonals of each unit square;",
            "- no torus wrap or boundary identification is used.",
            "",
            "## Boundary",
            "",
            "This oracle certifies event semantics only at `s=1,2`. It does not estimate a block-event",
            "probability near `pc`, validate a random sampler, or produce a critical-probability bound.",
            "",
            "## Reproduction",
            "",
            "```bash",
            "python scripts/pc_block_event_exact.py --format json",
            "python scripts/pc_block_event_exact.py --format markdown",
            "python -m unittest tests.test_pc_block_event_exact",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = build_artifact()
    rendered = (
        json.dumps(artifact, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(artifact)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
