#!/usr/bin/env python3
"""Linear-time evaluator for the frozen open-boundary two-cell event.

This module evaluates a deterministic site configuration.  It deliberately
does not provide a random sampler: PRNG output alone would not establish the
independence assumptions required by a rigorous confidence statement.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence, Tuple

import pc_block_event_exact as exact


GRAPHS = ("square", "matching")


class DisjointSet:
    def __init__(self, active: Sequence[bool]) -> None:
        self.parent = [index if value else -1 for index, value in enumerate(active)]
        self.size = [1 if value else 0 for value in active]

    def find(self, item: int) -> int:
        parent = self.parent[item]
        if parent < 0:
            raise ValueError("find called on a closed site")
        while parent != self.parent[parent]:
            parent = self.parent[parent]
        while item != parent:
            next_item = self.parent[item]
            self.parent[item] = parent
            item = next_item
        return parent

    def union(self, first: int, second: int) -> None:
        root_first = self.find(first)
        root_second = self.find(second)
        if root_first == root_second:
            return
        if self.size[root_first] < self.size[root_second]:
            root_first, root_second = root_second, root_first
        self.parent[root_second] = root_first
        self.size[root_first] += self.size[root_second]


@dataclass(frozen=True)
class LinearDecision:
    success: bool
    reason: str
    left_largest: Tuple[int, ...]
    right_largest: Tuple[int, ...]
    full_edge_checks: int
    half_edge_checks: int


def site_count(s: int) -> int:
    if s < 1:
        raise ValueError("s must be positive")
    return 2 * s * s


def flags_from_mask(mask: int, count: int) -> Tuple[bool, ...]:
    if count < 0 or mask < 0 or mask >= 1 << count:
        raise ValueError("mask is outside the requested site count")
    return tuple(bool(mask & (1 << index)) for index in range(count))


def _previous_neighbors(x: int, y: int, width: int, graph: str) -> Iterable[Tuple[int, int]]:
    steps = ((-1, 0), (0, -1))
    if graph == "matching":
        steps += ((-1, -1), (1, -1))
    for dx, dy in steps:
        other_x = x + dx
        other_y = y + dy
        if 0 <= other_x < width and other_y >= 0:
            yield other_x, other_y


def _unique_largest(
    dsu: DisjointSet,
    flags: Sequence[bool],
    indexes: Sequence[int],
) -> Tuple[Tuple[int, ...], str]:
    roots = {dsu.find(index) for index in indexes if flags[index]}
    if not roots:
        return (), "empty"
    maximum = max(dsu.size[root] for root in roots)
    largest_roots = [root for root in roots if dsu.size[root] == maximum]
    if len(largest_roots) != 1:
        return (), "tie"
    selected_root = largest_roots[0]
    selected = tuple(index for index in indexes if flags[index] and dsu.find(index) == selected_root)
    return selected, "unique"


def evaluate_open_flags(s: int, graph: str, open_flags: Sequence[bool]) -> LinearDecision:
    """Evaluate the frozen block event in O(s^2) time and memory."""

    count = site_count(s)
    if graph not in GRAPHS:
        raise ValueError("unknown graph: %s" % graph)
    if len(open_flags) != count:
        raise ValueError("expected %d site flags, got %d" % (count, len(open_flags)))
    flags = tuple(bool(value) for value in open_flags)
    width = 2 * s
    full = DisjointSet(flags)
    halves = DisjointSet(flags)
    full_edge_checks = 0
    half_edge_checks = 0

    for y in range(s):
        for x in range(width):
            index = y * width + x
            for other_x, other_y in _previous_neighbors(x, y, width, graph):
                other = other_y * width + other_x
                full_edge_checks += 1
                same_half = (x < s) == (other_x < s)
                if same_half:
                    half_edge_checks += 1
                if flags[index] and flags[other]:
                    full.union(index, other)
                    if same_half:
                        halves.union(index, other)

    left_list = [index for index in range(count) if index % width < s]
    right_list = [index for index in range(count) if index % width >= s]
    left_selected, left_status = _unique_largest(halves, flags, left_list)
    if not left_selected:
        return LinearDecision(False, "left_%s" % left_status, (), (), full_edge_checks, half_edge_checks)
    right_selected, right_status = _unique_largest(halves, flags, right_list)
    if not right_selected:
        return LinearDecision(
            False,
            "right_%s" % right_status,
            left_selected,
            (),
            full_edge_checks,
            half_edge_checks,
        )
    connected = full.find(left_selected[0]) == full.find(right_selected[0])
    return LinearDecision(
        connected,
        "success" if connected else "largest_clusters_disconnected",
        left_selected,
        right_selected,
        full_edge_checks,
        half_edge_checks,
    )


def _exact_indexes(vertices: Sequence[Tuple[int, int]], selected: Sequence[Tuple[int, int]]) -> Tuple[int, ...]:
    lookup = {vertex: index for index, vertex in enumerate(vertices)}
    return tuple(sorted(lookup[vertex] for vertex in selected))


def build_artifact() -> dict[str, Any]:
    comparisons = 0
    mismatches = []
    for graph in GRAPHS:
        for s in (1, 2):
            vertices = exact.ordered_vertices(s)
            for mask in range(1 << len(vertices)):
                expected = exact.decide_event(mask, s, graph)
                observed = evaluate_open_flags(s, graph, flags_from_mask(mask, len(vertices)))
                expected_left = _exact_indexes(vertices, expected.left_largest)
                expected_right = _exact_indexes(vertices, expected.right_largest)
                comparisons += 1
                if (
                    expected.success != observed.success
                    or expected.reason != observed.reason
                    or expected_left != observed.left_largest
                    or expected_right != observed.right_largest
                ):
                    mismatches.append({"graph": graph, "s": s, "mask": mask})

    large_controls = {}
    for graph in GRAPHS:
        s = 64
        decision = evaluate_open_flags(s, graph, (True,) * site_count(s))
        large_controls[graph] = {
            "s": s,
            "success": decision.success,
            "left_largest_size": len(decision.left_largest),
            "right_largest_size": len(decision.right_largest),
            "full_edge_checks": decision.full_edge_checks,
            "half_edge_checks": decision.half_edge_checks,
        }

    witness_before = evaluate_open_flags(2, "square", flags_from_mask(6, 8))
    witness_after = evaluate_open_flags(2, "square", flags_from_mask(22, 8))
    assert comparisons == 520
    assert not mismatches
    assert witness_before.success
    assert witness_after.reason == "left_tie"
    assert all(control["success"] for control in large_controls.values())

    return {
        "schema": "matching-one/pc-block-event-linear/v1",
        "issue": 112,
        "status": "deterministic_evaluator_validated",
        "complexity": {"time": "O(s^2)", "memory": "O(s^2)"},
        "differential_oracle": {
            "sizes": [1, 2],
            "graphs": list(GRAPHS),
            "configuration_graph_pairs": comparisons,
            "mismatch_count": len(mismatches),
        },
        "large_all_open_controls": large_controls,
        "nonmonotonicity_witness": {
            "before_mask": 6,
            "before_reason": witness_before.reason,
            "after_mask": 22,
            "after_reason": witness_after.reason,
        },
        "claim_boundary": {
            "proved": "deterministic evaluator agrees with the exact oracle on all s=1,2 cases",
            "not_proved": (
                "random-sampler independence, production event probabilities, timing guarantees, "
                "or a new critical-probability bound"
            ),
        },
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    oracle = artifact["differential_oracle"]
    lines = [
        "# Linear open-boundary block-event evaluator",
        "",
        "The deterministic evaluator uses two disjoint-set forests: one for connectivity in the",
        "full rectangle and one with cross-half edges suppressed for largest-cluster selection.",
        "It has `O(s^2)` time and memory complexity.",
        "",
        "## Differential validation",
        "",
        "- configuration/graph pairs compared: `%d`;" % oracle["configuration_graph_pairs"],
        "- exact-oracle mismatches: `%d`;" % oracle["mismatch_count"],
        "- covered sizes: `s=1,2`; graphs: square NN and matching NN+diagonals.",
        "",
        "## Deterministic large controls",
        "",
        "| graph | s | selected sites per half | full edge checks | half edge checks | result |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for graph in GRAPHS:
        row = artifact["large_all_open_controls"][graph]
        lines.append(
            "| %s | %d | %d | %d | %d | %s |"
            % (
                graph,
                row["s"],
                row["left_largest_size"],
                row["full_edge_checks"],
                row["half_edge_checks"],
                "success" if row["success"] else "failure",
            )
        )
    lines.extend(
        [
            "",
            "The frozen square-graph nonmonotonicity witness is also preserved: mask `6` succeeds,",
            "while opening one additional site to form mask `22` fails with `left_tie`.",
            "",
            "## Boundary",
            "",
            "No stochastic sampler is included. This result does not establish independent trials,",
            "an event probability near `p_c`, a runtime guarantee, or a new certified threshold bound.",
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
