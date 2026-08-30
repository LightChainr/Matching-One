#!/usr/bin/env python3
"""Second exact Hall compression for the corrected P334 reservoir.

After quotienting simultaneous translation, the relative translation phase of
the flat face remains as an N-fold family of source twins.  The corrected
reservoir translates that face over the entire finite HNF group, so all N
twins have exactly the same target neighbourhood.  This scorer replaces them
by one demand-N node and solves the resulting capacitated Hall problem.

The frozen corrected-reservoir observable is not changed.  Combined, MM-only
and YN-only flows are scored from the same target sets, and every failed flow
returns a deterministic minimum-cut certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

from p334_tm_configuration_cross_switch import translation_permutations
from p334_tm_corrected_reservoir_scan import rows_for_order
from p334_tm_translation_orbit_hall import (
    descriptor,
    hard_rows,
    inverse_to_origin,
    normalize_source,
    normalize_target,
    translate_face,
    transverse_reservoir_targets,
)


SCHEMA = "p334-tm-coarse-reservoir-hall-v1"


class _Edge:
    __slots__ = ("to", "reverse", "capacity")

    def __init__(self, to: int, reverse: int, capacity: int):
        self.to = to
        self.reverse = reverse
        self.capacity = capacity


class _Dinic:
    def __init__(self, nodes: int):
        self.graph: list[list[_Edge]] = [[] for _ in range(nodes)]

    def add_edge(self, source: int, target: int, capacity: int) -> None:
        forward = _Edge(target, len(self.graph[target]), capacity)
        reverse = _Edge(source, len(self.graph[source]), 0)
        self.graph[source].append(forward)
        self.graph[target].append(reverse)

    def maximum_flow(self, source: int, sink: int) -> int:
        total = 0
        nodes = len(self.graph)
        while True:
            level = [-1] * nodes
            level[source] = 0
            queue = deque([source])
            while queue:
                node = queue.popleft()
                for edge in self.graph[node]:
                    if edge.capacity and level[edge.to] < 0:
                        level[edge.to] = level[node] + 1
                        queue.append(edge.to)
            if level[sink] < 0:
                return total
            cursors = [0] * nodes

            def send(node: int, available: int) -> int:
                if node == sink:
                    return available
                while cursors[node] < len(self.graph[node]):
                    edge = self.graph[node][cursors[node]]
                    if edge.capacity and level[edge.to] == level[node] + 1:
                        pushed = send(edge.to, min(available, edge.capacity))
                        if pushed:
                            edge.capacity -= pushed
                            self.graph[edge.to][edge.reverse].capacity += pushed
                            return pushed
                    cursors[node] += 1
                return 0

            while True:
                pushed = send(source, 10**18)
                if not pushed:
                    break
                total += pushed

    def residual_reachable(self, source: int) -> set[int]:
        reached = {source}
        queue = deque([source])
        while queue:
            node = queue.popleft()
            for edge in self.graph[node]:
                if edge.capacity and edge.to not in reached:
                    reached.add(edge.to)
                    queue.append(edge.to)
        return reached


def coarse_source_key(source, inverses):
    replica, coexit, flat = source
    return (
        replica,
        translate_face(coexit, inverses[coexit[1]]),
        translate_face(flat, inverses[flat[1]]),
    )


def _class_payload(key) -> dict[str, Any]:
    replica, coexit, flat = key
    return {
        "replica": replica,
        "coexit_orbit_representative": list(coexit),
        "flat_orbit_representative": list(flat),
    }


def _digest(value: Any) -> str:
    rendered = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(rendered.encode("utf-8")).hexdigest()


def _compress_indices(indices: Sequence[int]) -> list[str]:
    if not indices:
        return []
    output = []
    start = previous = indices[0]
    for value in indices[1:]:
        if value == previous + 1:
            previous = value
            continue
        output.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = value
    output.append(str(start) if start == previous else f"{start}-{previous}")
    return output


def _filter_targets(targets: Iterable[tuple], channel: str) -> frozenset[tuple]:
    if channel == "combined":
        return frozenset(targets)
    if channel not in ("MM", "YN"):
        raise ValueError(f"unknown target channel: {channel}")
    return frozenset(target for target in targets if target[0] == channel)


def capacitated_hall(
    class_keys: Sequence[tuple],
    neighborhoods: Sequence[frozenset[tuple]],
    demand: int,
    channel: str,
) -> dict[str, Any]:
    filtered = [_filter_targets(row, channel) for row in neighborhoods]
    targets = sorted({target for row in filtered for target in row})
    target_index = {target: index for index, target in enumerate(targets)}
    classes = len(class_keys)
    total_demand = demand * classes
    source = 0
    class_start = 1
    target_start = class_start + classes
    sink = target_start + len(targets)
    network = _Dinic(sink + 1)
    infinite = total_demand + 1
    for class_index, row in enumerate(filtered):
        class_node = class_start + class_index
        network.add_edge(source, class_node, demand)
        for target in row:
            network.add_edge(
                class_node,
                target_start + target_index[target],
                infinite,
            )
    for index in range(len(targets)):
        network.add_edge(target_start + index, sink, 1)

    flow = network.maximum_flow(source, sink)
    result: dict[str, Any] = {
        "channel": channel,
        "coarse_classes": classes,
        "demand_per_class": demand,
        "total_demand": total_demand,
        "reachable_targets": len(targets),
        "maximum_flow": flow,
        "Hall_deficiency": total_demand - flow,
        "saturates": flow == total_demand,
        "minimum_class_degree": min(map(len, filtered), default=0),
        "maximum_class_degree": max(map(len, filtered), default=0),
    }
    if flow == total_demand:
        result["minimum_cut_certificate"] = {
            "class_count": 0,
            "class_index_ranges": [],
            "demand": 0,
            "neighbor_target_count": 0,
            "Hall_deficiency": 0,
            "status": "empty_cut_only_at_zero_deficiency",
        }
        return result

    reachable = network.residual_reachable(source)
    cut_indices = [
        index
        for index in range(classes)
        if class_start + index in reachable
    ]
    cut_neighbors = sorted(
        {
            target
            for index in cut_indices
            for target in filtered[index]
        }
    )
    cut_deficiency = demand * len(cut_indices) - len(cut_neighbors)
    if cut_deficiency != total_demand - flow:
        raise AssertionError("residual minimum cut does not reproduce the Hall deficiency")
    cut_payloads = [_class_payload(class_keys[index]) for index in cut_indices]
    result["minimum_cut_certificate"] = {
        "class_count": len(cut_indices),
        "class_index_ranges": _compress_indices(cut_indices),
        "class_descriptor_sha256": _digest(cut_payloads),
        "neighbor_target_sha256": _digest(cut_neighbors),
        "replica_histogram": dict(
            sorted(Counter(payload["replica"] for payload in cut_payloads).items())
        ),
        "demand": demand * len(cut_indices),
        "neighbor_target_count": len(cut_neighbors),
        "Hall_deficiency": cut_deficiency,
        "status": "exact_residual_minimum_cut",
    }
    return result


def coarse_row_audit(row, *, verify_all_twins: bool) -> dict[str, Any]:
    n, matrix, geometry, carrier, marks, line, lower_layer, faces = row
    permutations = translation_permutations(geometry)
    inverses = inverse_to_origin(permutations, n)
    raw_sources = [
        (replica, coexit, flat)
        for replica in range(4)
        for coexit in faces["D"]
        for flat in faces["F"]
    ]
    orbit_sources = sorted(
        {normalize_source(source, inverses) for source in raw_sources}
    )
    twins: dict[tuple, list[tuple]] = defaultdict(list)
    for source in orbit_sources:
        twins[coarse_source_key(source, inverses)].append(source)
    class_keys = sorted(twins)
    twin_sizes = [len(twins[key]) for key in class_keys]
    if any(size != n for size in twin_sizes):
        raise AssertionError("relative-phase twin classes must all have size N")

    neighborhoods = []
    verified_twins = 0
    for key in class_keys:
        members = sorted(twins[key])

        def neighbors(source):
            return frozenset(
                normalize_target(target, inverses)
                for target in transverse_reservoir_targets(
                    marks,
                    line,
                    source,
                    permutations,
                    n,
                    transport=True,
                )
            )

        reference = neighbors(members[0])
        if verify_all_twins:
            for member in members[1:]:
                if neighbors(member) != reference:
                    raise AssertionError("relative phase changed the corrected reservoir")
                verified_twins += 1
        neighborhoods.append(reference)

    channel_rows = {
        channel: capacitated_hall(class_keys, neighborhoods, n, channel)
        for channel in ("combined", "MM", "YN")
    }
    if not channel_rows["combined"]["saturates"]:
        pure_channel = "none_combined_failure"
    elif channel_rows["MM"]["saturates"] and channel_rows["YN"]["saturates"]:
        pure_channel = "both"
    elif channel_rows["MM"]["saturates"]:
        pure_channel = "MM"
    elif channel_rows["YN"]["saturates"]:
        pure_channel = "YN"
    else:
        pure_channel = "genuinely_mixed"

    return {
        **descriptor(n, matrix, carrier, line, lower_layer, faces),
        "source_compression": {
            "raw_sources": len(raw_sources),
            "translation_orbit_sources": len(orbit_sources),
            "coarse_twin_classes": len(class_keys),
            "twin_class_size": n,
            "raw_to_coarse_factor": n * n,
            "class_order_sha256": _digest([_class_payload(key) for key in class_keys]),
            "all_twin_sizes_exact": True,
            "neighborhood_equality": (
                "exhaustively_verified"
                if verify_all_twins
                else "exact_by_translation_reindexing"
            ),
            "extra_twins_compared": verified_twins,
        },
        "channel_flows": channel_rows,
        "pure_channel_classification": pure_channel,
    }


def known_n8_smith_gate():
    return next(
        row
        for row in hard_rows(8)
        if row[0] == 8
        and row[1] == ((2, 0), (0, 4))
        and row[3] == "matching"
        and row[5] == (1, 0)
        and row[6] == 4
    )


def build_result() -> dict[str, Any]:
    n6_rows = [
        coarse_row_audit(row, verify_all_twins=True)
        for row in rows_for_order(6)
    ]
    n8_row = coarse_row_audit(
        known_n8_smith_gate(),
        verify_all_twins=False,
    )
    if len(n6_rows) != 4:
        raise AssertionError("the frozen N6 gate must contain four rows")
    if not all(row["channel_flows"]["combined"]["saturates"] for row in n6_rows):
        raise AssertionError("a corrected N6 row lost combined saturation")
    if not all(row["pure_channel_classification"] == "YN" for row in n6_rows):
        raise AssertionError("the N6 rows must retain the YN-only injection")
    if not n8_row["channel_flows"]["combined"]["saturates"]:
        raise AssertionError("the corrected N8 Smith gate lost combined saturation")
    if n8_row["pure_channel_classification"] != "MM":
        raise AssertionError("the N8 Smith gate must retain the MM-only injection")

    return {
        "schema": SCHEMA,
        "parent_commit": "4bb75176c56558084c8397917995026e54420b9f",
        "observable": "unchanged one-carrier occupied-to-vacant exchange plus one transverse output-mark release",
        "theorem": {
            "source_bijection": "S_orbit = replicas x (D/Q) x (F/Q) x Q_relative",
            "twin_statement": "the corrected target neighborhood is independent of Q_relative because the flat-face translation loop reindexes g*h over the whole group",
            "capacitated_Hall": "saturation iff |union_(c in A) R(c)| >= N|A| for every coarse class set A",
            "deficiency": "max_A [N|A|-|union R(A)|]",
            "lifting": "integral coarse flow expands to relative-phase orbit matching, then the existing translation-orbit theorem lifts it to the raw injection",
            "Smith_boundary": "only finiteness and regular translation are used; Q need not be cyclic",
        },
        "N6_rows": n6_rows,
        "N8_Smith_2_4_gate": n8_row,
        "summary": {
            "N6_rows": 4,
            "N6_combined_saturated": 4,
            "N6_YN_only_saturated": 4,
            "N6_MM_only_saturated": 0,
            "N8_combined_saturated": True,
            "N8_MM_only_saturated": True,
            "N8_YN_only_saturated": False,
            "decision": "the exact second compression is general, while the successful pure target channel changes between the minimal and Smith gates",
        },
        "scientific_boundary": (
            "This proves the twin-class reduction and returns exact Hall cuts. "
            "It does not prove corrected-reservoir saturation for every HNF; "
            "that requires a uniform bound on target-neighborhood overlap or a constructive channel rule."
        ),
    }


def render_markdown(result: dict[str, Any]) -> str:
    n8 = result["N8_Smith_2_4_gate"]
    lines = [
        "# P334 coarse twin-class Hall certificate",
        "",
        "The corrected reservoir is unchanged.  Relative flat-face phase produces exactly `N` source twins with one common neighborhood, so each twin class becomes one demand-`N` flow node.",
        "",
        "| gate | raw | orbit | coarse | combined | MM | YN | pure channel |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for index, row in enumerate(result["N6_rows"]):
        compression = row["source_compression"]
        flows = row["channel_flows"]
        lines.append(
            f"| N6 row {index} | {compression['raw_sources']} | "
            f"{compression['translation_orbit_sources']} | "
            f"{compression['coarse_twin_classes']} | "
            f"{flows['combined']['maximum_flow']}/{flows['combined']['total_demand']} | "
            f"{flows['MM']['maximum_flow']}/{flows['MM']['total_demand']} | "
            f"{flows['YN']['maximum_flow']}/{flows['YN']['total_demand']} | "
            f"{row['pure_channel_classification']} |"
        )
    compression = n8["source_compression"]
    flows = n8["channel_flows"]
    lines.append(
        f"| N8 Smith-(2,4) | {compression['raw_sources']} | "
        f"{compression['translation_orbit_sources']} | "
        f"{compression['coarse_twin_classes']} | "
        f"{flows['combined']['maximum_flow']}/{flows['combined']['total_demand']} | "
        f"{flows['MM']['maximum_flow']}/{flows['MM']['total_demand']} | "
        f"{flows['YN']['maximum_flow']}/{flows['YN']['total_demand']} | "
        f"{n8['pure_channel_classification']} |"
    )
    lines.extend(
        [
            "",
            "Every failed channel stores a residual minimum-cut certificate with compressed class indices, exact demand, target-neighborhood size, deficiency and deterministic hashes.",
            "",
            result["scientific_boundary"],
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = build_result()
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    args.markdown.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()
