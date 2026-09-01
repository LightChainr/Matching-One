#!/usr/bin/env python3
"""Produce the fixed N25 distance-two #537 typed-carrier witness."""

from __future__ import annotations

import argparse
from collections import deque
import csv
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
KERNEL = ROOT / "experiments" / "p537-landing-matrix-preflight-20260901" / "kernel.tsv"
N = 25
ROW_X, ROW_Y, ROW_Z = 0, 6, 2
ROW_OCCUPIED_OFF_Z = (1, 3, 4, 5, 7, 9, 10, 12, 15, 16, 17)


class DSU:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.weight = [1] * size

    def root(self, vertex: int) -> int:
        while self.parent[vertex] != vertex:
            self.parent[vertex] = self.parent[self.parent[vertex]]
            vertex = self.parent[vertex]
        return vertex

    def join(self, left: int, right: int) -> None:
        left, right = self.root(left), self.root(right)
        if left == right:
            return
        if self.weight[left] < self.weight[right]:
            left, right = right, left
        self.parent[right] = left
        self.weight[left] += self.weight[right]


def geometry() -> tuple[list[tuple[int, int]], list[list[int]], list[list[int]], list[list[int]], list[int]]:
    """Reconstruct the axis (5,0) quotient and the 5x5 row-major map."""

    a, b = 5, 0

    def quotient_key(x: int, y: int) -> int:
        return N * ((a * x + b * y) % N) + ((-b * x + a * y) % N)

    index = [-1] * (N * N)
    representatives = [(0, 0)]
    index[quotient_key(0, 0)] = 0
    cursor = 0
    while cursor < len(representatives):
        x, y = representatives[cursor]
        cursor += 1
        for dx, dy in ((1, 0), (0, 1)):
            key = quotient_key(x + dx, y + dy)
            if index[key] < 0:
                index[key] = len(representatives)
                representatives.append((x + dx, y + dy))
    if len(representatives) != N:
        raise AssertionError("axis quotient does not have 25 vertices")

    neighbors: list[list[int]] = []
    diagonals: list[list[int]] = []
    for x, y in representatives:
        neighbors.append(
            [
                index[quotient_key(x, y + 1)],
                index[quotient_key(x + 1, y)],
                index[quotient_key(x, y - 1)],
                index[quotient_key(x - 1, y)],
            ]
        )
        diagonals.append(
            [
                index[quotient_key(x + 1, y + 1)],
                index[quotient_key(x - 1, y + 1)],
                index[quotient_key(x - 1, y - 1)],
                index[quotient_key(x + 1, y - 1)],
            ]
        )

    edge_id = [[-1] * 4 for _ in range(N)]
    next_edge = 0
    for vertex in range(N):
        for direction in (0, 1):
            other = neighbors[vertex][direction]
            edge_id[vertex][direction] = edge_id[other][direction + 2] = next_edge
            next_edge += 1
    if next_edge != 2 * N:
        raise AssertionError("wrong physical edge count")
    row_major = [index[quotient_key(cell % 5, cell // 5)] for cell in range(N)]
    return representatives, neighbors, diagonals, edge_id, row_major


def evaluate(
    occupied: list[bool], neighbors: list[list[int]], diagonals: list[list[int]]
) -> tuple[int, list[int], list[int], list[int]]:
    black, white = DSU(N), DSU(N)
    k = sum(occupied)
    edges = 0
    faces = 0
    for vertex in range(N):
        if occupied[vertex]:
            for direction in (0, 1):
                other = neighbors[vertex][direction]
                if occupied[other]:
                    edges += 1
                    black.join(vertex, other)
        else:
            for direction in (0, 1):
                other = neighbors[vertex][direction]
                if not occupied[other]:
                    white.join(vertex, other)
                other = diagonals[vertex][direction]
                if not occupied[other]:
                    white.join(vertex, other)
        if (
            occupied[vertex]
            and occupied[neighbors[vertex][1]]
            and occupied[neighbors[vertex][0]]
            and occupied[diagonals[vertex][0]]
        ):
            faces += 1
    black_root = [black.root(v) if occupied[v] else -1 for v in range(N)]
    white_root = [white.root(v) if not occupied[v] else -1 for v in range(N)]
    black_components = len(set(black_root) - {-1})
    white_components = len(set(white_root) - {-1})
    q = black_components - white_components - (k - edges + faces)
    degree = [
        sum(occupied[neighbors[v][d]] for d in range(4)) if occupied[v] else 0
        for v in range(N)
    ]
    return q, black_root, white_root, degree


def carrier_key(
    occupied: list[bool],
    black_root: list[int],
    neighbors: list[list[int]],
    edge_id: list[list[int]],
    centers: tuple[int, ...],
    bits: int,
) -> tuple[int, list[int]]:
    canonical: dict[int, int] = {}
    labels: list[int] = []
    key = 0
    for center in centers:
        for direction in range(4):
            vertex = neighbors[center][direction]
            identity = black_root[vertex] if occupied[vertex] else N + edge_id[center][direction]
            if identity not in canonical:
                canonical[identity] = len(canonical)
            label = canonical[identity]
            labels.append(label)
            key |= label << (bits * (len(labels) - 1))
    return key, labels


def partition(roots: list[int], occupied: list[bool], ports: list[int], black: bool) -> list[int]:
    canonical: dict[int, int] = {}
    result: list[int] = []
    for port in ports:
        if occupied[port] != black:
            result.append(-1)
            continue
        root = roots[port]
        if root not in canonical:
            canonical[root] = len(canonical)
        result.append(canonical[root])
    return result


def shortest_distance(left: int, right: int, neighbors: list[list[int]]) -> int:
    queue = deque([(left, 0)])
    seen = {left}
    while queue:
        vertex, distance = queue.popleft()
        if vertex == right:
            return distance
        for other in neighbors[vertex]:
            if other not in seen:
                seen.add(other)
                queue.append((other, distance + 1))
    raise AssertionError("connected quotient reported no path")


def read_kernel() -> dict[int, int]:
    with KERNEL.open(newline="") as handle:
        rows = (row for row in handle if row.strip() and not row.startswith("#"))
        reader = csv.DictReader(rows, delimiter="\t")
        key_field = "key" if "key" in (reader.fieldnames or ()) else "packed_key"
        return {int(row[key_field]): int(row["g16"]) for row in reader}


def produce() -> dict[str, object]:
    representatives, neighbors, diagonals, edge_id, row_major = geometry()
    x, y, z = (row_major[cell] for cell in (ROW_X, ROW_Y, ROW_Z))
    if (x, y, z) != (0, 4, 3):
        raise AssertionError("row-major to quotient mapping drifted")
    occupied_vertices = sorted(row_major[cell] for cell in ROW_OCCUPIED_OFF_Z)
    occupied = [vertex in occupied_vertices for vertex in range(N)]
    if occupied[x] or occupied[y] or occupied[z]:
        raise AssertionError("x, y, and z must be off in the frozen background")

    z_cardinal = neighbors[z]
    z_corner = [diagonals[z][0], diagonals[z][3], diagonals[z][2], diagonals[z][1]]
    arm_mask = sum(int(occupied[z_cardinal[d]]) << d for d in range(4))
    corner_mask = sum(int(occupied[z_corner[d]]) << d for d in range(4))
    distances = {
        "x_y": shortest_distance(x, y, neighbors),
        "x_z": shortest_distance(x, z, neighbors),
        "y_z": shortest_distance(y, z, neighbors),
    }

    kernel = read_kernel()
    nodes: list[dict[str, int]] = []
    joint_labels: list[list[int]] = []
    black_partitions: list[list[int]] = []
    white_cut: list[int] | None = None
    source_masks: list[int] = []
    global_contacts: list[bool] = []
    degree_branches: list[bool] = []
    for state in (0, 1):
        occupied[z] = bool(state)
        q, black_root, white_root, degree = evaluate(occupied, neighbors, diagonals)
        bell, _ = carrier_key(occupied, black_root, neighbors, edge_id, (x, y), 3)
        joint, labels = carrier_key(occupied, black_root, neighbors, edge_id, (x, y, z), 4)
        g16 = kernel.get(bell, 0)
        nodes.append(
            {
                "state": state,
                "q": q,
                "rank_index_q_plus_1": q + 1,
                "E": q * q,
                "bell": bell,
                "joint_C": joint,
                "g16": g16,
            }
        )
        joint_labels.append(labels)
        black_partitions.append(partition(black_root, occupied, z_cardinal, True))
        if state == 1:
            white_cut = partition(white_root, occupied, z_cardinal, False)
        source_masks.append(
            sum(
                int(occupied[neighbors[center][direction]]) << (4 * side + direction)
                for side, center in enumerate((x, y))
                for direction in range(4)
            )
        )
        landing_roots = {
            black_root[port] for port in z_cardinal if occupied[port]
        }
        global_contacts.append(
            any(
                occupied[neighbors[center][direction]]
                and black_root[neighbors[center][direction]] in landing_roots
                for center in (x, y)
                for direction in range(4)
            )
        )
        z_ports = set(z_cardinal)
        degree_branches.append(
            any(
                occupied[vertex]
                and vertex not in z_ports
                and black_root[vertex] in landing_roots
                and degree[vertex] >= 3
                for vertex in range(N)
            )
        )
    occupied[z] = False

    expected_nodes = [
        {"state": 0, "q": -1, "rank_index_q_plus_1": 0, "E": 1, "bell": 274568, "joint_C": 21990249529872, "g16": 8},
        {"state": 1, "q": 0, "rank_index_q_plus_1": 1, "E": 0, "bell": 8256, "joint_C": 3298535014656, "g16": 0},
    ]
    if nodes != expected_nodes or distances != {"x_y": 2, "x_z": 2, "y_z": 2}:
        raise AssertionError({"nodes": nodes, "distances": distances})
    terminal_incidence = [len(set(labels[:8]) & set(labels[8:])) for labels in joint_labels]
    if terminal_incidence != [2, 1] or not global_contacts[0] or degree_branches[0]:
        raise AssertionError("typed-carrier contact contract drifted")

    before, after = nodes
    stats = {
        "count": 1,
        "sum_q0": before["q"],
        "sum_E0": before["E"],
        "sum_a16_0": before["g16"],
        "sum_q0_a16_0": before["q"] * before["g16"],
        "sum_E0_a16_0": before["E"] * before["g16"],
        "sum_q1": after["q"],
        "sum_E1": after["E"],
        "sum_a16_1": after["g16"],
        "sum_q1_a16_1": after["q"] * after["g16"],
        "sum_E1_a16_1": after["E"] * after["g16"],
    }
    return {
        "schema": "matching-one/p537-one-defect-witness/v1",
        "status": "fixed_nonadjacent_typed_carrier_diagonal_edge",
        "N": N,
        "geometry": {"id": "axis", "a": 5, "b": 0},
        "fixed_contract": {
            "sampling": "none",
            "search": "none; reconstruct the frozen row-major configuration",
            "row_major_order": "cell=5*row+column on the 5x5 axis fundamental square",
        },
        "transition_id": "axis-N25:row-x0:y6:z2:fixed-nonadjacent:0to1",
        "k_minus": len(occupied_vertices),
        "vertices": {
            "x": x,
            "y": y,
            "z": z,
            "x_coordinate": list(representatives[x]),
            "y_coordinate": list(representatives[y]),
            "z_coordinate": list(representatives[z]),
        },
        "row_major_configuration": {
            "x": ROW_X,
            "y": ROW_Y,
            "z": ROW_Z,
            "occupied_off_z": list(ROW_OCCUPIED_OFF_Z),
            "row_to_internal_vertex": row_major,
        },
        "occupied_background_vertices": occupied_vertices,
        "pairwise_NN_distances": {**distances, "minimum": min(distances.values())},
        "collar": {
            "arm_order": ["N", "E", "S", "W"],
            "arm_mask": arm_mask,
            "alternating_four_arm": arm_mask in (5, 10),
            "outer_join_identity_applicable": arm_mask in (5, 10),
            "corner_mask_NE_SE_SW_NW": corner_mask,
            "outer_black_join_J_B": 0,
            "outer_white_join_J_W": 0,
        },
        "outer_C": {
            "meaning": "global first-occurrence carrier partition on x4+y4+z4 ports",
            "before_key": before["joint_C"],
            "before_labels": joint_labels[0],
            "after_key": after["joint_C"],
            "after_labels": joint_labels[1],
            "terminal_incidence_before": terminal_incidence[0],
            "terminal_incidence_after": terminal_incidence[1],
        },
        "outer_B": {
            "meaning": "global NN-black partition of the four z-cardinal ports; -1 means vacant",
            "before": black_partitions[0],
            "after": black_partitions[1],
        },
        "outer_W": {
            "meaning": "global matching-white partition in the off-z cut; -1 means occupied",
            "cut_partition": white_cut,
        },
        "source_ports": {
            "occupied_mask_before": source_masks[0],
            "occupied_mask_after": source_masks[1],
            "source_component": "diag1",
            "bell_before": before["bell"],
            "bell_after": after["bell"],
            "bell_from_joint_before": before["bell"],
            "bell_from_joint_after": after["bell"],
            "g16_before": before["g16"],
            "g16_after": after["g16"],
        },
        "off_port": {
            "measured_in_state": 0,
            "source_touches_global_black_landing_component": global_contacts[0],
            "black_landing_component_has_degree3_vertex_outside_z_ports": degree_branches[0],
        },
        "carrier_scope": {
            "classification": "joint_incidence_typed_carrier",
            "distance_at_most_one": False,
            "joint_terminal_incidence_before_after": terminal_incidence,
            "annular_separation_certified": False,
            "reason": "all marked centers are NN-distance two, but the common carrier loses one source/thermal terminal incidence",
        },
        "nodes": nodes,
        "edge": {
            "physical_move": "occupy z with every other site fixed",
            "changes_landing_rank": True,
            "changes_source_Bell": True,
            "changes_source_value_g16": True,
            "changes_joint_component_map": True,
        },
        "positive_signed_statistics": {
            "positive_count_each_state": 1,
            "signed_fields": ["q", "E=q^2", "a16=g16", "q*a16", "E*a16"],
            "normalization": "scorer converts a16 to a=g16/(16*N)",
        },
        "sufficient_statistics": stats,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(produce(), indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
