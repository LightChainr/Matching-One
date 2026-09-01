#!/usr/bin/env python3
"""Independent topology and certified-sign checks for the #537 witness."""

from __future__ import annotations

import csv
import json
from fractions import Fraction as F
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
WITNESS = ROOT / "results" / "p537-one-defect-gate-20260901" / "witness.json"
RESULT = ROOT / "results" / "p537-one-defect-gate-20260901" / "result.json"
KERNEL = ROOT / "experiments" / "p537-landing-matrix-preflight-20260901" / "kernel.tsv"
N = 25


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


def geometry() -> tuple[list[list[int]], list[list[int]], list[list[int]]]:
    a, b = 5, 0

    def mod(value: int) -> int:
        return value % N

    def quotient_key(x: int, y: int) -> int:
        return N * mod(a * x + b * y) + mod(-b * x + a * y)

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
        raise AssertionError("wrong quotient size")
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
            reverse = direction + 2
            if edge_id[vertex][direction] >= 0 or edge_id[other][reverse] >= 0:
                raise AssertionError("bad reciprocal edge")
            edge_id[vertex][direction] = edge_id[other][reverse] = next_edge
            next_edge += 1
    if next_edge != 2 * N:
        raise AssertionError("wrong physical edge count")
    return neighbors, diagonals, edge_id


def evaluate(occupied: list[bool], neighbors: list[list[int]], diagonals: list[list[int]]) -> tuple[int, list[int], list[int]]:
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
    occupied_root = [black.root(v) if occupied[v] else -1 for v in range(N)]
    vacant_root = [white.root(v) if not occupied[v] else -1 for v in range(N)]
    black_components = len({occupied_root[v] for v in range(N) if occupied[v]})
    white_components = len({vacant_root[v] for v in range(N) if not occupied[v]})
    q = black_components - white_components - (k - edges + faces)
    return q, occupied_root, vacant_root


def carrier_keys(
    occupied: list[bool],
    occupied_root: list[int],
    neighbors: list[list[int]],
    edge_id: list[list[int]],
    centers: tuple[int, ...],
    bits: int,
) -> int:
    canonical: dict[int, int] = {}
    result = 0
    for side, center in enumerate(centers):
        for direction in range(4):
            vertex = neighbors[center][direction]
            identity = occupied_root[vertex] if occupied[vertex] else N + edge_id[center][direction]
            if identity not in canonical:
                canonical[identity] = len(canonical)
            result |= canonical[identity] << (bits * (4 * side + direction))
    return result


def read_kernel() -> dict[int, int]:
    with KERNEL.open(newline="") as handle:
        rows = (row for row in handle if row.strip() and not row.startswith("#"))
        reader = csv.DictReader(rows, delimiter="\t")
        key_field = "key" if "key" in (reader.fieldnames or ()) else "packed_key"
        return {int(row[key_field]): int(row["g16"]) for row in reader}


class P537OneDefectGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.witness = json.loads(WITNESS.read_text())
        cls.result = json.loads(RESULT.read_text())

    def test_independent_python_topology_reconstructs_literal_edge(self) -> None:
        neighbors, diagonals, edge_id = geometry()
        free_sites = self.witness["free_site_order"]
        mask = int(self.witness["background_mask"])
        occupied = [False] * N
        for bit, vertex in enumerate(free_sites):
            occupied[vertex] = bool((mask >> bit) & 1)
        x = int(self.witness["vertices"]["x"])
        y = int(self.witness["vertices"]["y"])
        z = int(self.witness["vertices"]["z"])
        self.assertFalse(occupied[x])
        self.assertFalse(occupied[y])
        self.assertFalse(occupied[z])

        reconstructed = []
        kernel = read_kernel()
        for state in (0, 1):
            occupied[z] = bool(state)
            q, occupied_root, _ = evaluate(occupied, neighbors, diagonals)
            bell = carrier_keys(occupied, occupied_root, neighbors, edge_id, (x, y), 3)
            joint = carrier_keys(occupied, occupied_root, neighbors, edge_id, (x, y, z), 4)
            reconstructed.append((q, bell, joint, kernel.get(bell, 0)))
        self.assertEqual(reconstructed, [(-1, 9240712, 23090870354448, 4), (0, 6848576, 92359816642816, 0)])

    def test_machine_result_has_allocation_robust_certified_sign(self) -> None:
        self.assertEqual(
            self.result["status"],
            "allocation_robust_physical_diagonal_edge_nonzero",
        )
        edge = self.result["edge_weight_C4_orbit_pooled"]
        source = edge["source_midpoint_part"]
        counterterm = edge["root_counterterm_part"]
        full = edge["full"]
        self.assertLess(F(source["upper"]), 0)
        self.assertGreater(F(counterterm["lower"]), 0)
        self.assertLess(F(full["upper"]), 0)
        self.assertTrue(edge["allocation_robust"])
        self.assertEqual(self.result["topology"]["delta_a_exact"], "-1/100")

    def test_scope_is_contact_and_does_not_promote_separated_sector(self) -> None:
        topology = self.result["topology"]
        self.assertFalse(topology["collar"]["alternating_four_arm"])
        self.assertFalse(topology["collar"]["outer_join_identity_applicable"])
        self.assertTrue(topology["off_port"]["source_touches_global_black_landing_component"])
        decision = self.result["stop_decision"]
        self.assertEqual(decision["blanket_full_graph_two_independent_defect_route"], "falsified")
        self.assertIn("open", decision["separated_sector"])


if __name__ == "__main__":
    unittest.main()
