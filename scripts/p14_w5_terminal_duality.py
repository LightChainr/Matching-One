#!/usr/bin/env python3
"""Exact two-orbit W5 terminal distribution and planar-duality certificate."""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb, gcd
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "terminal-reliability" / "p14-w5-terminal-duality.json"

# Cyclic positions q_0,q_1,q_2,q_3 on the disk boundary.
Q = (0, 2, 1, 3)
HUB = 4
RIM = tuple(tuple(sorted((Q[i], Q[(i + 1) % 4]))) for i in range(4))
SPOKE = tuple(tuple(sorted((HUB, Q[i]))) for i in range(4))
EDGE_NAMES = tuple(f"R_{i}" for i in range(4)) + tuple(f"S_{i}" for i in range(4))


def canonical_partition(groups: Iterable[Iterable[int]]) -> tuple[tuple[int, ...], ...]:
    return tuple(sorted((tuple(sorted(group)) for group in groups), key=lambda group: group[0]))


def partition_key(partition: Sequence[Sequence[int]]) -> str:
    return "|".join("".join(str(value) for value in group) for group in partition)


def terminal_partition(
    vertices: Sequence[int], edges: Sequence[tuple[int, int]], terminals: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
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
    groups: dict[int, list[int]] = {}
    for terminal in terminals:
        groups.setdefault(find(terminal), []).append(terminal)
    return canonical_partition(groups.values())


def all_partitions(values: Sequence[int]) -> list[tuple[tuple[int, ...], ...]]:
    result: set[tuple[tuple[int, ...], ...]] = set()

    def recurse(index: int, groups: list[list[int]]) -> None:
        if index == len(values):
            result.add(canonical_partition(groups))
            return
        value = values[index]
        for group_index in range(len(groups)):
            groups[group_index].append(value)
            recurse(index + 1, groups)
            groups[group_index].pop()
        groups.append([value])
        recurse(index + 1, groups)
        groups.pop()

    recurse(0, [])
    return sorted(result, key=partition_key)


def primal_partition(bits: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    edges = [edge for edge, is_open in zip(RIM + SPOKE, bits) if is_open]
    return terminal_partition(tuple(range(5)), edges, Q)


def spherical_transform(bits: Sequence[int]) -> tuple[int, ...]:
    """Complement-dual under O->h and F_i->q_-i (an involution)."""

    rim = tuple(1 - bits[4 + ((-j) % 4)] for j in range(4))
    spoke = tuple(1 - bits[(-j) % 4] for j in range(4))
    return rim + spoke


def relative_dual_partition(bits: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    """Partition of the four split-outer-face boundary terminals B_i."""

    # Internal face vertices F_i are 0..3; boundary terminals B_i are 4..7.
    edges: list[tuple[int, int]] = []
    for i in range(4):
        if not bits[i]:
            edges.append((i, 4 + i))  # R_i*=(F_i,B_i)
        if not bits[4 + i]:
            edges.append(((i - 1) % 4, i))  # S_i*=(F_(i-1),F_i)
    raw = terminal_partition(tuple(range(8)), edges, tuple(range(4, 8)))
    return canonical_partition(tuple(value - 4 for value in group) for group in raw)


def open_counts(bits: Sequence[int]) -> tuple[int, int]:
    return sum(bits[:4]), sum(bits[4:])


def line_power_polynomial(counts: Sequence[Sequence[int]]) -> list[int]:
    """Substitute s=1-r into a bivariate degree-(4,4) Bernstein table."""

    coefficients = [0] * 9
    for a in range(5):
        for b in range(5):
            count = counts[a][b]
            # r^a (1-r)^(4-a) s^b (1-s)^(4-b), s=1-r.
            r_power = a + 4 - b
            one_minus_power = 4 - a + b
            for tail in range(one_minus_power + 1):
                coefficients[r_power + tail] += count * (-1) ** tail * comb(one_minus_power, tail)
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def subtract(left: Sequence[int], right: Sequence[int]) -> list[int]:
    length = max(len(left), len(right))
    result = [0] * length
    for index in range(length):
        result[index] = (left[index] if index < len(left) else 0) - (
            right[index] if index < len(right) else 0
        )
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    common = 0
    for value in result:
        common = gcd(common, abs(value))
    if common > 1:
        result = [value // common for value in result]
    return result


def evaluate(poly: Sequence[int], value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def degrees(vertices: Sequence[int], edges: Sequence[tuple[int, int]]) -> list[int]:
    result = {vertex: 0 for vertex in vertices}
    for left, right in edges:
        result[left] += 1
        result[right] += 1
    return [result[vertex] for vertex in vertices]


def build_certificate() -> dict[str, object]:
    partitions = all_partitions(tuple(range(4)))
    tables = {partition_key(partition): [[0] * 5 for _ in range(5)] for partition in partitions}
    rows: list[dict[str, object]] = []
    relative_outputs: dict[str, set[str]] = defaultdict(set)
    spherical_outputs: dict[str, set[str]] = defaultdict(set)

    for bits in product((0, 1), repeat=8):
        primal = primal_partition(bits)
        transformed = spherical_transform(bits)
        spherical = primal_partition(transformed)
        relative = relative_dual_partition(bits)
        a, b = open_counts(bits)
        primal_key = partition_key(primal)
        spherical_key = partition_key(spherical)
        relative_key = partition_key(relative)
        tables[primal_key][a][b] += 1
        relative_outputs[primal_key].add(relative_key)
        spherical_outputs[primal_key].add(spherical_key)
        rows.append(
            {
                "bits_R0_R1_R2_R3_S0_S1_S2_S3": "".join(map(str, bits)),
                "open_rim": a,
                "open_spoke": b,
                "primal_partition": primal_key,
                "spherical_dual_bits": "".join(map(str, transformed)),
                "spherical_dual_partition": spherical_key,
                "relative_dual_partition": relative_key,
            }
        )

    row_by_bits = {row["bits_R0_R1_R2_R3_S0_S1_S2_S3"]: row for row in rows}
    for row in rows:
        transformed = row_by_bits[row["spherical_dual_bits"]]
        if transformed["spherical_dual_bits"] != row["bits_R0_R1_R2_R3_S0_S1_S2_S3"]:
            raise AssertionError("spherical complement-duality is not involutive")
        a, b = int(row["open_rim"]), int(row["open_spoke"])
        if (int(transformed["open_rim"]), int(transformed["open_spoke"])) != (4 - b, 4 - a):
            raise AssertionError("edge-orbit count map failed")

    for a in range(5):
        for b in range(5):
            if sum(table[a][b] for table in tables.values()) != comb(4, a) * comb(4, b):
                raise AssertionError("bivariate Bernstein cell total failed")
    if any(len(outputs) != 1 for outputs in relative_outputs.values()):
        raise AssertionError("disk-relative partition map did not close")

    # Frozen witness order: minimum max(open edges), then bit strings.
    witnesses: list[tuple[int, str, str, dict[str, object], dict[str, object]]] = []
    for left, right in combinations(rows, 2):
        if left["primal_partition"] != right["primal_partition"]:
            continue
        if left["spherical_dual_partition"] == right["spherical_dual_partition"]:
            continue
        left_bits = str(left["bits_R0_R1_R2_R3_S0_S1_S2_S3"])
        right_bits = str(right["bits_R0_R1_R2_R3_S0_S1_S2_S3"])
        witnesses.append(
            (
                max(int(left["open_rim"]) + int(left["open_spoke"]), int(right["open_rim"]) + int(right["open_spoke"])),
                left_bits,
                right_bits,
                left,
                right,
            )
        )
    witness = min(witnesses, key=lambda item: item[:3])

    all_key, none_key = "0123", "0|1|2|3"
    all_line = line_power_polynomial(tables[all_key])
    none_line = line_power_polynomial(tables[none_key])
    defect_line = subtract(all_line, none_line)
    half_defect = evaluate(defect_line, Fraction(1, 2))

    canonical_rows = json.dumps(rows, sort_keys=True, separators=(",", ":"))
    relative_map = {key: next(iter(outputs)) for key, outputs in sorted(relative_outputs.items())}
    sphere_transition_support = {
        key: sorted(outputs) for key, outputs in sorted(spherical_outputs.items())
    }

    primal_edges = RIM + SPOKE
    relative_edges = tuple((i, 4 + i) for i in range(4)) + tuple(
        ((i - 1) % 4, i) for i in range(4)
    )
    return {
        "schema": "matching-one/p14-w5-terminal-duality/v1",
        "parent_commit": "2cc9466b49ed7b54fae85ef354e08f68163eda1f",
        "protocol_commit": "58a2bb3",
        "source_candidate": "4:5:0111111011",
        "edge_order": list(EDGE_NAMES),
        "conventions": {
            "cyclic_terminal_order": list(Q),
            "rim_edges": [list(edge) for edge in RIM],
            "spoke_edges": [list(edge) for edge in SPOKE],
            "rim_probability": "r",
            "spoke_probability": "s",
        },
        "primal_partition_distribution": {
            "basis": "count[a][b] * r^a(1-r)^(4-a)s^b(1-s)^(4-b)",
            "partition_count_including_zero": len(tables),
            "configuration_count": len(rows),
            "bivariate_bernstein_counts": tables,
            "identically_zero_partitions": [key for key, table in tables.items() if not any(map(any, table))],
        },
        "spherical_dual": {
            "vertices": ["O", "F_0", "F_1", "F_2", "F_3"],
            "edge_bijection": {
                "R_i_star": "(O,F_i), mapped by O->h,F_i->q_-i to S_-i",
                "S_i_star": "(F_(i-1),F_i), mapped to R_-i",
            },
            "parameter_map": "(r,s)->(1-s,1-r)",
            "fixed_line": "r+s=1",
            "configuration_map_is_involution": True,
            "fixed_line_weight_exponents": "(a+4-b,4-a+b), invariant under (a,b)->(4-b,4-a)",
            "terminal_partition_is_function_of_primal_partition": all(
                len(outputs) == 1 for outputs in spherical_outputs.values()
            ),
            "partition_transition_support": sphere_transition_support,
            "smallest_nonclosure_witness": {
                "max_open_edge_count": witness[0],
                "left": witness[3],
                "right": witness[4],
                "mechanism": "the extra open edge closes an internal triangle without changing the primal terminal partition, but changes the spherical-dual terminal partition",
            },
        },
        "disk_relative_dual": {
            "description": "internal F_i cycle plus one boundary leaf B_i at each F_i",
            "vertices": [f"F_{i}" for i in range(4)] + [f"B_{i}" for i in range(4)],
            "edge_bijection": {
                "R_i_star": "(F_i,B_i)",
                "S_i_star": "(F_(i-1),F_i)",
            },
            "primal_vertex_count": 5,
            "dual_vertex_count": 8,
            "primal_degree_multiset": sorted(degrees(tuple(range(5)), primal_edges)),
            "dual_degree_multiset": sorted(degrees(tuple(range(8)), relative_edges)),
            "primal_terminal_degrees": [3, 3, 3, 3],
            "dual_terminal_degrees": [1, 1, 1, 1],
            "boundary_terminal_preserving_isomorphic_to_primal": False,
            "relative_partition_is_function_of_primal_partition": True,
            "exact_planar_complement_map": relative_map,
        },
        "natural_line_scalar_diagnostic": {
            "line": "s=1-r",
            "P_all_minus_P_none_power_coefficients_low_to_high": defect_line,
            "value_at_r_equals_s_equals_one_half": str(half_defect),
            "identically_zero": defect_line == [0],
        },
        "configuration_rows_sha256": sha256(canonical_rows.encode("utf-8")).hexdigest(),
        "configurations": rows,
        "decision": "FAIL: W5 spherical self-duality exchanges edge orbits but is not the disk-relative boundary duality; the latter closes the planar partition complement but is an eight-vertex cycle-with-leaves cell, not W5.",
        "next_exact_object": "the joint primal/disk-relative-dual boundary connectivity (equivalently an alternating-boundary medial/Temperley-Lieb state) together with an explicit periodic tiling and comparison map",
        "claim_boundary": "exact finite-cell obstruction only; no threshold, critical-polynomial identity, periodic construction, stochastic comparison, or rigorous bound",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    certificate = build_certificate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
