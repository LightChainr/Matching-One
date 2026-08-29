#!/usr/bin/env python3
"""Exact tiny-graph oracles for a Q=1 matching derivative defect.

The primary object is the normalized Laurent FK/interface amplitude

    Xi(Q,p) = E_p[Q**d(A)],

where d(A) is a declared primal/complement cluster defect.  Xi(1,p)=1,
and its first log-Q derivative at Q=1 is exactly E_p[d(A)].
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb
from pathlib import Path
from typing import Callable, Iterable


Edge = tuple[int, int]


def component_count(vertex_count: int, edges: Iterable[Edge]) -> int:
    parent = list(range(vertex_count))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for a, b in edges:
        union(a, b)
    return len({find(v) for v in range(vertex_count)})


def induced_component_count(
    vertex_count: int, edges: Iterable[Edge], selected_mask: int
) -> int:
    selected = [v for v in range(vertex_count) if (selected_mask >> v) & 1]
    if not selected:
        return 0
    selected_set = set(selected)
    parent = {v: v for v in selected}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        if a in selected_set and b in selected_set:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra
    return len({find(v) for v in selected})


def bernstein_to_power(coefficients: list[int], degree: int) -> list[int]:
    """Expand sum_n coefficients[n] p^n (1-p)^(degree-n)."""
    power = [0] * (degree + 1)
    for n, coefficient in enumerate(coefficients):
        for j in range(degree - n + 1):
            power[n + j] += coefficient * ((-1) ** j) * comb(degree - n, j)
    return power


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def summarize_defects(defects: list[int], sizes: list[int], degree: int) -> dict:
    histograms: dict[int, Counter[int]] = defaultdict(Counter)
    tangent_bernstein = [0] * (degree + 1)
    for defect, size in zip(defects, sizes):
        histograms[size][defect] += 1
        tangent_bernstein[size] += defect

    denominator = 1 << degree
    mean = Fraction(sum(defects), denominator)
    second = Fraction(sum(d * d for d in defects), denominator)
    variance = second - mean * mean
    return {
        "configuration_count": len(defects),
        "defect_histogram_by_size": {
            str(n): {str(d): count for d, count in sorted(histograms[n].items())}
            for n in sorted(histograms)
        },
        "first_log_Q_tangent_bernstein_coefficients": tangent_bernstein,
        "first_log_Q_tangent_power_coefficients": bernstein_to_power(
            tangent_bernstein, degree
        ),
        "p_half": {
            "mean_defect": fraction_text(mean),
            "second_log_Q_tangent_variance": fraction_text(variance),
        },
    }


def edge_pair_oracle(
    primal_vertices: int,
    primal_edges: list[Edge],
    complement_vertices: int,
    complement_edges: list[Edge],
    local_term: Callable[[int], int],
) -> dict:
    degree = len(primal_edges)
    if len(complement_edges) != degree:
        raise ValueError("edge-complement pairing requires equal edge counts")
    full_mask = (1 << degree) - 1
    defects: list[int] = []
    sizes: list[int] = []
    for mask in range(1 << degree):
        occupied = [edge for i, edge in enumerate(primal_edges) if (mask >> i) & 1]
        complement = [
            edge for i, edge in enumerate(complement_edges) if not ((mask >> i) & 1)
        ]
        size = mask.bit_count()
        defect = (
            component_count(primal_vertices, occupied)
            - component_count(complement_vertices, complement)
            - local_term(size)
        )
        defects.append(defect)
        sizes.append(size)
    result = summarize_defects(defects, sizes, degree)
    result["configurationwise_zero"] = all(d == 0 for d in defects)
    result["complement_mask"] = full_mask
    return result


def site_pair_oracle(
    vertex_count: int, primal_edges: list[Edge], matching_edges: list[Edge]
) -> dict:
    full_mask = (1 << vertex_count) - 1
    defects: list[int] = []
    sizes: list[int] = []
    for mask in range(1 << vertex_count):
        defects.append(
            induced_component_count(vertex_count, primal_edges, mask)
            - induced_component_count(vertex_count, matching_edges, full_mask ^ mask)
        )
        sizes.append(mask.bit_count())
    return summarize_defects(defects, sizes, vertex_count)


def build_oracle() -> dict:
    triangle = [(0, 1), (1, 2), (2, 0)]
    triangle_dual = [(0, 1), (0, 1), (0, 1)]
    path4 = [(0, 1), (1, 2), (2, 3)]
    cycle4 = [(0, 1), (1, 2), (2, 3), (3, 0)]
    complete4 = [(i, j) for i in range(4) for j in range(i + 1, 4)]

    return {
        "schema": "matching-one.q1-matching-derivative-defect.v1",
        "exact_object": {
            "amplitude": "Xi(Q,p)=E_p[Q^d]",
            "normalization": "Xi(1,p)=1",
            "first_tangent": "(Q d/dQ) log Xi | Q=1 = E_p[d]",
            "second_tangent": "(Q d/dQ)^2 log Xi | Q=1 = Var_p(d)",
        },
        "edge_FK_planar_dual_control": edge_pair_oracle(
            3,
            triangle,
            2,
            triangle_dual,
            local_term=lambda occupied_edges: 2 - occupied_edges,
        ),
        "edge_FK_nondual_obstruction": edge_pair_oracle(
            3,
            triangle,
            4,
            path4,
            local_term=lambda occupied_edges: 2 - occupied_edges,
        ),
        "site_matching_C4_to_K4": site_pair_oracle(4, cycle4, complete4),
        "S_D_projection_boundary": {
            "formal_exchange": [[0, 1], [1, 0]],
            "projectors_twice_P_plus": [[1, 1], [1, 1]],
            "projectors_twice_P_minus": [[1, -1], [-1, 1]],
            "parity_requires": "a physical two-way identification J with J^2=I",
            "otherwise": "S and D are covectors on a doubled observable pair, not intrinsic fields",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    result = build_oracle()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
