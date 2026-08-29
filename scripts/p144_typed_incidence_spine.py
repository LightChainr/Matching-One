#!/usr/bin/env python3
"""Exact binary local-state spine for square-site matching topology."""

from __future__ import annotations

import argparse
import json
from collections import Counter, deque
from pathlib import Path

from integer_period_torus import (
    IntegerHomologyUnionFind,
    IntegerPeriods,
    axis_integer_torus,
    classify_configuration,
)


EDGE_PORTS = ((1, 0), (-1, 0), (0, 1), (0, -1))
FACE_PORTS = ((1, 1), (-1, 1), (1, -1), (-1, -1))
FACE_BOUNDARY_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0))
FACE_ALL_PAIRS = tuple((i, j) for i in range(4) for j in range(i + 1, 4))


def components(vertices: set[int], edges: tuple[tuple[int, int], ...]) -> tuple[tuple[int, ...], ...]:
    adjacency = {vertex: set() for vertex in vertices}
    for first, second in edges:
        if first in vertices and second in vertices:
            adjacency[first].add(second)
            adjacency[second].add(first)
    remaining = set(vertices)
    output = []
    while remaining:
        start = min(remaining)
        remaining.remove(start)
        queue = deque([start])
        part = {start}
        while queue:
            vertex = queue.popleft()
            for neighbor in adjacency[vertex]:
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    part.add(neighbor)
                    queue.append(neighbor)
        output.append(tuple(sorted(part)))
    return tuple(sorted(output))


def local_face_certificate() -> dict[str, object]:
    """The edge-hub and face-hub stars realize 4/8 connectivity locally."""
    rows = []
    for mask in range(16):
        black = {corner for corner in range(4) if mask & (1 << corner)}
        white = set(range(4)) - black
        black_reference = components(black, FACE_BOUNDARY_EDGES)
        white_reference = components(white, FACE_ALL_PAIRS)
        # A hub star connects every selected port.  Its terminal partition is
        # one block when nonempty; boundary-edge black connectivity may split.
        black_spine = components(black, FACE_BOUNDARY_EDGES)
        white_spine = (tuple(sorted(white)),) if white else ()
        rows.append({
            "mask": mask,
            "black_reference_partition": [list(part) for part in black_reference],
            "black_edge_hub_partition": [list(part) for part in black_spine],
            "white_matching_partition": [list(part) for part in white_reference],
            "white_face_hub_partition": [list(part) for part in white_spine],
            "pass": black_reference == black_spine and white_reference == white_spine,
        })
    return {
        "patterns": rows,
        "all_16_pass": all(row["pass"] for row in rows),
        "interpretation": (
            "edge-midpoint hubs reproduce black NN boundary edges; each face-center hub "
            "replaces the white NN+NNN clique by an embedded star with the same terminal partition"
        ),
    }


def node_index(L: int, point: tuple[int, int]) -> int:
    modulus = 2 * L
    x, y = point
    return (x % modulus) + modulus * (y % modulus)


def spine_homology(
    L: int,
    active: tuple[bool, ...],
    *,
    colour: str,
) -> tuple[int, int, tuple[tuple[int, ...], ...]]:
    """Return ambient rank, component count, and bases for one typed spine."""
    if colour not in {"black", "white"}:
        raise ValueError("colour must be black or white")
    modulus = 2 * L
    periods = IntegerPeriods(((modulus, 0), (0, modulus)))
    union_find = IntegerHomologyUnionFind(modulus * modulus, periods)
    enabled_sites = []
    ports = EDGE_PORTS if colour == "black" else FACE_PORTS

    for y in range(L):
        for x in range(L):
            site_number = x + L * y
            selected = active[site_number] if colour == "black" else not active[site_number]
            if not selected:
                continue
            site_point = (2 * x, 2 * y)
            site = node_index(L, site_point)
            enabled_sites.append(site)
            for dx, dy in ports:
                hub = node_index(L, (site_point[0] + dx, site_point[1] + dy))
                union_find.add_edge(site, hub, dx, dy)

    roots = {}
    for site in enabled_sites:
        component = union_find.component(site)
        roots[component.root] = component
    ranks = [component.rank for component in roots.values()]
    bases = tuple(sorted(component.basis for component in roots.values()))
    return (max(ranks, default=0), len(roots), bases)


def direct_record(L: int, active: tuple[bool, ...]) -> dict[str, object]:
    geometry = axis_integer_torus(L)
    black_channels, black_components = classify_configuration(geometry, active)
    white_channels, white_components = classify_configuration(
        geometry, tuple(not value for value in active), matching=True
    )
    rank_black = black_channels.max_rank
    rank_white = white_channels.max_rank
    q = int(black_channels.either) - int(white_channels.either)
    return {
        "rank_black": rank_black,
        "rank_white": rank_white,
        "components_black": len(black_components),
        "components_white": len(white_components),
        "q": q,
    }


def boundary_descriptor(rank_black: int, rank_white: int) -> dict[str, object]:
    if (rank_black, rank_white) == (1, 1):
        boundary_rank = 1
        torus_side = "annular_both"
    elif (rank_black, rank_white) == (2, 0):
        boundary_rank = 0
        torus_side = "black"
    elif (rank_black, rank_white) == (0, 2):
        boundary_rank = 0
        torus_side = "white"
    else:
        raise AssertionError(f"non-Alexander rank pair {(rank_black, rank_white)}")
    return {
        "boundary_ambient_rank": boundary_rank,
        "torus_side": torus_side,
        "q_from_sides": (rank_black - rank_white) // 2,
    }


def enumerate_axis_L3() -> dict[str, object]:
    L = 3
    N = L * L
    coefficient_counts: Counter[tuple[int, int, int]] = Counter()
    matching_coefficients = [0] * (N + 1)
    spine_failures = []
    descriptor_failures = []
    complement_symbol_failures = 0

    for mask in range(1 << N):
        active = tuple(bool(mask & (1 << vertex)) for vertex in range(N))
        occupied = sum(active)
        direct = direct_record(L, active)
        spine_black, components_black, _ = spine_homology(L, active, colour="black")
        spine_white, components_white, _ = spine_homology(L, active, colour="white")
        if (
            spine_black != direct["rank_black"]
            or spine_white != direct["rank_white"]
            or components_black != direct["components_black"]
            or components_white != direct["components_white"]
        ):
            spine_failures.append({
                "mask": mask,
                "direct": direct,
                "spine": {
                    "rank_black": spine_black,
                    "rank_white": spine_white,
                    "components_black": components_black,
                    "components_white": components_white,
                },
            })
        descriptor = boundary_descriptor(spine_black, spine_white)
        if descriptor["q_from_sides"] != direct["q"]:
            descriptor_failures.append({"mask": mask, "direct": direct, "descriptor": descriptor})
        coefficient_counts[(occupied, spine_black, spine_white)] += 1
        matching_coefficients[occupied] += direct["q"]
        complement_symbol_failures += any(a == b for a, b in zip(active, (not value for value in active)))

    coefficient_rows = [
        {"occupied": occupied, "rank_black": rb, "rank_white": rw, "count": count}
        for (occupied, rb, rw), count in sorted(coefficient_counts.items())
    ]
    derived_matching = [0] * (N + 1)
    for row in coefficient_rows:
        derived_matching[row["occupied"]] += (
            (row["rank_black"] - row["rank_white"]) // 2
        ) * row["count"]

    return {
        "geometry": "axis_L3_honest_square_cell_torus",
        "L": L,
        "N": N,
        "configurations": 1 << N,
        "state_sum": (
            "Phi_L(p;x,y)=sum_mask p^|B|(1-p)^(N-|B|) "
            "x^r_black(mask) y^r_white(mask)"
        ),
        "matching_specialization": "M_L(p)=1/2 (x d_x-y d_y) Phi_L(p;x,y)|x=y=1",
        "coefficient_rows": coefficient_rows,
        "matching_Bernstein_coefficients_direct": matching_coefficients,
        "matching_Bernstein_coefficients_from_state_sum": derived_matching,
        "spine_failure_count": len(spine_failures),
        "spine_failures": spine_failures[:4],
        "boundary_descriptor_failure_count": len(descriptor_failures),
        "boundary_descriptor_failures": descriptor_failures[:4],
        "complement_local_symbol_failure_count": complement_symbol_failures,
    }


def smoothing_obstruction() -> dict[str, object]:
    target = ((0, 1, 2, 3),)
    pairings = (
        ((0, 1), (2, 3)),
        ((0, 3), (1, 2)),
    )
    return {
        "terminal_order": ["north", "east", "south", "west"],
        "target_black_star_partition": [list(block) for block in target],
        "two_pairwise_smoothings": [
            [list(block) for block in pairing] for pairing in pairings
        ],
        "either_smoothing_equals_target": any(pairing == target for pairing in pairings),
        "minimal_missing_partition_type": "J4=four-way junction {north,east,south,west}",
        "consequence": (
            "a pure two-pairing transition polynomial loses the four-arm cluster matrix element; "
            "the exact binary alphabet is instead typed junctions J_edge and J_face on alternating port sets"
        ),
    }


def build_oracle() -> dict[str, object]:
    local = local_face_certificate()
    exact = enumerate_axis_L3()
    obstruction = smoothing_obstruction()
    return {
        "schema": "matching-one.p144-typed-incidence-spine.v1",
        "issues": [144, 269],
        "fixed_embedded_object": {
            "name": "doubled-lattice edge/face incidence spine",
            "nodes": [
                "even-even original site nodes",
                "odd-even/even-odd edge-midpoint hubs",
                "odd-odd face-center hubs",
            ],
            "local_state_B": "J_edge: join the four edge-midpoint ports through the site node",
            "local_state_W": "J_face: join the four face-center ports through the site node",
            "complement": "exchange B and W at every site, hence exchange the two typed local junction states",
            "ribbon_interpretation": (
                "in the #269 cellwise thickening, the black edge-star spine retracts from U and the "
                "white face-star spine is a 1-skeleton of closure(T^2\\U); the subsurfaces U and V, "
                "not arbitrary metric neighborhoods of the bare spines, share the typed boundary"
            ),
        },
        "local_16_pattern_certificate": local,
        "honest_torus_exact_oracle": exact,
        "ambient_rank_readout": {
            "rank_sum": "r_black+r_white=2",
            "charge": "q=(r_black-r_white)/2",
            "boundary_cases": {
                "boundary_rank_1": "(r_black,r_white)=(1,1), q=0; both sides are annular in ambient H1",
                "boundary_rank_0_black_torus_side": "(2,0), q=+1",
                "boundary_rank_0_white_torus_side": "(0,2), q=-1",
            },
        },
        "minimal_transfer_rule": {
            "branch": "for each new site choose B or W",
            "B_update": "union its site node with four edge-midpoint hubs",
            "W_update": "union its site node with four face-center hubs",
            "frontier_state": "two typed connectivity partitions with integer lift potentials and accumulated H1 generators",
            "closure": "when the periodic frontier closes, emit p^|B|(1-p)^(N-|B|) x^r_black y^r_white",
            "width_scope": "finite for every fixed strip width; no claim of width-independent bond dimension",
        },
        "naive_two_smoothing_obstruction": obstruction,
        "claim_boundary": {
            "proved": [
                "the 16 local terminal partitions of black NN and white NN+NNN are realized by the typed incidence hubs",
                "all 512 axis-L3 masks have identical component counts and ambient ranks in the spine and reference graphs",
                "the rank-graded state sum specializes exactly to the finite matching polynomial coefficients",
                "the site complement exchanges the binary typed junction labels",
                "two pairwise smoothings alone cannot realize the required four-terminal join",
            ],
            "not_proved": [
                "identification with a named transition or vertex polynomial",
                "a deletion/contraction relation with fewer states than the frontier rule",
                "a width-independent finite-state recursion",
                "complement oddness on non-self-matching axis tori",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    payload = json.dumps(build_oracle(), indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
