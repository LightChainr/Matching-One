#!/usr/bin/env python3
"""Exact carrier classification for direct ambient-H1 births.

For a rank-zero occupied set S and a vacant site v, lift every occupied
neighbour of v to the copy adjacent to a fixed lift of v.  Within each old
occupied component, those adjacent lifts belong to deck translates of one
chosen lifted component.  Differences of their deck addresses generate
exactly the homology created by adding v.

This bounded oracle exhausts small Gaussian square tori.  It certifies the
theta/two-carrier dichotomy; it does not estimate an arm exponent.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import Any

try:
    from p334_birth_age_collision_review_20260830 import enumerate_states
except ModuleNotFoundError:  # imported as scripts.p337_... from the repo root
    from scripts.p334_birth_age_collision_review_20260830 import enumerate_states


DEFAULT_OUTPUT = Path("results/exact-direct-birth-arm-topology/latest.json")
STEPS = ((1, 0), (-1, 0), (0, 1), (0, -1))


def quotient_graph(a: int, b: int) -> tuple[int, list[tuple[int, int]], list[list[tuple[int, int, int]]]]:
    n = a * a + b * b

    def key(x: int, y: int) -> tuple[int, int]:
        return ((a * x + b * y) % n, (-b * x + a * y) % n)

    representatives = [(0, 0)]
    ids = {key(0, 0): 0}
    for x, y in representatives:
        for dx, dy in STEPS:
            point = key(x + dx, y + dy)
            if point not in ids:
                ids[point] = len(representatives)
                representatives.append((x + dx, y + dy))
    if len(representatives) != n:
        raise AssertionError("quotient traversal did not find every vertex")
    neighbours = [
        [(ids[key(x + dx, y + dy)], dx, dy) for dx, dy in STEPS]
        for x, y in representatives
    ]
    return n, representatives, neighbours


def vector_rank(vectors: list[tuple[int, int]]) -> int:
    nonzero = [(x, y) for x, y in vectors if x or y]
    if not nonzero:
        return 0
    x0, y0 = nonzero[0]
    if any(x0 * y - y0 * x for x, y in nonzero[1:]):
        return 2
    return 1


def carrier_descriptor(
    a: int,
    b: int,
    mask: int,
    birth_site: int,
) -> dict[str, Any]:
    """Return the exact universal-cover address data at one vacant site.

    The caller is responsible for checking that the old occupied set has
    ambient rank zero.  Lift consistency below is itself a check of that fact.
    """
    n, representatives, neighbours = quotient_graph(a, b)
    occupied = {u for u in range(n) if (mask >> u) & 1}
    if birth_site in occupied:
        raise ValueError("birth_site is already occupied")

    component: dict[int, int] = {}
    lift_coordinate: dict[int, tuple[int, int]] = {}
    component_id = 0
    for start in occupied:
        if start in component:
            continue
        component[start] = component_id
        lift_coordinate[start] = (0, 0)
        queue = deque([start])
        while queue:
            u = queue.popleft()
            x, y = lift_coordinate[u]
            for other, dx, dy in neighbours[u]:
                if other not in occupied:
                    continue
                candidate = (x + dx, y + dy)
                if other not in component:
                    component[other] = component_id
                    lift_coordinate[other] = candidate
                    queue.append(other)
                else:
                    if component[other] != component_id or lift_coordinate[other] != candidate:
                        raise AssertionError("old component has nonzero ambient homology")
        component_id += 1

    incidence: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for neighbour, dx, dy in neighbours[birth_site]:
        if neighbour not in occupied:
            continue
        x, y = lift_coordinate[neighbour]
        # The adjacent lift is (dx,dy); the chosen component lift puts the
        # same quotient vertex at (x,y).  Their difference is a deck vector.
        address = (dx - x, dy - y)
        incidence[component[neighbour]].append(
            {
                "vertex": neighbour,
                "step": [dx, dy],
                "deck_address": list(address),
            }
        )

    groups = []
    all_differences: list[tuple[int, int]] = []
    for cid, items in sorted(incidence.items()):
        addresses = [tuple(item["deck_address"]) for item in items]
        origin = addresses[0]
        differences = [(x - origin[0], y - origin[1]) for x, y in addresses[1:]]
        affine_rank = vector_rank(differences)
        all_differences.extend(differences)
        multiplicities = sorted(Counter(addresses).values())
        groups.append(
            {
                "component": cid,
                "incidences": len(items),
                "distinct_deck_addresses": len(set(addresses)),
                "address_multiplicities": multiplicities,
                "affine_rank": affine_rank,
                "items": items,
            }
        )

    created_rank = vector_rank(all_differences)
    rank_two_groups = [group for group in groups if group["affine_rank"] == 2]
    if created_rank != 2:
        topology = "not_direct_rank2"
        arm_lower_bound = 0
    elif rank_two_groups:
        # Three affinely independent addresses are three distinct lifted
        # clusters, each forced to reach a nonzero deck translate.
        topology = "one_carrier_theta"
        arm_lower_bound = 3
    else:
        rank_one = [group for group in groups if group["affine_rank"] == 1]
        if len(rank_one) != 2 or any(group["incidences"] != 2 for group in rank_one):
            raise AssertionError("degree-four dichotomy failed")
        directions = []
        for group in rank_one:
            p, q = (tuple(item["deck_address"]) for item in group["items"])
            directions.append((q[0] - p[0], q[1] - p[1]))
        if vector_rank(directions) != 2:
            raise AssertionError("two carrier directions are not independent")
        topology = "two_carrier_figure_eight"
        arm_lower_bound = 4

    return {
        "topology": topology,
        "created_rank": created_rank,
        "occupied_degree": sum(len(items) for items in incidence.values()),
        "occupied_arm_lower_bound": arm_lower_bound,
        "groups": groups,
        "birth_site": birth_site,
        "birth_representative": list(representatives[birth_site]),
    }


def signature(descriptor: dict[str, Any]) -> str:
    parts = []
    for group in descriptor["groups"]:
        parts.append(
            f"{group['incidences']}i/{group['distinct_deck_addresses']}a/"
            f"r{group['affine_rank']}/m{','.join(map(str, group['address_multiplicities']))}"
        )
    return "+".join(sorted(parts))


def geometry_census(a: int, b: int) -> dict[str, Any]:
    n, states = enumerate_states(a, b)
    topology_counts: Counter[str] = Counter()
    signature_counts: Counter[str] = Counter()
    direct_edges = 0
    for mask, (old_rank, _) in enumerate(states):
        if old_rank != 0:
            continue
        for vertex in range(n):
            if (mask >> vertex) & 1:
                continue
            if states[mask | (1 << vertex)][0] != 2:
                continue
            descriptor = carrier_descriptor(a, b, mask, vertex)
            if descriptor["created_rank"] != 2:
                raise AssertionError("cover address rank disagrees with exhaustive homology rank")
            direct_edges += 1
            topology_counts[descriptor["topology"]] += 1
            signature_counts[signature(descriptor)] += 1
    return {
        "a": a,
        "b": b,
        "N": n,
        "direct_edges": direct_edges,
        "topology_counts": dict(sorted(topology_counts.items())),
        "signature_counts": dict(sorted(signature_counts.items())),
    }


def build_result() -> dict[str, Any]:
    geometries = [geometry_census(a, b) for a, b in ((2, 1), (2, 2), (3, 0), (3, 1), (3, 2), (4, 0), (4, 1))]
    expected = {
        (2, 1): (0, 0, 0),
        (2, 2): (40, 40, 0),
        (3, 0): (45, 36, 9),
        (3, 1): (80, 80, 0),
        (3, 2): (793, 793, 0),
        (4, 0): (4624, 4288, 336),
        (4, 1): (8823, 8704, 119),
    }
    for row in geometries:
        total, theta, figure_eight = expected[row["a"], row["b"]]
        if row["direct_edges"] != total:
            raise AssertionError("direct-edge reference count changed")
        if row["topology_counts"].get("one_carrier_theta", 0) != theta:
            raise AssertionError("theta reference count changed")
        if row["topology_counts"].get("two_carrier_figure_eight", 0) != figure_eight:
            raise AssertionError("figure-eight reference count changed")

    # Smallest figure-eight in this Gaussian census: the four neighbours of
    # the origin on the 3x3 torus.  Horizontal and vertical predecessor edges
    # are separate components and carry independent deck differences.
    witness = carrier_descriptor(3, 0, 30, 0)
    if witness["topology"] != "two_carrier_figure_eight":
        raise AssertionError("N=9 figure-eight witness changed")
    witness["old_mask"] = 30
    witness["occupied_vertices"] = [1, 2, 3, 4]

    return {
        "schema": "p337-direct-birth-arm-topology-v1",
        "status": "exact_finite_volume_plus_deterministic_lemma",
        "geometries": geometries,
        "minimal_figure_eight_witness": witness,
        "theorem": {
            "cover_address_identity": "new ambient H1 is generated by within-component deck-address differences",
            "degree_four_dichotomy": ["one_carrier_theta", "two_carrier_figure_eight"],
            "theta_occupied_arm_lower_bound": 3,
            "figure_eight_occupied_arm_lower_bound": 4,
            "annular_separator_condition": "on an embedded square annulus, planar site-matching separation supplies one vacant NN+NNN arm between consecutive lifted occupied carriers",
        },
        "claim_boundary": (
            "Direct birth implies at least a polychromatic six-arm event on embedded annuli; "
            "the converse, probability comparability, six-arm exponent transfer, and a nonzero "
            "scaling amplitude are not proved."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--stdout", action="store_true")
    args = parser.parse_args()
    result = build_result()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    if args.stdout:
        sys.stdout.write(text)


if __name__ == "__main__":
    main()
