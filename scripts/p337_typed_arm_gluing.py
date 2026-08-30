#!/usr/bin/env python3
"""Exact fixtures for globally typed six/eight-arm surgery.

The mathematical gluing statement is in
notes/p337-typed-arm-gluing-20260830.md.  This script certifies its two local
templates, a contractible ordinary-six-arm counterexample, and the smallest
exact theta/figure-eight quotient witnesses.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import deque
from pathlib import Path
from typing import Iterable

try:
    from p334_birth_age_collision_review_20260830 import enumerate_states
    from p337_direct_birth_arm_topology import carrier_descriptor
except ModuleNotFoundError:
    from scripts.p334_birth_age_collision_review_20260830 import enumerate_states
    from scripts.p337_direct_birth_arm_topology import carrier_descriptor


DEFAULT_OUTPUT = Path("results/exact-typed-arm-gluing/latest.json")
NN = ((1, 0), (-1, 0), (0, 1), (0, -1))


def mod_point(point: tuple[int, int], side: int) -> tuple[int, int]:
    return point[0] % side, point[1] % side


def occupied_components(side: int, occupied: set[tuple[int, int]]) -> list[set[tuple[int, int]]]:
    unseen = set(occupied)
    components = []
    while unseen:
        start = unseen.pop()
        component = {start}
        queue = deque([start])
        while queue:
            x, y = queue.popleft()
            for dx, dy in NN:
                other = mod_point((x + dx, y + dy), side)
                if other in unseen:
                    unseen.remove(other)
                    component.add(other)
                    queue.append(other)
        components.append(component)
    return components


def ambient_rank(side: int, occupied: set[tuple[int, int]]) -> int:
    """Ambient H1 rank of one square-torus site set, using exact lifts."""
    occupied = {mod_point(point, side) for point in occupied}
    lifted: dict[tuple[int, int], tuple[int, int]] = {}
    cycle_vectors: list[tuple[int, int]] = []
    for start in occupied:
        if start in lifted:
            continue
        lifted[start] = (0, 0)
        queue = deque([start])
        while queue:
            point = queue.popleft()
            x, y = lifted[point]
            for dx, dy in NN:
                other = mod_point((point[0] + dx, point[1] + dy), side)
                if other not in occupied:
                    continue
                candidate = (x + dx, y + dy)
                if other not in lifted:
                    lifted[other] = candidate
                    queue.append(other)
                else:
                    vector = (candidate[0] - lifted[other][0], candidate[1] - lifted[other][1])
                    if vector != (0, 0):
                        cycle_vectors.append(vector)
    if not cycle_vectors:
        return 0
    x0, y0 = cycle_vectors[0]
    return 2 if any(x0 * y - y0 * x for x, y in cycle_vectors[1:]) else 1


def matching_path_is_valid(
    path: list[tuple[int, int]],
    occupied: set[tuple[int, int]],
) -> bool:
    if any(point in occupied for point in path) or len(path) != len(set(path)):
        return False
    return all(
        max(abs(a[0] - b[0]), abs(a[1] - b[1])) == 1
        for a, b in zip(path, path[1:])
    )


def ordinary_six_arm_counterexample(radius: int = 4, side: int = 13) -> dict:
    """Three open and three matching-vacant arms with no homology birth."""
    if side <= 2 * radius + 2:
        raise ValueError("the arm box must embed in the torus")
    open_arms = [
        [(x, 0) for x in range(1, radius + 1)],
        [(0, y) for y in range(1, radius + 1)],
        [(-x, 0) for x in range(1, radius + 1)],
    ]
    vacant_arms = [
        [(d, d) for d in range(1, radius + 1)],
        [(-d, d) for d in range(1, radius + 1)],
        [(0, -d) for d in range(1, radius + 1)],
    ]
    occupied_integer = {point for arm in open_arms for point in arm}
    if not all(matching_path_is_valid(path, occupied_integer) for path in vacant_arms):
        raise AssertionError("vacant separator fixture failed")
    occupied = {mod_point(point, side) for point in occupied_integer}
    before_rank = ambient_rank(side, occupied)
    before_components = len(occupied_components(side, occupied))
    after = occupied | {(0, 0)}
    after_rank = ambient_rank(side, after)
    after_components = len(occupied_components(side, after))
    if (before_rank, before_components, after_rank, after_components) != (0, 3, 0, 1):
        raise AssertionError("ordinary six-arm counterexample changed")
    return {
        "side": side,
        "radius": radius,
        "open_arms": [[[x, y] for x, y in arm] for arm in open_arms],
        "vacant_matching_arms": [[[x, y] for x, y in arm] for arm in vacant_arms],
        "before_rank": before_rank,
        "before_components": before_components,
        "after_adding_origin_rank": after_rank,
        "after_adding_origin_components": after_components,
        "typed_failure": "occupied landings have three component IDs and no deck-address determinant",
    }


def local_template(kind: str) -> dict:
    """Return a fixed 5x5 surgery block; the central site is not sampled."""
    block = {(x, y) for x in range(-2, 3) for y in range(-2, 3)}
    centre = (0, 0)
    if kind == "theta":
        gates = [(3, 0), (0, 3), (-3, 0)]
        paths = [
            [(2, 0), (1, 0)],
            [(0, 2), (0, 1)],
            [(-2, 0), (-1, 0)],
        ]
        separators = [
            [(1, 1), (2, 1), (2, 2)],
            [(-1, 1), (-2, 1), (-2, 2)],
            [(0, -1), (0, -2)],
        ]
    elif kind == "figure_eight":
        gates = [(3, 0), (0, 3), (-3, 0), (0, -3)]
        paths = [
            [(2, 0), (1, 0)],
            [(0, 2), (0, 1)],
            [(-2, 0), (-1, 0)],
            [(0, -2), (0, -1)],
        ]
        separators = [
            [(1, 1), (2, 1), (2, 2)],
            [(-1, 1), (-2, 1), (-2, 2)],
            [(-1, -1), (-2, -1), (-2, -2)],
            [(1, -1), (2, -1), (2, -2)],
        ]
    else:
        raise ValueError("kind must be theta or figure_eight")

    open_sites = {point for path in paths for point in path}
    closed_sites = block - open_sites - {centre}
    if sum(len(component) for component in _finite_components(open_sites)) != len(open_sites):
        raise AssertionError("template component accounting failed")
    if len(_finite_components(open_sites)) != len(paths):
        raise AssertionError("local corridors touch each other")
    if not all(matching_path_is_valid(path, open_sites) for path in separators):
        raise AssertionError("local vacant separators failed")
    if len(open_sites) + len(closed_sites) != 24:
        raise AssertionError("all noncentral block sites must be forced")
    return {
        "kind": kind,
        "block": "[-2,2]^2 minus the birth site",
        "gates": [list(point) for point in gates],
        "open_corridors": [[[x, y] for x, y in path] for path in paths],
        "vacant_separator_seeds": [[[x, y] for x, y in path] for path in separators],
        "forced_open": len(open_sites),
        "forced_closed": len(closed_sites),
        "forced_total": 24,
        "bernoulli_weight": f"p^{len(open_sites)}*(1-p)^{len(closed_sites)}",
        "uniform_finite_energy_bound": "eta^24 for p in [eta,1-eta]",
    }


def _finite_components(points: Iterable[tuple[int, int]]) -> list[set[tuple[int, int]]]:
    unseen = set(points)
    answer = []
    while unseen:
        start = unseen.pop()
        component = {start}
        queue = deque([start])
        while queue:
            x, y = queue.popleft()
            for dx, dy in NN:
                other = x + dx, y + dy
                if other in unseen:
                    unseen.remove(other)
                    component.add(other)
                    queue.append(other)
        answer.append(component)
    return answer


def quotient_fixtures() -> dict:
    n10, states10 = enumerate_states(3, 1)
    theta_mask, theta_site = 122, 7
    theta = carrier_descriptor(3, 1, theta_mask, theta_site)
    if states10[theta_mask][0] != 0 or states10[theta_mask | (1 << theta_site)][0] != 2:
        raise AssertionError("N10 theta fixture rank changed")
    if theta["topology"] != "one_carrier_theta":
        raise AssertionError("N10 theta fixture type changed")

    n9, states9 = enumerate_states(3, 0)
    figure_mask, figure_site = 30, 0
    figure = carrier_descriptor(3, 0, figure_mask, figure_site)
    if states9[figure_mask][0] != 0 or states9[figure_mask | 1][0] != 2:
        raise AssertionError("N9 figure-eight fixture rank changed")
    if figure["topology"] != "two_carrier_figure_eight":
        raise AssertionError("N9 figure-eight fixture type changed")
    return {
        "theta": {
            "generator": [3, 1],
            "N": n10,
            "old_mask": theta_mask,
            "birth_site": theta_site,
            "descriptor": theta,
        },
        "figure_eight": {
            "generator": [3, 0],
            "N": n9,
            "old_mask": figure_mask,
            "birth_site": figure_site,
            "descriptor": figure,
        },
    }


def build_result() -> dict:
    return {
        "schema": "p337-typed-arm-gluing-v1",
        "status": "deterministic_gluing_with_exact_fixtures",
        "ordinary_six_arm_counterexample": ordinary_six_arm_counterexample(),
        "local_surgery_templates": {
            "theta": local_template("theta"),
            "figure_eight": local_template("figure_eight"),
        },
        "small_quotient_fixtures": quotient_fixtures(),
        "minimal_typed_fields": {
            "theta": [
                "ambient_rank_zero",
                "common_occupied_component_id_for_three_landings",
                "two_relative_deck_address_vectors_with_nonzero_determinant",
                "fixed_cyclic_landing_word_with_three_matching_vacant_separators",
            ],
            "figure_eight": [
                "ambient_rank_zero",
                "two_pair_occupied_component_partition",
                "one_nonzero_deck_difference_per_pair_with_nonzero_mutual_determinant",
                "fixed_cyclic_landing_word_with_four_matching_vacant_separators",
            ],
        },
        "finite_energy": {
            "theta_fixed_gate_lower_factor": "p^6*(1-p)^18",
            "figure_eight_fixed_gate_lower_factor": "p^8*(1-p)^16",
            "uniform_factor_on_p_in_eta_1_minus_eta": "eta^24",
            "scale_dependence": "none; the surgery block has 24 sampled sites",
        },
        "claim_boundary": (
            "The reverse comparison holds for globally typed, landing-separated events. "
            "Ordinary six arms do not determine torus homology. No arm exponent, "
            "universality transfer, or asymptotic figure-eight suppression is claimed."
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
