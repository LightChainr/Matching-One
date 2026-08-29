#!/usr/bin/env python3
"""Tiny exact reference for a fixed-root, landing-marked pivotal H4 counter.

The unmarked count is only a Russo regression.  The new observable is the
axis-minus-diagonal landing mark attached to a pivotal root.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

import mpmath as mp

from exact_pivotal_russo import total_pivotal_mass, wrapping_event
from integer_period_torus import IntegerTorusGeometry, axis_integer_torus


PRIMAL_STEPS = ((1, 0), (-1, 0), (0, 1), (0, -1))
MATCHING_STEPS = PRIMAL_STEPS + (
    (1, 1), (1, -1), (-1, 1), (-1, -1),
)


def _sector(x: int, y: int, registry_shift: int) -> int:
    angle = math.atan2(y, x)
    base = math.floor((angle + math.pi / 8) / (math.pi / 4)) % 8
    return (base - registry_shift) % 8


def _components(
    geometry: IntegerTorusGeometry,
    active: Sequence[bool],
    radius: int,
    *,
    matching: bool,
    enabled_value: bool,
    registry_shift: int,
) -> list[int]:
    steps = MATCHING_STEPS if matching else PRIMAL_STEPS
    points = [
        (x, y)
        for y in range(-radius, radius + 1)
        for x in range(-radius, radius + 1)
        if (x, y) != (0, 0)
    ]
    vertices = {point: geometry.vertex(point) for point in points}
    if len(set(vertices.values())) != len(vertices):
        raise ValueError("local annulus is not injective in the quotient")
    enabled = {
        point for point, vertex in vertices.items()
        if bool(active[vertex]) == enabled_value
    }
    unseen = set(enabled)
    masks: list[int] = []
    inner_steps = set(steps)
    while unseen:
        start = unseen.pop()
        stack = [start]
        component = {start}
        while stack:
            x, y = stack.pop()
            for dx, dy in steps:
                neighbour = (x + dx, y + dy)
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    component.add(neighbour)
                    stack.append(neighbour)
        if not any(point in inner_steps for point in component):
            continue
        mask = 0
        for x, y in component:
            if max(abs(x), abs(y)) == radius:
                mask |= 1 << _sector(x, y, registry_shift)
        if mask:
            masks.append(mask)
    return masks


def _distinct_pair(masks: Sequence[int], first: int, second: int) -> bool:
    return any(
        i != j and masks[i] & (1 << first) and masks[j] & (1 << second)
        for i in range(len(masks)) for j in range(len(masks))
    )


def landing_mark(
    geometry: IntegerTorusGeometry,
    active: Sequence[bool],
    radius: int,
    *,
    open_matching: bool,
    registry_shift: int = 0,
) -> dict[str, int]:
    """Classify two open and two complementary closed landing components."""

    opened = _components(
        geometry, active, radius, matching=open_matching, enabled_value=True,
        registry_shift=registry_shift,
    )
    closed = _components(
        geometry, active, radius, matching=not open_matching, enabled_value=False,
        registry_shift=registry_shift,
    )
    axis = (
        _distinct_pair(opened, 0, 4) and _distinct_pair(closed, 2, 6)
    ) or (
        _distinct_pair(opened, 2, 6) and _distinct_pair(closed, 0, 4)
    )
    diagonal = (
        _distinct_pair(opened, 1, 5) and _distinct_pair(closed, 3, 7)
    ) or (
        _distinct_pair(opened, 3, 7) and _distinct_pair(closed, 1, 5)
    )
    return {
        "axis": int(axis),
        "diagonal": int(diagonal),
        "both": int(axis and diagonal),
        "landed": int(axis or diagonal),
        "h4": int(axis) - int(diagonal),
    }


def pivotal_contribution(
    geometry: IntegerTorusGeometry,
    active: list[bool],
    radius: int,
    *,
    matching: bool,
    registry_shift: int = 0,
) -> dict[str, int]:
    root = geometry.vertex((0, 0))
    active[root] = False
    without = wrapping_event(
        geometry, active, matching=matching, channel="cross"
    )
    active[root] = True
    with_root = wrapping_event(
        geometry, active, matching=matching, channel="cross"
    )
    active[root] = False
    pivotal = with_root - without
    if pivotal not in (0, 1):
        raise AssertionError("cross event is not monotone")
    mark = {
        "axis": 0, "diagonal": 0, "both": 0, "landed": 0, "h4": 0,
    }
    if pivotal:
        mark = landing_mark(
            geometry, active, radius, open_matching=matching,
            registry_shift=registry_shift,
        )
    return {"pivotal": pivotal, **mark}


def marked_pair(
    geometry: IntegerTorusGeometry,
    black: Sequence[bool],
    radius: int,
    registry_shift: int = 0,
) -> dict[str, dict[str, int]]:
    root = geometry.vertex((0, 0))
    primal_active = list(black)
    primal_active[root] = False
    matching_active = [not value for value in black]
    matching_active[root] = False
    return {
        "primal": pivotal_contribution(
            geometry, primal_active, radius, matching=False,
            registry_shift=registry_shift,
        ),
        "matching": pivotal_contribution(
            geometry, matching_active, radius, matching=True,
            registry_shift=registry_shift,
        ),
    }


def _transform_axis(
    geometry: IntegerTorusGeometry,
    active: Sequence[bool],
    transform,
) -> list[bool]:
    output = [False] * geometry.n
    for vertex, (x, y) in enumerate(geometry.coordinates):
        output[geometry.vertex(transform(x, y))] = bool(active[vertex])
    return output


def exact_axis_l4_oracle() -> dict[str, object]:
    """Enumerate 2^15 fixed-root configurations and run symmetry controls."""

    geometry = axis_integer_torus(4)
    root = geometry.vertex((0, 0))
    others = [vertex for vertex in range(geometry.n) if vertex != root]
    totals = {
        side: {key: 0 for key in ("pivotal", "axis", "diagonal", "both", "landed", "h4")}
        for side in ("primal", "matching")
    }
    violations = {"rotation_C4": 0, "reflection": 0, "registry_pi_over_4": 0}
    configurations = 1 << len(others)
    for mask in range(configurations):
        active = [False] * geometry.n
        for bit, vertex in enumerate(others):
            active[vertex] = bool((mask >> bit) & 1)
        result = marked_pair(geometry, active, 1)
        rotated = marked_pair(
            geometry, _transform_axis(geometry, active, lambda x, y: (-y, x)), 1
        )
        reflected = marked_pair(
            geometry, _transform_axis(geometry, active, lambda x, y: (x, -y)), 1
        )
        shifted = marked_pair(geometry, active, 1, registry_shift=1)
        for side in ("primal", "matching"):
            for key in totals[side]:
                totals[side][key] += result[side][key]
            if result[side] != rotated[side]:
                violations["rotation_C4"] += 1
            if result[side] != reflected[side]:
                violations["reflection"] += 1
            if (
                result[side]["h4"] != -shifted[side]["h4"]
                or result[side]["landed"] != shifted[side]["landed"]
            ):
                violations["registry_pi_over_4"] += 1

    mp.mp.dps = 60
    primal_mass = total_pivotal_mass(
        geometry, mp.mpf("0.5"), matching=False, channel="cross"
    )
    matching_mass = total_pivotal_mass(
        geometry, mp.mpf("0.5"), matching=True, channel="cross"
    )
    fixed_root_primal = mp.mpf(geometry.n * totals["primal"]["pivotal"]) / configurations
    fixed_root_matching = mp.mpf(geometry.n * totals["matching"]["pivotal"]) / configurations
    return {
        "schema": "matching-one/marked-pivotal-h4-exact/v1",
        "geometry": "axis-L4",
        "N": geometry.n,
        "radius": 1,
        "probability": "0.5",
        "fixed_root_configurations": configurations,
        "totals": totals,
        "russo_control": {
            "fixed_root_primal_mass": mp.nstr(fixed_root_primal, 40),
            "all_site_primal_mass": mp.nstr(primal_mass, 40),
            "primal_difference": mp.nstr(fixed_root_primal - primal_mass, 20),
            "fixed_root_matching_mass": mp.nstr(fixed_root_matching, 40),
            "all_site_matching_mass": mp.nstr(matching_mass, 40),
            "matching_difference": mp.nstr(fixed_root_matching - matching_mass, 20),
        },
        "symmetry_violations": violations,
        "interpretation": (
            "R=1 is a schema/orientation oracle only; total pivotal mass is a "
            "Russo control, while h4 is the new marked observable."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-axis-l4", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not args.exact_axis_l4:
        raise SystemExit("select --exact-axis-l4")
    payload = exact_axis_l4_oracle()
    text = json.dumps(payload, indent=2) + "\n"
    print(text, end="")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
