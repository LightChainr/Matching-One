#!/usr/bin/env python3
"""Exact local complement-odd pivotal-H4 readout on the C4 control.

The observable is half the difference between the fixed-root pivotal H4
landing mark of a configuration and its occupation complement.  It is a
local, complement-odd readout intended to complement the global wrapping
readout in the self-matching ``(t, lambda)`` tangent plane.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Sequence

from c4_self_matching_exact import c4_self_matching_torus
from integer_period_torus import IntegerTorusGeometry, classify_configuration


Point = tuple[int, int]


def _active(mask: int, n: int) -> list[bool]:
    return [bool(mask & (1 << vertex)) for vertex in range(n)]


def _cross(geometry: IntegerTorusGeometry, active: Sequence[bool]) -> int:
    return int(classify_configuration(geometry, tuple(active))[0].cross)


def _sector(x: int, y: int, registry_shift: int = 0) -> int:
    angle = math.atan2(y, x)
    base = math.floor((angle + math.pi / 8) / (math.pi / 4)) % 8
    return (base - registry_shift) % 8


def _local_adjacent(first: Point, second: Point) -> bool:
    x, y = first
    dx, dy = second[0] - x, second[1] - y
    if abs(dx) + abs(dy) == 1:
        return True
    # Checkerboard diagonals join the even-parity corners of each square.
    return abs(dx) == abs(dy) == 1 and (x + y) % 2 == 0


def _component_masks(
    geometry: IntegerTorusGeometry,
    active: Sequence[bool],
    radius: int,
    enabled_value: bool,
    registry_shift: int,
) -> list[int]:
    points = [
        (x, y)
        for y in range(-radius, radius + 1)
        for x in range(-radius, radius + 1)
        if (x, y) != (0, 0)
    ]
    vertices = {point: geometry.vertex(point) for point in points}
    if len(set(vertices.values())) != len(vertices):
        raise ValueError("local annulus is not injective in the quotient")
    adjacency = {
        point: [other for other in points if _local_adjacent(point, other)]
        for point in points
    }
    root_neighbours = {
        point for point in points if _local_adjacent((0, 0), point)
    }
    unseen = {
        point
        for point, vertex in vertices.items()
        if bool(active[vertex]) == enabled_value
    }
    masks: list[int] = []
    while unseen:
        start = unseen.pop()
        stack = [start]
        component = {start}
        while stack:
            point = stack.pop()
            for neighbour in adjacency[point]:
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    component.add(neighbour)
                    stack.append(neighbour)
        if not (component & root_neighbours):
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
        for i in range(len(masks))
        for j in range(len(masks))
    )


def landing_h4(
    geometry: IntegerTorusGeometry,
    active: Sequence[bool],
    radius: int = 1,
    registry_shift: int = 0,
) -> int:
    opened = _component_masks(geometry, active, radius, True, registry_shift)
    closed = _component_masks(geometry, active, radius, False, registry_shift)
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
    return int(axis) - int(diagonal)


def pivotal_h4(
    geometry: IntegerTorusGeometry,
    active: Sequence[bool],
    radius: int = 1,
    registry_shift: int = 0,
) -> int:
    root = geometry.vertex((0, 0))
    work = list(active)
    work[root] = False
    without = _cross(geometry, work)
    work[root] = True
    with_root = _cross(geometry, work)
    work[root] = False
    pivotal = with_root - without
    if pivotal not in (0, 1):
        raise AssertionError("cross event is not monotone")
    return pivotal * landing_h4(geometry, work, radius, registry_shift)


def _transform(
    geometry: IntegerTorusGeometry,
    active: Sequence[bool],
    transform,
) -> list[bool]:
    output = [False] * geometry.n
    for vertex, (x, y) in enumerate(geometry.coordinates):
        output[geometry.vertex(transform(x, y))] = bool(active[vertex])
    return output


def _fraction(value: Fraction) -> str:
    return str(value)


def exact_n10_report() -> dict[str, object]:
    geometry = c4_self_matching_torus(3, 1)
    n = geometry.n
    scale = 1 << n
    totals = {"global_t": 0, "global_lambda": 0, "local_t": 0, "local_lambda": 0}
    local_values: Counter[int] = Counter()
    violations = {"complement_odd": 0, "rotation_C4": 0, "registry_pi_over_4": 0}
    full_mask = scale - 1
    for mask in range(scale):
        black = _active(mask, n)
        white = _active(full_mask ^ mask, n)
        global_twice = _cross(geometry, black) - _cross(geometry, white)
        black_local = pivotal_h4(geometry, black)
        white_local = pivotal_h4(geometry, white)
        local_twice = black_local - white_local
        local_values[local_twice] += 1

        score_t = 2 * sum(1 if value else -1 for value in black)
        score_lambda = 2 * sum(
            (1 if vertex % 2 == 0 else -1) * (1 if value else -1)
            for vertex, value in enumerate(black)
        )
        totals["global_t"] += global_twice * score_t
        totals["global_lambda"] += global_twice * score_lambda
        totals["local_t"] += local_twice * score_t
        totals["local_lambda"] += local_twice * score_lambda

        complement_local_twice = pivotal_h4(geometry, white) - pivotal_h4(geometry, black)
        if complement_local_twice != -local_twice:
            violations["complement_odd"] += 1
        transformed = _transform(geometry, black, lambda x, y: (-y, x))
        transformed_white = [not value for value in transformed]
        transformed_twice = pivotal_h4(geometry, transformed) - pivotal_h4(
            geometry, transformed_white
        )
        if transformed_twice != local_twice:
            violations["rotation_C4"] += 1
        shifted_twice = pivotal_h4(geometry, black, registry_shift=1) - pivotal_h4(
            geometry, white, registry_shift=1
        )
        if shifted_twice != -local_twice:
            violations["registry_pi_over_4"] += 1

    # Both observables are defined as half differences, so divide accumulated
    # score products by 2 * 2^N.
    response = [
        [Fraction(totals["global_t"], 2 * scale), Fraction(totals["global_lambda"], 2 * scale)],
        [Fraction(totals["local_t"], 2 * scale), Fraction(totals["local_lambda"], 2 * scale)],
    ]
    determinant = response[0][0] * response[1][1] - response[0][1] * response[1][0]
    passed = (
        response[0] == [Fraction(15, 8), Fraction(5, 4)]
        and determinant != 0
        and not any(violations.values())
    )
    return {
        "schema": "matching-one/c4-local-odd-pivotal/v1",
        "geometry": {"a": 3, "b": 1, "N": n, "configurations": scale, "radius": 1},
        "observable_rows": ["global_cross_half_difference", "local_pivotal_h4_half_difference"],
        "parameter_columns": ["t", "lambda"],
        "response_matrix": [[_fraction(value) for value in row] for row in response],
        "determinant": _fraction(determinant),
        "local_twice_observable_counts": {str(key): value for key, value in sorted(local_values.items())},
        "symmetry_violations": violations,
        "passed": passed,
        "interpretation": (
            "The exact local complement-odd row is linearly independent of the global wrapping row. "
            "This is an oracle and conditioning gate, not a large-N exponent measurement."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    payload = exact_n10_report()
    text = json.dumps(payload, indent=2) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
