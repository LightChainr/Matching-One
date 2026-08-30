#!/usr/bin/env python3
"""Exact mapping gate for the separated axis/diagonal arm observer.

This script contains no Monte Carlo data.  It verifies that the frozen R=6
anchors and their local eight-neighbour rings are distinct in each production
quotient, and records the exact two-orbit scalar/spin-4 inverse.
"""

from __future__ import annotations

import argparse
import json
from math import hypot
from pathlib import Path
from typing import Any

from integer_period_torus import IntegerPeriods, determinant


MATRICES = {
    "N325_first": ((17, -6), (6, 17)),
    "N325_second": ((18, -1), (1, 18)),
    "N425_first": ((16, -13), (13, 16)),
    "N425_second": ((19, -8), (8, 19)),
}
LOCAL_RING = (
    (1, 0), (1, 1), (0, 1), (-1, 1),
    (-1, 0), (-1, -1), (0, -1), (1, -1),
)


def rotate(point: tuple[int, int], turns: int) -> tuple[int, int]:
    x, y = point
    for _ in range(turns % 4):
        x, y = -y, x
    return x, y


def quotient_distance(periods: IntegerPeriods, point: tuple[int, int]) -> float:
    """Shortest Euclidean lift, with an exact finite search for these small P."""
    (a, b), (c, d) = periods.matrix
    return min(
        hypot(point[0] + a * u + b * v, point[1] + c * u + d * v)
        for u in range(-3, 4)
        for v in range(-3, 4)
    )


def geometry_gate(matrix: tuple[tuple[int, int], tuple[int, int]], radius: int) -> dict[str, Any]:
    periods = IntegerPeriods(matrix)
    origin = periods.quotient_key((0, 0))
    axis = [rotate((radius, 0), turn) for turn in range(4)]
    diagonal = [rotate((radius, radius), turn) for turn in range(4)]
    axis_keys = [periods.quotient_key(point) for point in axis]
    diagonal_keys = [periods.quotient_key(point) for point in diagonal]

    ring_injective = True
    root_excluded = True
    ring_keys: list[tuple[int, int]] = []
    for anchor in axis + diagonal:
        keys = [
            periods.quotient_key((anchor[0] + dx, anchor[1] + dy))
            for dx, dy in LOCAL_RING
        ]
        ring_injective &= len(set(keys)) == len(LOCAL_RING)
        root_excluded &= origin not in keys and periods.quotient_key(anchor) != origin
        ring_keys.extend(keys)

    # Translation covariance is tested on every quotient representative key:
    # q(root+d) is q(root)+q(d), hence the relative anchor key is root-free.
    translations = [(x, y) for x in range(periods.order) for y in (0,)]
    translation_checks = 0
    for x, y in translations:
        root = (x, y)
        root_key = periods.quotient_key(root)
        for displacement in axis + diagonal:
            lhs = periods.quotient_key((root[0] + displacement[0], root[1] + displacement[1]))
            displacement_key = periods.quotient_key(displacement)
            rhs = (
                (root_key[0] + displacement_key[0]) % periods.order,
                (root_key[1] + displacement_key[1]) % periods.order,
            )
            translation_checks += 1
            if lhs != rhs:
                raise AssertionError("translation covariance failed")

    result = {
        "N": periods.order,
        "determinant": determinant(matrix),
        "matrix": matrix,
        "radius": radius,
        "axis_orbit": axis,
        "diagonal_orbit": diagonal,
        "axis_orbit_size": len(set(axis_keys)),
        "diagonal_orbit_size": len(set(diagonal_keys)),
        "cross_orbit_disjoint": set(axis_keys).isdisjoint(diagonal_keys),
        "local_ring_injective": ring_injective,
        "source_root_excluded_from_all_rings": root_excluded,
        "translation_checks": translation_checks,
        "minimum_axis_quotient_distance": min(quotient_distance(periods, point) for point in axis),
        "minimum_diagonal_quotient_distance": min(quotient_distance(periods, point) for point in diagonal),
        "typed_complex": "axis_landing + i*diagonal_landing",
        "spatial_h4": "axis_landing - diagonal_landing",
    }
    assert result["axis_orbit_size"] == 4
    assert result["diagonal_orbit_size"] == 4
    assert result["cross_orbit_disjoint"]
    assert result["local_ring_injective"]
    assert result["source_root_excluded_from_all_rings"]
    return result


def build_artifact(radius: int = 6) -> dict[str, Any]:
    response = ((1, 1), (1, -1))
    det = response[0][0] * response[1][1] - response[0][1] * response[1][0]
    assert det == -2
    return {
        "schema": "matching-one/separated-two-orbit-observer-gate/v1",
        "data_class": "exact quotient mapping and character algebra only",
        "radius": radius,
        "rotation_stream": "one counter-derived common C4 turn per replica and pre-insertion k",
        "response_matrix": response,
        "response_determinant": det,
        "response_rank": 2,
        "scalar_projection": "(axis+diagonal)/2",
        "spin4_projection": "(axis-diagonal)/2",
        "geometries": {
            name: geometry_gate(matrix, radius) for name, matrix in MATRICES.items()
        },
        "boundary": (
            "The gate authorizes the observer semantics and removes the one-orbit alias. "
            "It does not assert a nonzero coupling, continuum field identity, or scaling law."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius", type=int, default=6)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = build_artifact(args.radius)
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
