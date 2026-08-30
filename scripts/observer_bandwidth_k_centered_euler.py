#!/usr/bin/env python3
"""Exact K-centered Euler/Walsh control on a periodic square-cell torus."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Mapping

from observer_bandwidth_product_walsh import (
    expectation,
    fraction,
    popcount,
    source_fixture,
    walsh_degree_coefficients,
)


def site(x: int, y: int, length: int) -> int:
    return (y % length) * length + (x % length)


def square_torus_incidence(length: int) -> tuple[list[tuple[int, int]], list[tuple[int, ...]]]:
    if length < 3:
        raise ValueError("distinct-edge/four-vertex square-cell contract requires L>=3")
    edges = set()
    faces = []
    for y in range(length):
        for x in range(length):
            origin = site(x, y, length)
            for neighbor in (site(x + 1, y, length), site(x, y + 1, length)):
                edges.add(tuple(sorted((origin, neighbor))))
            face = (
                origin,
                site(x + 1, y, length),
                site(x, y + 1, length),
                site(x + 1, y + 1, length),
            )
            if len(set(face)) != 4:
                raise AssertionError("face vertices are not distinct")
            faces.append(face)
    return sorted(edges), faces


def falling(value: int, order: int) -> int:
    result = 1
    for offset in range(order):
        result *= value - offset
    return result


def conditional_mean(n: int, k: int) -> Fraction:
    return (
        Fraction(k)
        - 2 * n * Fraction(falling(k, 2), falling(n, 2))
        + n * Fraction(falling(k, 4), falling(n, 4))
    )


def euler_observer_values(length: int) -> tuple[dict[int, Fraction], dict[str, int]]:
    n = length * length
    edges, faces = square_torus_incidence(length)
    values = {}
    for mask in range(1 << n):
        occupied_edges = sum(
            ((mask >> left) & 1) and ((mask >> right) & 1)
            for left, right in edges
        )
        full_faces = sum(
            all((mask >> vertex) & 1 for vertex in face) for face in faces
        )
        values[mask] = Fraction(popcount(mask) - occupied_edges + full_faces)
    return values, {"sites": n, "edges": len(edges), "faces": len(faces)}


def k_center(values: Mapping[int, Fraction], n: int) -> dict[int, Fraction]:
    return {
        mask: values[mask] - conditional_mean(n, popcount(mask))
        for mask in range(1 << n)
    }


def conditional_checks(
    values: Mapping[int, Fraction], centered: Mapping[int, Fraction], n: int
) -> list[dict[str, object]]:
    rows = []
    for k in range(n + 1):
        masks = [mask for mask in range(1 << n) if popcount(mask) == k]
        empirical = sum(values[mask] for mask in masks) / len(masks)
        centered_sum = sum(centered[mask] for mask in masks)
        expected = conditional_mean(n, k)
        if empirical != expected or centered_sum != 0:
            raise AssertionError("conditional Euler centering failed")
        rows.append(
            {
                "k": k,
                "configurations": len(masks),
                "enumerated_mean": str(empirical),
                "formula_mean": str(expected),
                "centered_sum": str(centered_sum),
            }
        )
    return rows


def degree_one_projections(
    values: Mapping[int, Fraction], n: int, p: Fraction
) -> list[Fraction]:
    return [
        expectation(
            {
                mask: values[mask] * (Fraction((mask >> index) & 1) - p)
                for mask in range(1 << n)
            },
            n,
            p,
        )
        for index in range(n)
    ]


def build_report(manifest: Mapping[str, object]) -> dict[str, object]:
    length = int(manifest["L"])
    n = length * length
    values, incidence = euler_observer_values(length)
    if incidence != {"sites": n, "edges": 2 * n, "faces": n}:
        raise AssertionError("square-torus incidence count differs from contract")
    centered = k_center(values, n)
    conditionals = conditional_checks(values, centered, n)
    p_rows = []
    for p_value in manifest["p_values"]:
        p = fraction(p_value)
        original = degree_one_projections(values, n, p)
        residual = degree_one_projections(centered, n, p)
        if len(set(original)) != 1:
            raise AssertionError("Euler degree-one projection is not translation invariant")
        if any(residual):
            raise AssertionError("K-centered Euler retains a degree-one projection")
        coefficients = walsh_degree_coefficients(centered, source_fixture(n), n, p)
        if coefficients[1] != 0 or any(coefficients[j] for j in range(5, n + 1)):
            raise AssertionError("K-centered Euler left the degree-two-to-four envelope")
        p_rows.append(
            {
                "p": str(p),
                "original_site_projection": str(original[0]),
                "centered_site_projections": [str(value) for value in residual],
                "source_covariance_degree_coefficients": {
                    str(j): str(coefficients[j]) for j in range(1, n + 1)
                },
                "active_degrees": [j for j, value in coefficients.items() if value],
            }
        )
    return {
        "schema": manifest["schema"],
        "status": "exact_k_centered_degree_one_zero_verified",
        "L": length,
        "incidence": incidence,
        "enumerated_configurations": 1 << n,
        "conditional_rows": conditionals,
        "p_checks": p_rows,
        "boundary": manifest["boundary"],
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "analysis/observer_bandwidth_k_centered_euler_manifest.json",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_report(json.loads(args.manifest.read_text(encoding="utf-8")))
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
