#!/usr/bin/env python3
"""Exhaustive central-parity check for the C4 self-matching triangulation.

The checkerboard triangulation in ``c4_self_matching_cyclic_geometry.py`` is
its own site-matching graph.  This script lifts that finite graph to the
Gaussian period lattice, enumerates every occupation mask, and records exact
integer Bernstein coefficients for wrapping probabilities and their
black-minus-complement difference.

This is a tiny-quotient correctness gate, not a production simulator.  The
default ``(a,b)=(3,1)`` quotient has ten sites and 1024 configurations.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
from typing import Sequence

from integer_period_torus import (
    IntegerTorusGeometry,
    classify_configuration,
    gaussian_integer_torus,
)
from matched_torus_reference import Edge


CHANNELS = ("cross", "both", "either", "direction_0", "direction_1")


def c4_self_matching_torus(a: int, b: int) -> IntegerTorusGeometry:
    """Return the lifted checkerboard triangulation for odd Gaussian periods."""

    if a % 2 != 1 or b % 2 != 1:
        raise ValueError("checkerboard parity requires odd a and b")
    base = gaussian_integer_torus(a, b)
    diagonals: list[Edge] = []
    for source, (x, y) in enumerate(base.coordinates):
        # Period-vector coordinate sums are even for odd a,b, so x+y parity
        # is independent of the selected quotient representative.
        if (x + y) % 2 == 0:
            for dx, dy in ((1, 1), (1, -1)):
                diagonals.append(
                    Edge(source, base.vertex((x + dx, y + dy)), dx, dy)
                )
    edges = base.primal_edges + tuple(diagonals)
    if len(edges) != 3 * base.n:
        raise AssertionError("checkerboard triangulation must have 3N edges")
    return replace(
        base,
        name=f"c4-self-matching-{a}-{b}",
        primal_edges=edges,
        matching_edges=edges,
    )


def _active(mask: int, n: int) -> tuple[bool, ...]:
    return tuple(bool(mask & (1 << vertex)) for vertex in range(n))


def _fraction_payload(value: Fraction) -> dict[str, int | str]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "fraction": str(value),
    }


def enumerate_exact(a: int = 3, b: int = 1) -> dict[str, object]:
    """Enumerate the quotient and return exact per-occupancy coefficients."""

    geometry = c4_self_matching_torus(a, b)
    n = geometry.n
    if n > 22:
        raise ValueError(f"N={n} is too large for this exhaustive reference")

    r_coeff = {channel: [0] * (n + 1) for channel in CHANNELS}
    m_coeff = {channel: [0] * (n + 1) for channel in CHANNELS}
    pair_failures = {channel: 0 for channel in CHANNELS}
    invariant_rotation_failures = {channel: 0 for channel in CHANNELS[:3]}

    # The exact C4 action on cover coordinates is (x,y)->(-y,x).
    rotation = [geometry.vertex((-y, x)) for x, y in geometry.coordinates]

    cache: list[dict[str, int]] = []
    for mask in range(1 << n):
        wrapping, _ = classify_configuration(geometry, _active(mask, n))
        cache.append({channel: int(getattr(wrapping, channel)) for channel in CHANNELS})

    full_mask = (1 << n) - 1
    for mask, black in enumerate(cache):
        k = mask.bit_count()
        white = cache[full_mask ^ mask]
        rotated_mask = 0
        for vertex, target in enumerate(rotation):
            if mask & (1 << vertex):
                rotated_mask |= 1 << target
        rotated = cache[rotated_mask]
        for channel in CHANNELS:
            r_coeff[channel][k] += black[channel]
            m_coeff[channel][k] += black[channel] - white[channel]
            complement_q = white[channel] - black[channel]
            if black[channel] - white[channel] != -complement_q:
                pair_failures[channel] += 1
        for channel in invariant_rotation_failures:
            if black[channel] != rotated[channel]:
                invariant_rotation_failures[channel] += 1

    channel_results: dict[str, object] = {}
    for channel in CHANNELS:
        coeff = m_coeff[channel]
        r = r_coeff[channel]
        anti_palindromic = all(coeff[k] == -coeff[n - k] for k in range(n + 1))
        complement_coefficient_identity = all(
            coeff[k] == r[k] - r[n - k] for k in range(n + 1)
        )
        central = Fraction(sum(coeff), 2**n)
        one_third = Fraction(
            sum(coeff[k] * 2 ** (n - k) for k in range(n + 1)), 3**n
        )
        channel_results[channel] = {
            "R_bernstein_integer_coefficients": r,
            "M_bernstein_integer_coefficients": coeff,
            "coefficient_identity_Mk_equals_Rk_minus_RNminusk": (
                complement_coefficient_identity
            ),
            "M_coefficients_anti_palindromic": anti_palindromic,
            "configuration_complement_pair_failures": pair_failures[channel],
            "M_at_p_one_half": _fraction_payload(central),
            "M_at_p_one_third": _fraction_payload(one_third),
            "M_polynomial_identically_zero": all(value == 0 for value in coeff),
        }

    passed = all(
        result["coefficient_identity_Mk_equals_Rk_minus_RNminusk"]
        and result["M_coefficients_anti_palindromic"]
        and result["configuration_complement_pair_failures"] == 0
        and result["M_at_p_one_half"]["numerator"] == 0
        for result in channel_results.values()
    ) and all(value == 0 for value in invariant_rotation_failures.values())
    return {
        "schema": "matching-one/c4-self-matching-exact/v1",
        "geometry": {
            "a": a,
            "b": b,
            "N": n,
            "configurations": 2**n,
            "edges": len(geometry.primal_edges),
            "matching_edges_equal_primal_edges": (
                geometry.matching_edges == geometry.primal_edges
            ),
            "period_matrix": geometry.periods.matrix,
        },
        "c4_invariant_channel_rotation_failures": invariant_rotation_failures,
        "channels": channel_results,
        "passed": passed,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a", type=int, default=3)
    parser.add_argument("--b", type=int, default=1)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args(argv)
    payload = enumerate_exact(args.a, args.b)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
