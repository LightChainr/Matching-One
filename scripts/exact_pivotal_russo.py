#!/usr/bin/env python3
"""Exact finite-volume Russo regression for torus wrapping observables.

For a monotone primal wrapping event ``R_G`` and its matching-lattice partner,

    M(p) = R_G(p) - R_hat(1-p),

the chain rule and Russo's formula give

    M'(p) = sum_v P_p(v pivotal for R_G)
            + sum_v P_(1-p)(v pivotal for R_hat).

This module verifies the identity by independent subset enumeration.  For the
cross channel on N<=9 it also compares with the repository's exact
threshold-rank reconstruction.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Sequence

import mpmath as mp

from analyze_threshold_ranks import matching_derivative
from exact_matching_polynomial import bernstein_to_power, evaluate_power
from integer_period_torus import (
    IntegerTorusGeometry,
    axis_integer_torus,
    classify_configuration,
    diamond_integer_torus,
    gaussian_integer_torus,
)
from threshold_rank_nz import enumerate_exact


CHANNELS = ("cross", "either", "both", "direction_0", "direction_1")


def wrapping_event(
    geometry: IntegerTorusGeometry,
    active: Sequence[bool],
    *,
    matching: bool,
    channel: str,
) -> int:
    if channel not in CHANNELS:
        raise ValueError(f"unsupported wrapping channel {channel!r}")
    channels, _components = classify_configuration(
        geometry, active, matching=matching
    )
    return int(getattr(channels, channel))


def matching_bernstein_counts(
    geometry: IntegerTorusGeometry, channel: str
) -> List[int]:
    """Return exact Bernstein sums for the selected matching difference."""

    if geometry.n > 20:
        raise ValueError("reference subset enumeration is limited to N<=20")
    counts = [0] * (geometry.n + 1)
    for mask in range(1 << geometry.n):
        active = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
        white = [not value for value in active]
        counts[sum(active)] += wrapping_event(
            geometry, active, matching=False, channel=channel
        ) - wrapping_event(geometry, white, matching=True, channel=channel)
    return counts


def total_pivotal_mass(
    geometry: IntegerTorusGeometry,
    probability: mp.mpf,
    *,
    matching: bool,
    channel: str,
) -> mp.mpf:
    """Return ``sum_v P(v is pivotal)`` by exact enumeration of other sites."""

    p = mp.mpf(probability)
    if not 0 <= p <= 1:
        raise ValueError("probability must lie in [0,1]")
    if geometry.n > 16:
        raise ValueError("reference pivotal enumeration is limited to N<=16")

    q = 1 - p
    total = mp.mpf(0)
    for pivot in range(geometry.n):
        other_vertices = [vertex for vertex in range(geometry.n) if vertex != pivot]
        for mask in range(1 << (geometry.n - 1)):
            active = [False] * geometry.n
            occupied = 0
            for bit, vertex in enumerate(other_vertices):
                enabled = bool((mask >> bit) & 1)
                active[vertex] = enabled
                occupied += int(enabled)

            active[pivot] = False
            event_without = wrapping_event(
                geometry, active, matching=matching, channel=channel
            )
            active[pivot] = True
            event_with = wrapping_event(
                geometry, active, matching=matching, channel=channel
            )
            increment = event_with - event_without
            if increment not in (0, 1):
                raise AssertionError(
                    f"{channel} is not monotone at vertex {pivot}: {increment}"
                )
            weight = p**occupied * q ** (geometry.n - 1 - occupied)
            total += increment * weight
    return total


def threshold_rank_cross_derivative(
    geometry: IntegerTorusGeometry, probability: mp.mpf
) -> mp.mpf:
    """Reconstruct the cross-channel derivative from all N! rank paths."""

    if geometry.n > 9:
        raise ValueError("exact threshold-rank comparison is limited to N<=9")
    counts = enumerate_exact(geometry)
    return matching_derivative(
        geometry.n,
        counts.sample_count,
        counts.kminus,
        counts.kplus,
        mp.mpf(probability),
    )


def russo_audit(
    geometry: IntegerTorusGeometry,
    probability: mp.mpf,
    *,
    channel: str = "cross",
    include_threshold_rank: bool = True,
) -> Dict[str, object]:
    p = mp.mpf(probability)
    counts = matching_bernstein_counts(geometry, channel)
    coefficients = bernstein_to_power(counts)
    derivative_coefficients = [
        degree * value for degree, value in enumerate(coefficients)
    ][1:]
    analytic_derivative = evaluate_power(derivative_coefficients, p)
    primal_mass = total_pivotal_mass(
        geometry, p, matching=False, channel=channel
    )
    matching_mass = total_pivotal_mass(
        geometry, 1 - p, matching=True, channel=channel
    )
    pivotal_sum = primal_mass + matching_mass

    result: Dict[str, object] = {
        "geometry": geometry.name,
        "N": geometry.n,
        "channel": channel,
        "p": p,
        "analytic_matching_derivative": analytic_derivative,
        "primal_total_pivotal_mass": primal_mass,
        "matching_total_pivotal_mass_at_complement": matching_mass,
        "russo_pivotal_sum": pivotal_sum,
        "analytic_minus_pivotal": analytic_derivative - pivotal_sum,
    }
    if include_threshold_rank:
        if channel != "cross":
            raise ValueError("threshold-rank comparison is defined only for cross")
        threshold_derivative = threshold_rank_cross_derivative(geometry, p)
        result.update(
            {
                "threshold_rank_derivative": threshold_derivative,
                "threshold_rank_minus_pivotal": threshold_derivative - pivotal_sum,
            }
        )
    return result


def json_ready(value: object) -> object:
    if isinstance(value, mp.mpf):
        return mp.nstr(value, 60)
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_ready(item) for item in value]
    return value


def geometry_from_args(args: argparse.Namespace) -> IntegerTorusGeometry:
    if args.geometry == "axis":
        if args.L is None:
            raise SystemExit("axis geometry requires --L")
        return axis_integer_torus(args.L)
    if args.geometry == "diamond":
        if args.L is None:
            raise SystemExit("diamond geometry requires --L")
        return diamond_integer_torus(args.L)
    if args.a is None or args.b is None:
        raise SystemExit("gaussian geometry requires --a and --b")
    return gaussian_integer_torus(args.a, args.b)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", choices=("axis", "diamond", "gaussian"), required=True)
    parser.add_argument("--L", type=int)
    parser.add_argument("--a", type=int)
    parser.add_argument("--b", type=int)
    parser.add_argument("--p", required=True)
    parser.add_argument("--channel", choices=CHANNELS, default="cross")
    parser.add_argument("--skip-threshold-rank", action="store_true")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    if args.dps < 40:
        raise SystemExit("--dps must be at least 40")
    mp.mp.dps = args.dps
    geometry = geometry_from_args(args)
    result = russo_audit(
        geometry,
        mp.mpf(args.p),
        channel=args.channel,
        include_threshold_rank=not args.skip_threshold_rank,
    )
    payload = json_ready(result)
    rendered = json.dumps(payload, indent=2) + "\n"
    print(rendered, end="")
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
