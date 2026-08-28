#!/usr/bin/env python3
"""Exact matching polynomial for a tiny primitive Gaussian square torus.

The square NN graph and its NN+NNN matching graph are built on the quotient
with periods `(a,b),(-b,a)`.  Every site configuration is enumerated.  For a
selected wrapping channel the script records

    q(C) = R_G(C) - R_hat(complement(C))

and sums q over configurations of each occupation count k.  Those integer
Bernstein sums define

    M(p) = sum_k b_k p^k (1-p)^(N-k).

The tool is intentionally bounded to tiny quotients and is an exact algebraic
frontier/reference, not a production Monte Carlo kernel.
"""
from __future__ import annotations

import argparse
import json
from math import comb
from pathlib import Path

from integer_period_torus import classify_configuration, gaussian_integer_torus


CHANNELS = ("cross", "both", "either", "direction_0", "direction_1")


def power_coefficients(bernstein_sums: list[int]) -> list[int]:
    n = len(bernstein_sums) - 1
    coefficients = [0] * (n + 1)
    for k, value in enumerate(bernstein_sums):
        for j in range(n - k + 1):
            coefficients[k + j] += value * ((-1) ** j) * comb(n - k, j)
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return coefficients


def enumerate_gaussian(a: int, b: int, channel: str = "either", max_n: int = 22) -> dict[str, object]:
    if channel not in CHANNELS:
        raise ValueError(f"channel must be one of {CHANNELS}")
    geometry = gaussian_integer_torus(a, b)
    n = geometry.n
    if n > max_n:
        raise ValueError(f"N={n} exceeds exhaustive max_n={max_n}")

    bernstein = [0] * (n + 1)
    primal_counts = [0] * (n + 1)
    matching_complement_counts = [0] * (n + 1)

    for mask in range(1 << n):
        active = tuple(bool(mask & (1 << vertex)) for vertex in range(n))
        complement = tuple(not value for value in active)
        k = bin(mask).count("1")
        primal, _ = classify_configuration(geometry, active, matching=False)
        matching, _ = classify_configuration(geometry, complement, matching=True)
        primal_value = int(getattr(primal, channel))
        matching_value = int(getattr(matching, channel))
        primal_counts[k] += primal_value
        matching_complement_counts[k] += matching_value
        bernstein[k] += primal_value - matching_value

    power = power_coefficients(bernstein)
    return {
        "schema": "matching-one/exact-gaussian-matching-polynomial/v1",
        "geometry": {
            "a": a,
            "b": b,
            "N": n,
            "period_matrix": geometry.periods.matrix,
            "configurations": 1 << n,
        },
        "channel": channel,
        "primal_wrap_counts_by_occupancy": primal_counts,
        "matching_complement_wrap_counts_by_occupancy": matching_complement_counts,
        "bernstein_sums": bernstein,
        "power_coefficients_ascending": power,
        "degree": len(power) - 1,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a", type=int, required=True)
    parser.add_argument("--b", type=int, required=True)
    parser.add_argument("--channel", choices=CHANNELS, default="either")
    parser.add_argument("--max-n", type=int, default=22)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    result = enumerate_gaussian(args.a, args.b, args.channel, args.max_n)
    print(
        f"gaussian ({args.a},{args.b}) N={result['geometry']['N']} "
        f"degree={result['degree']} channel={args.channel}"
    )
    print("Bernstein:", result["bernstein_sums"])
    print("Power:", result["power_coefficients_ascending"])
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
