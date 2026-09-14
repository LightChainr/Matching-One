#!/usr/bin/env python3
"""Exact tiny-torus controls for the rank-one projective slope harmonic.

The observable is defined only on rank-one configurations. For a primitive
unoriented winding line (a,b) on the square torus, report

    Z4 = ((a+i b)^4)/(a^2+b^2)^2.

Two exact weightings are returned:

1. p=1/2, where every configuration has equal weight;
2. the integral over p in [0,1].  A k-site configuration then has beta weight
   1 / ((N+1) * binom(N,k)), so this equals the plateau-duration weighting
   E[(T2-T1) Z4] without enumerating site orders.

The script also checks the configurationwise NN/complementary-matching
rank-one line identity. This is a finite control, not a continuum
extrapolation.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb, gcd
from pathlib import Path


def canonical_line(a: int, b: int) -> tuple[int, int]:
    divisor = gcd(abs(a), abs(b))
    if divisor == 0:
        raise ValueError("zero vector has no projective line")
    a //= divisor
    b //= divisor
    if a < 0 or (a == 0 and b < 0):
        a, b = -a, -b
    return a, b


def z4_components(line: tuple[int, int]) -> tuple[Fraction, Fraction]:
    a, b = line
    denominator = (a * a + b * b) ** 2
    real = Fraction(a**4 - 6 * a * a * b * b + b**4, denominator)
    imag = Fraction(4 * a * b * (a * a - b * b), denominator)
    return real, imag


def graph_rank_line(mask: int, L: int, *, matching: bool) -> tuple[int, tuple[int, int] | None]:
    N = L * L
    steps = [(1, 0), (0, 1)]
    if matching:
        steps += [(1, 1), (1, -1)]

    def vid(x: int, y: int) -> int:
        return (y % L) * L + (x % L)

    def xy(vertex: int) -> tuple[int, int]:
        return vertex % L, vertex // L

    occupied = [bool(mask >> vertex & 1) for vertex in range(N)]
    adjacency: list[list[tuple[int, tuple[int, int]]]] = [[] for _ in range(N)]
    for u in range(N):
        if not occupied[u]:
            continue
        x, y = xy(u)
        for dx, dy in steps:
            v = vid(x + dx, y + dy)
            if occupied[v]:
                adjacency[u].append((v, (dx, dy)))
                adjacency[v].append((u, (-dx, -dy)))

    seen = [False] * N
    gains: list[tuple[int, int]] = []
    for root in range(N):
        if not occupied[root] or seen[root]:
            continue
        seen[root] = True
        lift = {root: (0, 0)}
        stack = [root]
        while stack:
            u = stack.pop()
            ux, uy = lift[u]
            for v, (dx, dy) in adjacency[u]:
                candidate = (ux + dx, uy + dy)
                if not seen[v]:
                    seen[v] = True
                    lift[v] = candidate
                    stack.append(v)
                else:
                    vx, vy = lift[v]
                    defect = (candidate[0] - vx, candidate[1] - vy)
                    if defect != (0, 0):
                        if defect[0] % L or defect[1] % L:
                            raise AssertionError("non-period lift defect")
                        gain = (defect[0] // L, defect[1] // L)
                        if gain != (0, 0):
                            gains.append(gain)

    if not gains:
        return 0, None
    a, b = gains[0]
    for c, d in gains[1:]:
        if a * d - b * c != 0:
            return 2, None
    return 1, canonical_line(a, b)


def _fraction_map(values: dict[tuple[int, int], Fraction]) -> dict[str, str]:
    return {f"{a},{b}": str(value) for (a, b), value in sorted(values.items())}


def run(L: int) -> dict[str, object]:
    N = L * L
    if N > 16:
        raise ValueError("This exact control intentionally caps at L<=4.")
    full = (1 << N) - 1
    counts: Counter[tuple[int, int]] = Counter()
    cardinality_counts: Counter[int] = Counter()
    rank_one = 0
    z4_real_sum = Fraction(0)
    z4_imag_sum = Fraction(0)
    beta_rank_one = Fraction(0)
    beta_z4_real = Fraction(0)
    beta_z4_imag = Fraction(0)
    beta_line_weights: dict[tuple[int, int], Fraction] = defaultdict(Fraction)
    complement_violations = 0

    for mask in range(1 << N):
        rank4, line4 = graph_rank_line(mask, L, matching=False)
        if rank4 == 1:
            rank_one += 1
            assert line4 is not None
            counts[line4] += 1
            occupied_count = mask.bit_count()
            cardinality_counts[occupied_count] += 1
            real, imag = z4_components(line4)
            z4_real_sum += real
            z4_imag_sum += imag

            beta_weight = Fraction(1, (N + 1) * comb(N, occupied_count))
            beta_rank_one += beta_weight
            beta_z4_real += beta_weight * real
            beta_z4_imag += beta_weight * imag
            beta_line_weights[line4] += beta_weight

            rank8, line8 = graph_rank_line(full ^ mask, L, matching=True)
            if rank8 != 1 or line8 != line4:
                complement_violations += 1

    if complement_violations:
        raise AssertionError(f"rank-one complement slope violations: {complement_violations}")

    total = 1 << N
    conditional_real = z4_real_sum / rank_one
    conditional_imag = z4_imag_sum / rank_one
    unnormalized_real = z4_real_sum / total
    unnormalized_imag = z4_imag_sum / total
    gap_conditional_real = beta_z4_real / beta_rank_one
    gap_conditional_imag = beta_z4_imag / beta_rank_one

    return {
        "L": L,
        "N": N,
        "p": "1/2",
        "configurations_checked": total,
        "rank_one_configurations": rank_one,
        "rank_one_cardinality_counts": {str(k): count for k, count in sorted(cardinality_counts.items())},
        "slope_counts": {f"{a},{b}": count for (a, b), count in sorted(counts.items())},
        "z4_unnormalized_real": str(unnormalized_real),
        "z4_unnormalized_real_float": float(unnormalized_real),
        "z4_unnormalized_imag": str(unnormalized_imag),
        "z4_conditional_real": str(conditional_real),
        "z4_conditional_real_float": float(conditional_real),
        "z4_conditional_imag": str(conditional_imag),
        "integrated_rank_one_probability": str(beta_rank_one),
        "integrated_z4_real": str(beta_z4_real),
        "integrated_z4_imag": str(beta_z4_imag),
        "gap_duration_weighted_z4_real": str(gap_conditional_real),
        "gap_duration_weighted_z4_real_float": float(gap_conditional_real),
        "gap_duration_weighted_z4_imag": str(gap_conditional_imag),
        "integrated_line_weights": _fraction_map(beta_line_weights),
        "rank_one_NN_complement_matching_line_identity": True,
        "complement_violations": 0,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--L", type=int, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(args.L)
    text = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)
