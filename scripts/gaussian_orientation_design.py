#!/usr/bin/env python3
"""Generate exactly square torus orientation designs for the square lattice.

A Gaussian-integer period pair

    v1 = (a, b), v2 = (-b, a)

has modulus tau=i, N=a^2+b^2 lattice sites, and microscopic orientation
angle theta=atan2(b,a).  Integers with multiple sum-of-two-squares
representations give fixed-N/fixed-shape orientation comparisons.

This script is design-only: it does not simulate percolation.
"""

from __future__ import annotations

import argparse
import math
from collections import defaultdict
from dataclasses import dataclass


@dataclass(frozen=True)
class Orientation:
    a: int
    b: int

    @property
    def n_sites(self) -> int:
        return self.a * self.a + self.b * self.b

    @property
    def theta_rad(self) -> float:
        return math.atan2(self.b, self.a)

    @property
    def theta_deg(self) -> float:
        return math.degrees(self.theta_rad)

    @property
    def cos4(self) -> float:
        # Exact rational expression evaluated in float.
        a2 = self.a * self.a
        b2 = self.b * self.b
        return (a2 * a2 - 6 * a2 * b2 + b2 * b2) / (a2 + b2) ** 2


def primitive_orientations(max_a: int) -> dict[int, list[Orientation]]:
    by_n: dict[int, list[Orientation]] = defaultdict(list)
    for a in range(1, max_a + 1):
        for b in range(0, a + 1):
            if math.gcd(a, b) != 1:
                continue
            by_n[a * a + b * b].append(Orientation(a, b))
    return by_n


def magic_sequence(count: int) -> list[Orientation]:
    """Convergents to tan(pi/8)=sqrt(2)-1.

    Starting with (a,b)=(2,1) gives a coarse point; the sequence used for
    production begins at (5,2).  Recurrence corresponds to multiplication by
    1+sqrt(2) in the Pell representation.
    """

    result: list[Orientation] = []
    a, b = 2, 1
    while len(result) < count:
        a, b = 2 * a + b, a
        result.append(Orientation(a, b))
    return result


def print_fixed_n(max_a: int, max_n: int, min_reps: int) -> None:
    by_n = primitive_orientations(max_a)
    print("N     (a,b)       theta_deg      cos4theta")
    for n in sorted(by_n):
        reps = by_n[n]
        if n > max_n or len(reps) < min_reps:
            continue
        print(f"\nN={n}  reps={len(reps)}")
        for item in sorted(reps, key=lambda x: x.theta_rad):
            print(
                f"      ({item.a:3d},{item.b:3d})"
                f"   {item.theta_deg:12.7f}"
                f"   {item.cos4:+.10f}"
            )


def print_magic(count: int) -> None:
    print("(a,b)        N        theta_deg       cos4theta")
    for item in magic_sequence(count):
        print(
            f"({item.a:5d},{item.b:5d})"
            f"  {item.n_sites:8d}"
            f"  {item.theta_deg:13.9f}"
            f"  {item.cos4:+.12e}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="mode", required=True)

    fixed = sub.add_parser("fixed-n", help="list N with multiple primitive orientations")
    fixed.add_argument("--max-a", type=int, default=100)
    fixed.add_argument("--max-n", type=int, default=5000)
    fixed.add_argument("--min-reps", type=int, default=2)

    magic = sub.add_parser("magic", help="list pi/8 spin-4-zero approximants")
    magic.add_argument("--count", type=int, default=8)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.mode == "fixed-n":
        print_fixed_n(args.max_a, args.max_n, args.min_reps)
    else:
        print_magic(args.count)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
