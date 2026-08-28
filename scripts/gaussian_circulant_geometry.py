#!/usr/bin/env python3
"""Reference utilities for primitive Gaussian-integer square tori.

For coprime integers a,b the quotient with periods (a,b),(-b,a) has
N=a^2+b^2 vertices and is canonically represented by Z/NZ with

    j = a*x + b*y (mod N).

The NN square graph is the circulant with steps +/-a,+/-b; the site-matching
graph adds +/-(a+b), +/-(a-b).

This module is intentionally small and correctness-oriented.  Production
connectivity/homology code can reuse the formulas without depending on it.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class LiftedStep:
    residue: int
    dx: int
    dy: int
    name: str


@dataclass(frozen=True)
class GaussianTorus:
    a: int
    b: int

    def __post_init__(self) -> None:
        if self.a <= 0 or self.b < 0:
            raise ValueError("require a>0 and b>=0")
        if math.gcd(self.a, self.b) != 1:
            raise ValueError("this cyclic reference requires gcd(a,b)=1")

    @property
    def n(self) -> int:
        return self.a * self.a + self.b * self.b

    @property
    def theta(self) -> float:
        return math.atan2(self.b, self.a)

    def label(self, x: int, y: int) -> int:
        return (self.a * x + self.b * y) % self.n

    def primal_forward_steps(self) -> tuple[LiftedStep, ...]:
        return (
            LiftedStep(self.a % self.n, 1, 0, "+x"),
            LiftedStep(self.b % self.n, 0, 1, "+y"),
        )

    def matching_forward_steps(self) -> tuple[LiftedStep, ...]:
        return self.primal_forward_steps() + (
            LiftedStep((self.a + self.b) % self.n, 1, 1, "+x+y"),
            LiftedStep((self.a - self.b) % self.n, 1, -1, "+x-y"),
        )

    def winding_coordinates(self, dx: int, dy: int) -> tuple[int, int]:
        """Convert a closed lifted displacement into the period basis.

        Raises ValueError when (dx,dy) is not a period-lattice vector.
        """

        num_m = self.a * dx + self.b * dy
        num_n = -self.b * dx + self.a * dy
        if num_m % self.n or num_n % self.n:
            raise ValueError("displacement is not a closed torus winding")
        return num_m // self.n, num_n // self.n

    def edge_residues(self, matching: bool = False) -> set[int]:
        steps = self.matching_forward_steps() if matching else self.primal_forward_steps()
        out: set[int] = set()
        for step in steps:
            out.add(step.residue % self.n)
            out.add((-step.residue) % self.n)
        return out


def units(n: int) -> Iterable[int]:
    for value in range(1, n):
        if math.gcd(value, n) == 1:
            yield value


def best_multiplier(
    first: GaussianTorus,
    second: GaussianTorus,
    *,
    matching_weight: float = 1.0,
) -> tuple[int, float, int, int]:
    """Graph-only CRN relabeling score for two same-N orientations.

    Returns `(t, score, nn_directed_overlap, matching_directed_overlap)` where
    multiplication by the unit t maps first-geometry vertex labels into the
    second coupling.  This is a structural preselection, not an empirical
    covariance optimizer.
    """

    if first.n != second.n:
        raise ValueError("multiplier coupling requires equal N")
    n = first.n
    nn1 = first.edge_residues(False)
    nn2 = second.edge_residues(False)
    ma1 = first.edge_residues(True)
    ma2 = second.edge_residues(True)
    best: tuple[float, int, int, int] | None = None
    for t in units(n):
        mapped_nn = {(t * s) % n for s in nn1}
        mapped_ma = {(t * s) % n for s in ma1}
        nn_overlap = len(mapped_nn & nn2)
        ma_overlap = len(mapped_ma & ma2)
        score = nn_overlap + matching_weight * ma_overlap
        candidate = (score, -t, nn_overlap, ma_overlap)
        if best is None or candidate > best:
            best = candidate
    assert best is not None
    score, neg_t, nn_overlap, ma_overlap = best
    return -neg_t, score, nn_overlap, ma_overlap


def _self_test() -> None:
    examples = [GaussianTorus(8, 1), GaussianTorus(7, 4), GaussianTorus(12, 5)]
    for g in examples:
        # Both declared periods must vanish under the cyclic label.
        assert g.label(g.a, g.b) == 0
        assert g.label(-g.b, g.a) == 0
        # Fundamental period windings must recover unit basis vectors.
        assert g.winding_coordinates(g.a, g.b) == (1, 0)
        assert g.winding_coordinates(-g.b, g.a) == (0, 1)
        # Unit x/y moves have the advertised residues.
        assert g.label(1, 0) == g.a % g.n
        assert g.label(0, 1) == g.b % g.n

    g1, g2 = GaussianTorus(8, 1), GaussianTorus(7, 4)
    t, _score, nn_overlap, _matching_overlap = best_multiplier(g1, g2)
    assert math.gcd(t, g1.n) == 1
    assert nn_overlap >= 2  # at least one undirected direction can be aligned


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("a", type=int)
    parser.add_argument("b", type=int)
    parser.add_argument("--compare", nargs=2, type=int, metavar=("A2", "B2"))
    args = parser.parse_args()

    g = GaussianTorus(args.a, args.b)
    print(f"N={g.n} theta_deg={math.degrees(g.theta):.9f}")
    print("primal directed residues:", sorted(g.edge_residues(False)))
    print("matching directed residues:", sorted(g.edge_residues(True)))
    if args.compare:
        other = GaussianTorus(*args.compare)
        t, score, nn_overlap, matching_overlap = best_multiplier(g, other)
        print(
            f"best structural multiplier t={t} score={score:g} "
            f"nn_overlap={nn_overlap} matching_overlap={matching_overlap}"
        )
    return 0


if __name__ == "__main__":
    _self_test()
    raise SystemExit(main())
