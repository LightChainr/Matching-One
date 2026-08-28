#!/usr/bin/env python3
"""Exact cyclic covering maps for primitive Gaussian square tori.

For parent (a,b) with N=a^2+b^2 and primitive child (A,B) with
A^2+B^2=Q*N, search units t modulo N for which

    f(k) = t*k mod N

maps the child's square-lattice NN residues bijectively to the parent's NN
residues. The same map then covers the NN+NNN matching graph.

This is a correctness-first reference for Issue #67. Production coupling code
must preserve its graph/fiber contract.
"""

from __future__ import annotations

from dataclasses import dataclass
import argparse
import math
from typing import Iterable


@dataclass(frozen=True)
class GaussianPair:
    a: int
    b: int

    def __post_init__(self) -> None:
        if self.a <= 0 or self.b < 0:
            raise ValueError("require a>0 and b>=0")
        if math.gcd(self.a, self.b) != 1:
            raise ValueError("cyclic covering reference requires gcd(a,b)=1")

    @property
    def n(self) -> int:
        return self.a * self.a + self.b * self.b

    @property
    def nn_signed(self) -> tuple[int, int, int, int]:
        n = self.n
        return (self.a % n, (-self.a) % n, self.b % n, (-self.b) % n)

    @property
    def matching_signed(self) -> tuple[int, ...]:
        n = self.n
        raw = (
            self.a,
            -self.a,
            self.b,
            -self.b,
            self.a + self.b,
            -(self.a + self.b),
            self.a - self.b,
            -(self.a - self.b),
        )
        return tuple(value % n for value in raw)


def units(n: int) -> Iterable[int]:
    for value in range(1, n):
        if math.gcd(value, n) == 1:
            yield value


def signed_direction(value: int, geometry: GaussianPair) -> str:
    n = geometry.n
    lookup = {
        geometry.a % n: "+x",
        (-geometry.a) % n: "-x",
        geometry.b % n: "+y",
        (-geometry.b) % n: "-y",
    }
    if value % n not in lookup:
        raise ValueError("residue is not a parent NN direction")
    return lookup[value % n]


@dataclass(frozen=True)
class CoveringMap:
    parent: GaussianPair
    child: GaussianPair
    t: int

    def __post_init__(self) -> None:
        if self.child.n % self.parent.n:
            raise ValueError("child size must be an integer multiple of parent size")
        if math.gcd(self.t, self.parent.n) != 1:
            raise ValueError("covering multiplier t must be a unit modulo parent N")
        if set(self.map_child_nn_residues()) != set(self.parent.nn_signed):
            raise ValueError("t does not map child NN directions bijectively to parent")
        if set(self.map_child_matching_residues()) != set(self.parent.matching_signed):
            raise ValueError("t does not map child matching directions bijectively to parent")

    @property
    def degree(self) -> int:
        return self.child.n // self.parent.n

    @property
    def inverse_t(self) -> int:
        return pow(self.t, -1, self.parent.n)

    def map_vertex(self, child_label: int) -> int:
        return (self.t * child_label) % self.parent.n

    def map_child_nn_residues(self) -> tuple[int, ...]:
        n = self.parent.n
        return tuple((self.t * residue) % n for residue in self.child.nn_signed)

    def map_child_matching_residues(self) -> tuple[int, ...]:
        n = self.parent.n
        return tuple((self.t * residue) % n for residue in self.child.matching_signed)

    def direction_map(self) -> dict[str, str]:
        child_positive = (("+x", self.child.a), ("+y", self.child.b))
        return {
            name: signed_direction(self.t * residue, self.parent)
            for name, residue in child_positive
        }

    def fiber(self, parent_label: int) -> tuple[int, ...]:
        n = self.parent.n
        child_n = self.child.n
        base = (self.inverse_t * parent_label) % n
        return tuple((base + m * n) % child_n for m in range(self.degree))

    def inverse_coordinates(self, child_label: int) -> tuple[int, int]:
        """Return `(parent_label, kernel_index)` for one child cyclic label."""

        parent = self.map_vertex(child_label)
        n = self.parent.n
        base = (self.inverse_t * parent) % n
        difference = (child_label - base) % self.child.n
        if difference % n:
            raise AssertionError("fiber arithmetic lost divisibility by parent N")
        kernel_index = (difference // n) % self.degree
        return parent, kernel_index

    def verify_partition(self) -> None:
        seen: set[int] = set()
        for parent_label in range(self.parent.n):
            fiber = self.fiber(parent_label)
            if len(set(fiber)) != self.degree:
                raise AssertionError("fiber has duplicate child labels")
            for child_label in fiber:
                if child_label in seen:
                    raise AssertionError("fibers overlap")
                if self.map_vertex(child_label) != parent_label:
                    raise AssertionError("fiber element maps to wrong parent")
                seen.add(child_label)
        if seen != set(range(self.child.n)):
            raise AssertionError("fibers do not partition child vertices")

    def verify_edges(self) -> None:
        parent_nn = set(self.parent.nn_signed)
        parent_matching = set(self.parent.matching_signed)
        for child_label in range(self.child.n):
            parent_label = self.map_vertex(child_label)
            for residue in self.child.nn_signed:
                mapped_step = (
                    self.map_vertex(child_label + residue) - parent_label
                ) % self.parent.n
                if mapped_step not in parent_nn:
                    raise AssertionError("NN edge does not map to parent NN edge")
            for residue in self.child.matching_signed:
                mapped_step = (
                    self.map_vertex(child_label + residue) - parent_label
                ) % self.parent.n
                if mapped_step not in parent_matching:
                    raise AssertionError(
                        "matching edge does not map to parent matching edge"
                    )


def covering_units(parent: GaussianPair, child: GaussianPair) -> tuple[int, ...]:
    if child.n % parent.n:
        return ()
    parent_nn = set(parent.nn_signed)
    parent_matching = set(parent.matching_signed)
    found: list[int] = []
    for t in units(parent.n):
        mapped_nn = {(t * residue) % parent.n for residue in child.nn_signed}
        if mapped_nn != parent_nn:
            continue
        mapped_matching = {
            (t * residue) % parent.n for residue in child.matching_signed
        }
        if mapped_matching != parent_matching:
            continue
        found.append(t)
    return tuple(found)


def canonical_cover(parent: GaussianPair, child: GaussianPair) -> CoveringMap:
    candidates = covering_units(parent, child)
    if not candidates:
        raise ValueError("no cyclic Gaussian covering unit found")
    # Purely graph/arithmetic convention. Production variance pilots may compare
    # the finite candidate set, but must freeze that choice before evaluation.
    return CoveringMap(parent, child, min(candidates))


def _self_test() -> None:
    examples = (
        (GaussianPair(8, 1), GaussianPair(9, 7), 2),
        (GaussianPair(7, 4), GaussianPair(11, 3), 2),
        (GaussianPair(8, 1), GaussianPair(17, 6), 5),
        (GaussianPair(7, 4), GaussianPair(18, 1), 5),
        (GaussianPair(9, 2), GaussianPair(16, 13), 5),
        (GaussianPair(7, 6), GaussianPair(19, 8), 5),
    )
    for parent, child, degree in examples:
        candidates = covering_units(parent, child)
        assert len(candidates) == 4  # the four D4 direction identifications
        cover = canonical_cover(parent, child)
        assert cover.degree == degree
        assert set(cover.direction_map()) == {"+x", "+y"}
        assert len(set(cover.direction_map().values())) == 2
        cover.verify_partition()
        cover.verify_edges()
        for child_label in range(child.n):
            parent_label, kernel_index = cover.inverse_coordinates(child_label)
            assert child_label in cover.fiber(parent_label)
            assert 0 <= kernel_index < degree


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent_a", type=int)
    parser.add_argument("parent_b", type=int)
    parser.add_argument("child_a", type=int)
    parser.add_argument("child_b", type=int)
    args = parser.parse_args()

    parent = GaussianPair(args.parent_a, args.parent_b)
    child = GaussianPair(args.child_a, args.child_b)
    candidates = covering_units(parent, child)
    if not candidates:
        raise SystemExit("no covering units")
    print(f"parent N={parent.n} child N={child.n} degree={child.n // parent.n}")
    print("covering units:", " ".join(map(str, candidates)))
    for t in candidates:
        cover = CoveringMap(parent, child, t)
        print(f"t={t} direction_map={cover.direction_map()}")
    cover = canonical_cover(parent, child)
    cover.verify_partition()
    cover.verify_edges()
    print(f"canonical t={cover.t}; exact partition/edge checks PASS")
    return 0


if __name__ == "__main__":
    _self_test()
    raise SystemExit(main())
