#!/usr/bin/env python3
"""Exact cyclic covering maps for primitive Gaussian square tori.

For parent (a,b) with N=a^2+b^2 and primitive child (A,B) with
A^2+B^2=Q*N, search units t modulo N for which

    f(k) = t*k mod N

maps the child's square-lattice NN residues bijectively to the parent's NN
residues. The same map then covers the NN+NNN matching graph.

The local step map is a D4 signed permutation R, but the induced map on torus
winding coordinates is the integer matrix

    H = P_parent^{-1} R P_child,

where P(a,b) has columns (a,b),(-b,a). It satisfies |det H|=Q. Thus homology
rank/cross-wrapping is preserved, while individual winding directions need not
map by a mere D4 permutation.

This is a correctness-first reference for Issue #67. Production coupling code
must preserve its graph/fiber/homology contract.
"""

from __future__ import annotations

from dataclasses import dataclass
import argparse
import math
from typing import Iterable


Matrix2 = tuple[tuple[int, int], tuple[int, int]]
Vector2 = tuple[int, int]


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
    def period_matrix(self) -> Matrix2:
        # Columns are (a,b) and (-b,a).
        return ((self.a, -self.b), (self.b, self.a))

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


def direction_vector(name: str) -> Vector2:
    vectors = {
        "+x": (1, 0),
        "-x": (-1, 0),
        "+y": (0, 1),
        "-y": (0, -1),
    }
    return vectors[name]


def matmul(left: Matrix2, right: Matrix2) -> Matrix2:
    return (
        (
            left[0][0] * right[0][0] + left[0][1] * right[1][0],
            left[0][0] * right[0][1] + left[0][1] * right[1][1],
        ),
        (
            left[1][0] * right[0][0] + left[1][1] * right[1][0],
            left[1][0] * right[0][1] + left[1][1] * right[1][1],
        ),
    )


def matvec(matrix: Matrix2, vector: Vector2) -> Vector2:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1],
    )


def det(matrix: Matrix2) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


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
        self.verify_homology()

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

    def direction_matrix(self) -> Matrix2:
        """D4 map of lifted unit lattice steps from child to parent."""

        mapping = self.direction_map()
        ex = direction_vector(mapping["+x"])
        ey = direction_vector(mapping["+y"])
        # Columns are the images of child e_x and e_y.
        return ((ex[0], ey[0]), (ex[1], ey[1]))

    def homology_matrix(self) -> Matrix2:
        """Map child torus winding coordinates into parent winding coordinates."""

        a, b = self.parent.a, self.parent.b
        adj_parent: Matrix2 = ((a, b), (-b, a))
        numerator = matmul(
            matmul(adj_parent, self.direction_matrix()),
            self.child.period_matrix,
        )
        n = self.parent.n
        if any(value % n for row in numerator for value in row):
            raise AssertionError("cover does not induce an integral homology map")
        return (
            (numerator[0][0] // n, numerator[0][1] // n),
            (numerator[1][0] // n, numerator[1][1] // n),
        )

    def map_winding(self, child_winding: Vector2) -> Vector2:
        return matvec(self.homology_matrix(), child_winding)

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

    def verify_homology(self) -> None:
        direction = self.direction_matrix()
        if abs(det(direction)) != 1:
            raise AssertionError("local direction map is not D4/unimodular")
        homology = self.homology_matrix()
        if abs(det(homology)) != self.degree:
            raise AssertionError("homology map determinant does not equal cover degree")

        # Check the two child period generators directly in lifted coordinates.
        child_periods = (
            (self.child.a, self.child.b),
            (-self.child.b, self.child.a),
        )
        a, b = self.parent.a, self.parent.b
        for column, displacement in enumerate(child_periods):
            mapped = matvec(direction, displacement)
            numerators = (
                a * mapped[0] + b * mapped[1],
                -b * mapped[0] + a * mapped[1],
            )
            if numerators[0] % self.parent.n or numerators[1] % self.parent.n:
                raise AssertionError("mapped child period is not a parent period")
            winding = (
                numerators[0] // self.parent.n,
                numerators[1] // self.parent.n,
            )
            expected = (homology[0][column], homology[1][column])
            if winding != expected:
                raise AssertionError("homology column disagrees with lifted period")


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
        try:
            CoveringMap(parent, child, t)
        except (AssertionError, ValueError):
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
        assert len(candidates) == 4  # four D4 local direction identifications
        cover = canonical_cover(parent, child)
        assert cover.degree == degree
        assert set(cover.direction_map()) == {"+x", "+y"}
        assert len(set(cover.direction_map().values())) == 2
        assert abs(det(cover.homology_matrix())) == degree
        cover.verify_partition()
        cover.verify_edges()
        cover.verify_homology()
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
        print(
            f"t={t} direction_map={cover.direction_map()} "
            f"homology={cover.homology_matrix()}"
        )
    cover = canonical_cover(parent, child)
    cover.verify_partition()
    cover.verify_edges()
    cover.verify_homology()
    print(f"canonical t={cover.t}; exact partition/edge/homology checks PASS")
    return 0


if __name__ == "__main__":
    _self_test()
    raise SystemExit(main())
