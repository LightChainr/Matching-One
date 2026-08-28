#!/usr/bin/env python3
"""Exact winding-subgroup classification for periodic graph configurations.

The geometry reference code stores every edge displacement in coordinates of
the universal cover.  Closing an edge inside one union-find component then
produces an integer multiple of the two quotient periods.  This module keeps
up to two independent such winding vectors per component, which is sufficient
to distinguish no wrapping, one-dimensional (including spiral) wrapping, and
cross wrapping on a two-dimensional torus.

``direction_0`` and ``direction_1`` refer to the two quotient generators.  For
the axis geometry these are the horizontal and vertical periods.  For the
diamond geometry they are the two diagonal periods; using generator-relative
names avoids silently identifying a rotated quotient direction with a lattice
axis.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from typing import Iterable, Protocol, Sequence, Union


class DisplacedEdge(Protocol):
    i: int
    j: int
    dx: int
    dy: int


Winding = tuple[int, int]


def _primitive(vector: Winding) -> Winding:
    """Return a canonical primitive representative of a nonzero direction."""

    x, y = vector
    divisor = gcd(abs(x), abs(y))
    if divisor == 0:
        raise ValueError("the zero vector has no winding direction")
    x //= divisor
    y //= divisor
    if x < 0 or (x == 0 and y < 0):
        x, y = -x, -y
    return x, y


def _extend_basis(basis: list[Winding], vector: Winding) -> None:
    """Add ``vector`` if it increases the rational rank of a 2-D basis."""

    if vector == (0, 0) or len(basis) == 2:
        return
    vector = _primitive(vector)
    if not basis:
        basis.append(vector)
        return
    x0, y0 = basis[0]
    x1, y1 = vector
    if x0 * y1 != y0 * x1:
        basis.append(vector)


@dataclass(frozen=True)
class ComponentHomology:
    """Winding subgroup summary for one occupied connected component."""

    root: int
    size: int
    basis: tuple[Winding, ...]

    @property
    def rank(self) -> int:
        return len(self.basis)

    @property
    def direction_0(self) -> bool:
        return any(x != 0 for x, _ in self.basis)

    @property
    def direction_1(self) -> bool:
        return any(y != 0 for _, y in self.basis)

    @property
    def either(self) -> bool:
        return self.rank > 0

    @property
    def both(self) -> bool:
        return self.direction_0 and self.direction_1

    @property
    def cross(self) -> bool:
        return self.rank == 2


@dataclass(frozen=True)
class WrappingChannels:
    """Configuration-level wrapping events aggregated over its components."""

    max_rank: int
    direction_0: bool
    direction_1: bool
    either: bool
    both: bool
    cross: bool

    @property
    def horizontal(self) -> bool:
        """Alias for generator 0, meaningful as horizontal on axis tori."""

        return self.direction_0

    @property
    def vertical(self) -> bool:
        """Alias for generator 1, meaningful as vertical on axis tori."""

        return self.direction_1

    def as_dict(self) -> dict[str, Union[bool, int]]:
        return {
            "rank": self.max_rank,
            "direction_0": self.direction_0,
            "direction_1": self.direction_1,
            "either": self.either,
            "both": self.both,
            "cross": self.cross,
        }


class HomologyUnionFind:
    """Displacement-potential union-find retaining component winding bases."""

    def __init__(self, n: int, periods: tuple[int, int]) -> None:
        if n < 0:
            raise ValueError("n must be nonnegative")
        if periods[0] <= 0 or periods[1] <= 0:
            raise ValueError("periods must be positive")
        self.periods = periods
        self.parent = list(range(n))
        self.size = [1] * n
        self.delta_x = [0] * n
        self.delta_y = [0] * n
        self.basis: list[list[Winding]] = [[] for _ in range(n)]

    def find(self, x: int) -> tuple[int, int, int]:
        if self.parent[x] == x:
            return x, 0, 0
        parent = self.parent[x]
        root, parent_x, parent_y = self.find(parent)
        dx = self.delta_x[x] + parent_x
        dy = self.delta_y[x] + parent_y
        self.parent[x] = root
        self.delta_x[x] = dx
        self.delta_y[x] = dy
        return root, dx, dy

    def _winding(self, dx: int, dy: int) -> Winding:
        px, py = self.periods
        if dx % px or dy % py:
            raise ValueError(
                "closed-cycle displacement is not in the quotient period "
                f"lattice: ({dx}, {dy}) versus periods ({px}, {py})"
            )
        return dx // px, dy // py

    def add_edge(self, i: int, j: int, edge_dx: int, edge_dy: int) -> None:
        """Add an edge satisfying ``pos(j)=pos(i)+(edge_dx, edge_dy)``."""

        root_i, ix, iy = self.find(i)
        root_j, jx, jy = self.find(j)
        root_dx = ix + edge_dx - jx
        root_dy = iy + edge_dy - jy

        if root_i == root_j:
            _extend_basis(self.basis[root_i], self._winding(root_dx, root_dy))
            return

        if self.size[root_i] < self.size[root_j]:
            root_i, root_j = root_j, root_i
            root_dx, root_dy = -root_dx, -root_dy

        self.parent[root_j] = root_i
        self.delta_x[root_j] = root_dx
        self.delta_y[root_j] = root_dy
        self.size[root_i] += self.size[root_j]
        for vector in self.basis[root_j]:
            _extend_basis(self.basis[root_i], vector)
        self.basis[root_j] = []

    def component(self, vertex: int) -> ComponentHomology:
        root, _, _ = self.find(vertex)
        return ComponentHomology(root, self.size[root], tuple(self.basis[root]))


def component_homologies(
    active: Sequence[bool],
    edges: Iterable[DisplacedEdge],
    periods: tuple[int, int],
) -> tuple[ComponentHomology, ...]:
    """Return homology summaries for every active connected component."""

    union_find = HomologyUnionFind(len(active), periods)
    for edge in edges:
        if active[edge.i] and active[edge.j]:
            union_find.add_edge(edge.i, edge.j, edge.dx, edge.dy)

    roots: set[int] = set()
    for vertex, enabled in enumerate(active):
        if enabled:
            root, _, _ = union_find.find(vertex)
            roots.add(root)
    return tuple(
        ComponentHomology(root, union_find.size[root], tuple(union_find.basis[root]))
        for root in sorted(roots)
    )


def wrapping_channels(components: Iterable[ComponentHomology]) -> WrappingChannels:
    """Aggregate component subgroups into directional/either/both/cross events."""

    components = tuple(components)
    direction_0 = any(component.direction_0 for component in components)
    direction_1 = any(component.direction_1 for component in components)
    max_rank = max((component.rank for component in components), default=0)
    return WrappingChannels(
        max_rank=max_rank,
        direction_0=direction_0,
        direction_1=direction_1,
        either=direction_0 or direction_1,
        both=direction_0 and direction_1,
        cross=any(component.cross for component in components),
    )


def geometry_periods(geometry: object) -> tuple[int, int]:
    """Return universal-cover periods for the repository's reference geometries."""

    name = getattr(geometry, "name", None)
    length = getattr(geometry, "L", None)
    if not isinstance(length, int) or length <= 0:
        raise ValueError("geometry must expose a positive integer L")
    if name == "axis":
        return length, length
    if name == "diamond":
        return 2 * length, 2 * length
    raise ValueError(f"unsupported geometry {name!r}")


def classify_configuration(
    geometry: object,
    active: Sequence[bool],
    *,
    matching: bool = False,
) -> tuple[WrappingChannels, tuple[ComponentHomology, ...]]:
    """Classify occupied sites using primal or matching-lattice connectivity."""

    coordinates = getattr(geometry, "coordinates")
    if len(active) != len(coordinates):
        raise ValueError("active mask length does not match geometry")
    edges = getattr(geometry, "matching_edges" if matching else "primal_edges")
    components = component_homologies(active, edges, geometry_periods(geometry))
    return wrapping_channels(components), components


def matching_channel_differences(
    geometry: object, active: Sequence[bool]
) -> dict[str, int]:
    """Return coupled black-primal minus white-matching channel indicators."""

    black, _ = classify_configuration(geometry, active)
    white, _ = classify_configuration(
        geometry, [not value for value in active], matching=True
    )
    return {
        channel: int(getattr(black, channel)) - int(getattr(white, channel))
        for channel in (
            "cross",
            "both",
            "either",
            "direction_0",
            "direction_1",
        )
    }
