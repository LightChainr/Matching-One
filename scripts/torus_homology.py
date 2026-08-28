#!/usr/bin/env python3
"""Exact winding-subgroup classification for periodic graph configurations.

The geometry reference code stores every edge displacement in coordinates of
the universal cover.  Closing an edge inside one union-find component then
produces an integer linear combination of the two quotient period vectors.
This module keeps up to two independent such winding vectors per component,
which is sufficient to distinguish no wrapping, one-dimensional (including
spiral) wrapping, and cross wrapping on a two-dimensional torus.

Period lattices are represented by a 2x2 integer matrix P whose columns are
the chosen generators.  Closed-cycle displacements ``d`` convert to winding
coefficients by the exact identity

    P w = d    <=>    w = adj(P) d / det(P)

with integer arithmetic only: never a floating inverse.  Scalar axis/diamond
periods ``(px, py)`` remain valid input and are the diagonal special case
``P = diag(px, py)``.

``direction_0`` and ``direction_1`` refer to the two chosen generators.  For
the axis geometry these are the horizontal and vertical periods.  For the
diamond geometry they are the two diagonal periods; using generator-relative
names avoids silently identifying a rotated quotient direction with a lattice
axis.  Rank, ``either``, and ``cross`` are invariant under unimodular changes
of the period basis; the directional flags are not, because they are defined
relative to the declared generators.
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
# Columns are the two period generators: ((P00, P01), (P10, P11)).
PeriodMatrix = tuple[tuple[int, int], tuple[int, int]]
PeriodLike = Union[tuple[int, int], PeriodMatrix]


def determinant(period_matrix: PeriodMatrix) -> int:
    """Return det(P) for a 2x2 integer period matrix."""

    p00, p01 = period_matrix[0]
    p10, p11 = period_matrix[1]
    return p00 * p11 - p01 * p10


def adjugate(period_matrix: PeriodMatrix) -> PeriodMatrix:
    """Return the integer adjugate of a 2x2 matrix, adj(P) P = det(P) I."""

    p00, p01 = period_matrix[0]
    p10, p11 = period_matrix[1]
    return ((p11, -p01), (-p10, p00))


def matmul2(left: PeriodMatrix, right: PeriodMatrix) -> PeriodMatrix:
    """Return the 2x2 integer product ``left @ right``."""

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


def apply_matrix(matrix: PeriodMatrix, vector: Winding) -> Winding:
    """Return ``matrix @ vector`` with integer arithmetic."""

    x, y = vector
    return (
        matrix[0][0] * x + matrix[0][1] * y,
        matrix[1][0] * x + matrix[1][1] * y,
    )


def unimodular_inverse(matrix: PeriodMatrix) -> PeriodMatrix:
    """Return the inverse of a det +/-1 integer matrix."""

    det = determinant(matrix)
    if abs(det) != 1:
        raise ValueError("matrix is not unimodular (require det +/-1)")
    adj = adjugate(matrix)
    if det == 1:
        return adj
    return ((-adj[0][0], -adj[0][1]), (-adj[1][0], -adj[1][1]))


def coerce_period_matrix(periods: PeriodLike) -> PeriodMatrix:
    """Accept diagonal ``(px, py)`` or a 2x2 generator matrix with det != 0."""

    if not isinstance(periods, (tuple, list)) or len(periods) != 2:
        raise ValueError("periods must be a 2-vector or a 2x2 integer matrix")
    first = periods[0]
    if isinstance(first, int):
        px, py = int(periods[0]), int(periods[1])
        if px <= 0 or py <= 0:
            raise ValueError("periods must be positive")
        return ((px, 0), (0, py))
    if isinstance(first, (tuple, list)) and len(first) == 2:
        try:
            matrix: PeriodMatrix = (
                (int(periods[0][0]), int(periods[0][1])),
                (int(periods[1][0]), int(periods[1][1])),
            )
        except (TypeError, ValueError, IndexError) as exc:
            raise ValueError("period matrix entries must be integers") from exc
        if determinant(matrix) == 0:
            raise ValueError("period matrix must have nonzero integer determinant")
        return matrix
    raise ValueError("periods must be a 2-vector or a 2x2 integer matrix")


def winding_coefficients(dx: int, dy: int, periods: PeriodLike) -> Winding:
    """Convert a cover displacement into generator windings via adj(P)/det(P).

    Raises ValueError when ``(dx, dy)`` is not in the period lattice.
    """

    matrix = coerce_period_matrix(periods)
    det = determinant(matrix)
    numerator = apply_matrix(adjugate(matrix), (dx, dy))
    if numerator[0] % det or numerator[1] % det:
        raise ValueError(
            "closed-cycle displacement is not in the quotient period "
            f"lattice: ({dx}, {dy}) versus period matrix {matrix}"
        )
    return numerator[0] // det, numerator[1] // det


def reduce_by_period_matrix(x: int, y: int, periods: PeriodLike) -> Winding:
    """Reduce a cover point into the half-open parallelogram of ``P``."""

    matrix = coerce_period_matrix(periods)
    det = determinant(matrix)
    numerator = apply_matrix(adjugate(matrix), (x, y))
    winding = (numerator[0] // det, numerator[1] // det)
    shift = apply_matrix(matrix, winding)
    return x - shift[0], y - shift[1]


def fundamental_domain_sites(periods: PeriodLike) -> tuple[Winding, ...]:
    """Return the unique Z^2 representatives in the half-open parallelogram.

    There are exactly ``|det(P)|`` such points.  Order is ``(y, x)`` so the
    axis case ``P = diag(L, L)`` matches ``axis_geometry(L)``.
    """

    matrix = coerce_period_matrix(periods)
    det = determinant(matrix)
    p00, p01 = matrix[0]
    p10, p11 = matrix[1]
    xs = (0, p00, p01, p00 + p01)
    ys = (0, p10, p11, p10 + p11)
    sites: list[Winding] = []
    for x in range(min(xs), max(xs) + 1):
        for y in range(min(ys), max(ys) + 1):
            reduced = reduce_by_period_matrix(x, y, matrix)
            if reduced == (x, y):
                sites.append((x, y))
    if len(sites) != abs(det):
        raise RuntimeError(
            f"fundamental domain of {matrix} has {len(sites)} sites, "
            f"expected |det|={abs(det)}"
        )
    return tuple(sorted(sites, key=lambda site: (site[1], site[0])))


def gaussian_period_matrix(a: int, b: int) -> PeriodMatrix:
    """Period matrix with columns ``(a, b)`` and ``(-b, a)``."""

    if a == 0 and b == 0:
        raise ValueError("Gaussian generators must be nonzero")
    return ((a, -b), (b, a))


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


def qspan_contains(vector: Winding, basis: Sequence[Winding]) -> bool:
    """Return whether ``vector`` lies in the rational span of ``basis``."""

    if vector == (0, 0):
        return True
    if not basis:
        return False
    if len(basis) == 1:
        x0, y0 = basis[0]
        x1, y1 = vector
        return x0 * y1 == y0 * x1
    x0, y0 = basis[0]
    x1, y1 = basis[1]
    if x0 * y1 != y0 * x1:
        return True
    return qspan_contains(vector, basis[:1])


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

    def __init__(self, n: int, periods: PeriodLike) -> None:
        if n < 0:
            raise ValueError("n must be nonnegative")
        self.period_matrix = coerce_period_matrix(periods)
        self.determinant = determinant(self.period_matrix)
        self.adjugate = adjugate(self.period_matrix)
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
        numerator = apply_matrix(self.adjugate, (dx, dy))
        det = self.determinant
        if numerator[0] % det or numerator[1] % det:
            raise ValueError(
                "closed-cycle displacement is not in the quotient period "
                f"lattice: ({dx}, {dy}) versus period matrix "
                f"{self.period_matrix}"
            )
        return numerator[0] // det, numerator[1] // det

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
    periods: PeriodLike,
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


def geometry_period_matrix(geometry: object) -> PeriodMatrix:
    """Return the 2x2 generator matrix used by a reference geometry."""

    matrix = getattr(geometry, "period_matrix", None)
    if matrix is not None:
        return coerce_period_matrix(matrix)
    return coerce_period_matrix(geometry_periods(geometry))


def geometry_periods(geometry: object) -> tuple[int, int]:
    """Return universal-cover periods for the repository's reference geometries.

    Axis and diamond embeddings are diagonal in their stored coordinates, so
    this legacy helper still returns a positive pair.  General 2x2 quotients
    should use :func:`geometry_period_matrix`.
    """

    matrix = getattr(geometry, "period_matrix", None)
    if matrix is not None:
        coerced = coerce_period_matrix(matrix)
        if coerced[0][1] == 0 and coerced[1][0] == 0:
            px, py = coerced[0][0], coerced[1][1]
            if px > 0 and py > 0:
                return px, py
        raise ValueError(
            "geometry_periods only supports diagonal embeddings; "
            "use geometry_period_matrix for a general 2x2 period matrix"
        )
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
    components = component_homologies(
        active, edges, geometry_period_matrix(geometry)
    )
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


def exhaustive_channel_counts(geometry: object) -> dict[str, int]:
    """Return rank and directional occupancy counts over all 2^N masks."""

    n = int(getattr(geometry, "n"))
    if n > 24:
        raise ValueError(
            f"exhaustive homology enumeration would require 2^{n} configurations"
        )
    counts = {"rank0": 0, "rank1": 0, "rank2": 0, "d0": 0, "d1": 0}
    for mask in range(1 << n):
        active = [bool((mask >> vertex) & 1) for vertex in range(n)]
        channels, _ = classify_configuration(geometry, active)
        counts[f"rank{channels.max_rank}"] += 1
        counts["d0"] += int(channels.direction_0)
        counts["d1"] += int(channels.direction_1)
    return counts


def transport_basis(basis: Sequence[Winding], matrix: PeriodMatrix) -> tuple[Winding, ...]:
    """Apply an integer 2x2 map to a primitive winding basis and recanonicalize."""

    transported: list[Winding] = []
    for vector in basis:
        image = apply_matrix(matrix, vector)
        _extend_basis(transported, image)
    return tuple(transported)
