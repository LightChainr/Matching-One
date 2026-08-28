#!/usr/bin/env python3
"""Reference topology engine for arbitrary integer-period square tori.

The quotient is ``Z^2 / P Z^2``.  ``P`` is written in ordinary row-major
form, while its *columns* are the two period vectors.  Closed lifted
displacements are converted to period-basis windings exactly as

    winding = adj(P) * displacement / det(P).

No floating-point inverse or axis-alignment assumption is used.  The module
is a correctness reference for tiny exhaustive checks and for validating
faster Monte Carlo implementations.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from math import gcd
from typing import Dict, Iterable, Sequence, Tuple

from matched_torus_reference import Edge
from torus_homology import WrappingChannels, wrapping_channels


Vector = Tuple[int, int]
Matrix = Tuple[Tuple[int, int], Tuple[int, int]]


def _matrix_entries(matrix: Matrix) -> Tuple[int, int, int, int]:
    try:
        (a, b), (c, d) = matrix
    except (TypeError, ValueError) as exc:
        raise ValueError("matrix must be 2x2") from exc
    if not all(isinstance(value, int) for value in (a, b, c, d)):
        raise TypeError("matrix entries must be integers")
    return a, b, c, d


def determinant(matrix: Matrix) -> int:
    a, b, c, d = _matrix_entries(matrix)
    return a * d - b * c


def matrix_vector(matrix: Matrix, vector: Vector) -> Vector:
    a, b, c, d = _matrix_entries(matrix)
    x, y = vector
    return a * x + b * y, c * x + d * y


def matrix_product(left: Matrix, right: Matrix) -> Matrix:
    a, b, c, d = _matrix_entries(left)
    e, f, g, h = _matrix_entries(right)
    return ((a * e + b * g, a * f + b * h),
            (c * e + d * g, c * f + d * h))


def unimodular_inverse(matrix: Matrix) -> Matrix:
    a, b, c, d = _matrix_entries(matrix)
    det = determinant(matrix)
    if abs(det) != 1:
        raise ValueError("basis change must be unimodular (determinant +/-1)")
    return ((d // det, -b // det), (-c // det, a // det))


def _canonical_sign(vector: Vector) -> Vector:
    x, y = vector
    if x < 0 or (x == 0 and y < 0):
        return -x, -y
    return vector


def _rank(vectors: Iterable[Vector]) -> int:
    nonzero = [vector for vector in vectors if vector != (0, 0)]
    if not nonzero:
        return 0
    x0, y0 = nonzero[0]
    if any(x0 * y != y0 * x for x, y in nonzero[1:]):
        return 2
    return 1


@dataclass(frozen=True)
class IntegerPeriods:
    """A nonsingular integer period matrix with exact quotient operations."""

    matrix: Matrix

    def __post_init__(self) -> None:
        _matrix_entries(self.matrix)
        if self.det == 0:
            raise ValueError("period matrix must be nonsingular")

    @property
    def det(self) -> int:
        return determinant(self.matrix)

    @property
    def order(self) -> int:
        return abs(self.det)

    @property
    def adjugate(self) -> Matrix:
        a, b, c, d = _matrix_entries(self.matrix)
        return ((d, -b), (-c, a))

    def period_vector(self, winding: Vector) -> Vector:
        return matrix_vector(self.matrix, winding)

    def winding(self, displacement: Vector) -> Vector:
        """Return exact period-basis coordinates of a closed displacement."""

        numerator = matrix_vector(self.adjugate, displacement)
        if numerator[0] % self.det or numerator[1] % self.det:
            raise ValueError(
                "closed-cycle displacement is not in the quotient period "
                f"lattice: {displacement} versus {self.matrix}"
            )
        return numerator[0] // self.det, numerator[1] // self.det

    def quotient_key(self, point: Vector) -> Vector:
        """Return an exact, hashable key for a coset of the period lattice."""

        numerator = matrix_vector(self.adjugate, point)
        modulus = self.order
        return numerator[0] % modulus, numerator[1] % modulus

    def equivalent(self, first: Vector, second: Vector) -> bool:
        return self.quotient_key(first) == self.quotient_key(second)

    def change_basis(self, change: Matrix) -> "IntegerPeriods":
        """Return the same lattice with period basis ``P * change``."""

        unimodular_inverse(change)  # validate before constructing the result
        return IntegerPeriods(matrix_product(self.matrix, change))


@dataclass(frozen=True)
class IntegerTorusGeometry:
    """Finite square-lattice graph on an arbitrary integer-period quotient."""

    name: str
    periods: IntegerPeriods
    coordinates: Tuple[Vector, ...]
    primal_edges: Tuple[Edge, ...]
    matching_edges: Tuple[Edge, ...]
    L: int = 0
    physical_period: str = "integer-matrix"
    _vertex_by_key: Dict[Vector, int] = field(
        default_factory=dict, repr=False, compare=False
    )

    @property
    def n(self) -> int:
        return len(self.coordinates)

    def vertex(self, point: Vector) -> int:
        return self._vertex_by_key[self.periods.quotient_key(point)]


def _quotient_representatives(periods: IntegerPeriods) -> Tuple[Vector, ...]:
    """Enumerate the finite quotient by a deterministic generator BFS."""

    representatives: Dict[Vector, Vector] = {periods.quotient_key((0, 0)): (0, 0)}
    queue = deque([(0, 0)])
    generators = ((1, 0), (0, 1))
    while queue:
        x, y = queue.popleft()
        for dx, dy in generators:
            point = x + dx, y + dy
            key = periods.quotient_key(point)
            if key not in representatives:
                representatives[key] = point
                queue.append(point)
    if len(representatives) != periods.order:
        raise AssertionError(
            f"quotient enumeration found {len(representatives)} of "
            f"{periods.order} cosets"
        )
    return tuple(representatives[key] for key in sorted(representatives))


def _make_edges(
    coordinates: Sequence[Vector],
    vertex_by_key: Dict[Vector, int],
    periods: IntegerPeriods,
    vectors: Iterable[Vector],
) -> Tuple[Edge, ...]:
    edges = []
    for source, (x, y) in enumerate(coordinates):
        for dx, dy in vectors:
            target_key = periods.quotient_key((x + dx, y + dy))
            edges.append(Edge(source, vertex_by_key[target_key], dx, dy))
    return tuple(edges)


def integer_torus_geometry(
    matrix: Matrix,
    *,
    name: str = "integer-period",
    L: int = 0,
    primal_vectors: Sequence[Vector] = ((1, 0), (0, 1)),
    matching_vectors: Sequence[Vector] = (
        (1, 0), (0, 1), (1, 1), (-1, 1)
    ),
) -> IntegerTorusGeometry:
    periods = IntegerPeriods(matrix)
    coordinates = _quotient_representatives(periods)
    vertex_by_key = {
        periods.quotient_key(point): vertex
        for vertex, point in enumerate(coordinates)
    }
    return IntegerTorusGeometry(
        name=name,
        periods=periods,
        coordinates=coordinates,
        primal_edges=_make_edges(
            coordinates, vertex_by_key, periods, primal_vectors
        ),
        matching_edges=_make_edges(
            coordinates, vertex_by_key, periods, matching_vectors
        ),
        L=L,
        physical_period=str(matrix),
        _vertex_by_key=vertex_by_key,
    )


def axis_integer_torus(L: int) -> IntegerTorusGeometry:
    if L <= 0:
        raise ValueError("L must be positive")
    return integer_torus_geometry(((L, 0), (0, L)), name="axis", L=L)


def diamond_integer_torus(L: int) -> IntegerTorusGeometry:
    if L <= 0:
        raise ValueError("L must be positive")
    # Columns are (L,L) and (-L,L), in original square-lattice coordinates.
    return integer_torus_geometry(((L, -L), (L, L)), name="diamond", L=L)


def gaussian_integer_torus(a: int, b: int) -> IntegerTorusGeometry:
    if a <= 0 or b < 0 or gcd(a, b) != 1:
        raise ValueError("require a>0, b>=0, and gcd(a,b)=1")
    return integer_torus_geometry(
        ((a, -b), (b, a)), name=f"gaussian-{a}-{b}"
    )


@dataclass(frozen=True)
class IntegerComponentHomology:
    """Exact cycle generators and rational homology rank of one component."""

    root: int
    size: int
    generators: Tuple[Vector, ...]

    @property
    def rank(self) -> int:
        return _rank(self.generators)

    @property
    def basis(self) -> Tuple[Vector, ...]:
        basis = []
        for vector in self.generators:
            if not basis or _rank(basis + [vector]) > len(basis):
                divisor = gcd(abs(vector[0]), abs(vector[1]))
                primitive = vector[0] // divisor, vector[1] // divisor
                basis.append(_canonical_sign(primitive))
            if len(basis) == 2:
                break
        return tuple(basis)

    @property
    def direction_0(self) -> bool:
        return any(x != 0 for x, _ in self.generators)

    @property
    def direction_1(self) -> bool:
        return any(y != 0 for _, y in self.generators)

    @property
    def either(self) -> bool:
        return self.rank > 0

    @property
    def both(self) -> bool:
        return self.direction_0 and self.direction_1

    @property
    def cross(self) -> bool:
        return self.rank == 2


class IntegerHomologyUnionFind:
    """Potential union-find using exact general-matrix winding conversion."""

    def __init__(self, n: int, periods: IntegerPeriods) -> None:
        if n < 0:
            raise ValueError("n must be nonnegative")
        self.periods = periods
        self.parent = list(range(n))
        self.size = [1] * n
        self.delta_x = [0] * n
        self.delta_y = [0] * n
        self.generators = [set() for _ in range(n)]

    def find(self, vertex: int) -> Tuple[int, int, int]:
        if self.parent[vertex] == vertex:
            return vertex, 0, 0
        parent = self.parent[vertex]
        root, parent_x, parent_y = self.find(parent)
        dx = self.delta_x[vertex] + parent_x
        dy = self.delta_y[vertex] + parent_y
        self.parent[vertex] = root
        self.delta_x[vertex] = dx
        self.delta_y[vertex] = dy
        return root, dx, dy

    def add_edge(self, i: int, j: int, edge_dx: int, edge_dy: int) -> None:
        root_i, ix, iy = self.find(i)
        root_j, jx, jy = self.find(j)
        root_dx = ix + edge_dx - jx
        root_dy = iy + edge_dy - jy

        if root_i == root_j:
            winding = self.periods.winding((root_dx, root_dy))
            if winding != (0, 0):
                self.generators[root_i].add(_canonical_sign(winding))
            return

        if self.size[root_i] < self.size[root_j]:
            root_i, root_j = root_j, root_i
            root_dx, root_dy = -root_dx, -root_dy
        self.parent[root_j] = root_i
        self.delta_x[root_j] = root_dx
        self.delta_y[root_j] = root_dy
        self.size[root_i] += self.size[root_j]
        self.generators[root_i].update(self.generators[root_j])
        self.generators[root_j].clear()

    def component(self, vertex: int) -> IntegerComponentHomology:
        root, _, _ = self.find(vertex)
        return IntegerComponentHomology(
            root, self.size[root], tuple(sorted(self.generators[root]))
        )


def component_homologies(
    active: Sequence[bool],
    edges: Iterable[Edge],
    periods: IntegerPeriods,
) -> Tuple[IntegerComponentHomology, ...]:
    union_find = IntegerHomologyUnionFind(len(active), periods)
    for edge in edges:
        if active[edge.i] and active[edge.j]:
            union_find.add_edge(edge.i, edge.j, edge.dx, edge.dy)
    roots = set()
    for vertex, enabled in enumerate(active):
        if enabled:
            root, _, _ = union_find.find(vertex)
            roots.add(root)
    return tuple(
        IntegerComponentHomology(
            root,
            union_find.size[root],
            tuple(sorted(union_find.generators[root])),
        )
        for root in sorted(roots)
    )


def classify_configuration(
    geometry: IntegerTorusGeometry,
    active: Sequence[bool],
    *,
    matching: bool = False,
) -> Tuple[WrappingChannels, Tuple[IntegerComponentHomology, ...]]:
    if len(active) != geometry.n:
        raise ValueError("active mask length does not match geometry")
    edges = geometry.matching_edges if matching else geometry.primal_edges
    components = component_homologies(active, edges, geometry.periods)
    return wrapping_channels(components), components
