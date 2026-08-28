#!/usr/bin/env python3
"""Reference checker for the finite square/matching-lattice identity.

This is intentionally a correctness-first implementation, not a production
Monte Carlo kernel. It supports two periodic quotients of Z^2:

axis:
    periods (L, 0), (0, L), N=L^2 sites.

diamond:
    periods (L, L), (-L, L), N=2 L^2 sites. The two periods have physical
    length sqrt(2) L and are rotated by pi/4 relative to the lattice axes.

For black sites use nearest-neighbour square connectivity. For the white
complement use the NN+NNN matching lattice. The program checks

    E[N_black - N_white_matching] - N * chi(p)
      = P(black wraps) - P(white matching wraps),

with chi(p)=p-2p^2+p^4.

The equality is the Mertens-Ziff finite matching relation. Exact mode
enumerates all 2^N site configurations and is therefore limited to tiny
systems. Monte Carlo mode is only a reference implementation for testing
faster C/C++/GPU kernels.
"""

from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import mpmath as mp


@dataclass(frozen=True)
class Edge:
    i: int
    j: int
    dx: int
    dy: int


@dataclass(frozen=True)
class Geometry:
    name: str
    L: int
    coordinates: tuple[tuple[int, int], ...]
    primal_edges: tuple[Edge, ...]
    matching_edges: tuple[Edge, ...]
    physical_period: str

    @property
    def n(self) -> int:
        return len(self.coordinates)


class WrapUnionFind:
    """Union-find with lattice displacement potentials and wrap detection.

    delta_x[x], delta_y[x] store position(x)-position(parent(x)) in the
    universal-cover coordinate system used to define the geometry. Closing a
    cycle with non-zero accumulated displacement therefore detects non-trivial
    torus homology.
    """

    def __init__(self, n: int) -> None:
        self.parent = list(range(n))
        self.size = [1] * n
        self.delta_x = [0] * n
        self.delta_y = [0] * n
        self.wrap = [False] * n

    def find(self, x: int) -> tuple[int, int, int]:
        if self.parent[x] == x:
            return x, 0, 0

        parent = self.parent[x]
        root, px, py = self.find(parent)
        dx = self.delta_x[x]
        dy = self.delta_y[x]
        self.parent[x] = root
        self.delta_x[x] = dx + px
        self.delta_y[x] = dy + py
        return root, self.delta_x[x], self.delta_y[x]

    def add_edge(self, i: int, j: int, edge_dx: int, edge_dy: int) -> None:
        """Add an edge satisfying pos(j)=pos(i)+(edge_dx, edge_dy)."""

        ri, ix, iy = self.find(i)
        rj, jx, jy = self.find(j)

        # Required position of root j relative to root i.
        root_dx = ix + edge_dx - jx
        root_dy = iy + edge_dy - jy

        if ri == rj:
            if root_dx != 0 or root_dy != 0:
                self.wrap[ri] = True
            return

        if self.size[ri] >= self.size[rj]:
            self.parent[rj] = ri
            self.delta_x[rj] = root_dx
            self.delta_y[rj] = root_dy
            self.size[ri] += self.size[rj]
            self.wrap[ri] = self.wrap[ri] or self.wrap[rj]
        else:
            self.parent[ri] = rj
            self.delta_x[ri] = -root_dx
            self.delta_y[ri] = -root_dy
            self.size[rj] += self.size[ri]
            self.wrap[rj] = self.wrap[ri] or self.wrap[rj]


def _make_edges(
    coordinates: list[tuple[int, int]],
    ids: dict[tuple[int, int], int],
    period: int,
    vectors: Iterable[tuple[int, int]],
) -> tuple[Edge, ...]:
    edges: list[Edge] = []
    for i, (x, y) in enumerate(coordinates):
        for dx, dy in vectors:
            target = ((x + dx) % period, (y + dy) % period)
            edges.append(Edge(i=i, j=ids[target], dx=dx, dy=dy))
    return tuple(edges)


def axis_geometry(L: int) -> Geometry:
    if L <= 0:
        raise ValueError("L must be positive")
    coordinates = [(x, y) for y in range(L) for x in range(L)]
    ids = {coordinate: i for i, coordinate in enumerate(coordinates)}

    primal = _make_edges(
        coordinates,
        ids,
        L,
        vectors=((1, 0), (0, 1)),
    )
    matching = _make_edges(
        coordinates,
        ids,
        L,
        vectors=((1, 0), (0, 1), (1, 1), (1, -1)),
    )
    return Geometry(
        name="axis",
        L=L,
        coordinates=tuple(coordinates),
        primal_edges=primal,
        matching_edges=matching,
        physical_period=f"{L}",
    )


def diamond_geometry(L: int) -> Geometry:
    if L <= 0:
        raise ValueError("L must be positive")

    # u=x+y, v=y-x. The periods (L,L) and (-L,L) become independent
    # periods 2L in u and v. Valid Z^2 sites have u and v of equal parity.
    period = 2 * L
    coordinates = [
        (u, v)
        for u in range(period)
        for v in range(period)
        if (u - v) % 2 == 0
    ]
    ids = {coordinate: i for i, coordinate in enumerate(coordinates)}

    # +x and +y in original coordinates become (1,-1) and (1,1).
    primal = _make_edges(
        coordinates,
        ids,
        period,
        vectors=((1, -1), (1, 1)),
    )

    # Add one orientation from each NNN pair: x+y -> (2,0),
    # -x+y -> (0,2).
    matching = _make_edges(
        coordinates,
        ids,
        period,
        vectors=((1, -1), (1, 1), (2, 0), (0, 2)),
    )
    return Geometry(
        name="diamond",
        L=L,
        coordinates=tuple(coordinates),
        primal_edges=primal,
        matching_edges=matching,
        physical_period=f"sqrt(2)*{L}",
    )


def cluster_stats(active: list[bool], edges: tuple[Edge, ...]) -> tuple[int, bool]:
    uf = WrapUnionFind(len(active))
    for edge in edges:
        if active[edge.i] and active[edge.j]:
            uf.add_edge(edge.i, edge.j, edge.dx, edge.dy)

    roots: set[int] = set()
    for i, enabled in enumerate(active):
        if enabled:
            root, _, _ = uf.find(i)
            roots.add(root)

    return len(roots), any(uf.wrap[root] for root in roots)


def matching_polynomial(p: mp.mpf) -> mp.mpf:
    return p - 2 * p**2 + p**4


def configuration_observables(
    geometry: Geometry, active: list[bool]
) -> tuple[int, int, bool, bool]:
    black_clusters, black_wrap = cluster_stats(active, geometry.primal_edges)
    white = [not value for value in active]
    white_clusters, white_wrap = cluster_stats(white, geometry.matching_edges)
    return black_clusters, white_clusters, black_wrap, white_wrap


def exact_check(geometry: Geometry, p: mp.mpf) -> dict[str, mp.mpf | int | str]:
    if geometry.n > 24:
        raise ValueError(
            f"exact enumeration would require 2^{geometry.n} configurations; "
            "use N<=24 or Monte Carlo mode"
        )

    q = 1 - p
    cluster_difference = mp.mpf(0)
    wrapping_difference = mp.mpf(0)

    for mask in range(1 << geometry.n):
        occupied = mask.bit_count()
        weight = p**occupied * q ** (geometry.n - occupied)
        active = [bool((mask >> i) & 1) for i in range(geometry.n)]
        nb, nw, wb, ww = configuration_observables(geometry, active)
        cluster_difference += weight * (nb - nw)
        wrapping_difference += weight * (int(wb) - int(ww))

    matching_value = cluster_difference - geometry.n * matching_polynomial(p)
    return {
        "geometry": geometry.name,
        "L": geometry.L,
        "N": geometry.n,
        "physical_period": geometry.physical_period,
        "p": p,
        "matching_function_cluster_side": matching_value,
        "matching_function_wrapping_side": wrapping_difference,
        "difference": matching_value - wrapping_difference,
    }


def monte_carlo_check(
    geometry: Geometry, p: float, samples: int, seed: int
) -> dict[str, float | int | str]:
    if not 0.0 <= p <= 1.0:
        raise ValueError("p must lie in [0,1]")
    if samples <= 0:
        raise ValueError("samples must be positive")

    rng = random.Random(seed)
    cluster_values: list[float] = []
    wrapping_values: list[float] = []
    chi = p - 2.0 * p * p + p**4

    for _ in range(samples):
        active = [rng.random() < p for _ in range(geometry.n)]
        nb, nw, wb, ww = configuration_observables(geometry, active)
        cluster_values.append((nb - nw) - geometry.n * chi)
        wrapping_values.append(float(int(wb) - int(ww)))

    def mean_se(values: list[float]) -> tuple[float, float]:
        mean = math.fsum(values) / len(values)
        if len(values) == 1:
            return mean, math.nan
        variance = math.fsum((value - mean) ** 2 for value in values) / (
            len(values) - 1
        )
        return mean, math.sqrt(variance / len(values))

    cluster_mean, cluster_se = mean_se(cluster_values)
    wrap_mean, wrap_se = mean_se(wrapping_values)
    paired = [a - b for a, b in zip(cluster_values, wrapping_values)]
    paired_mean, paired_se = mean_se(paired)

    return {
        "geometry": geometry.name,
        "L": geometry.L,
        "N": geometry.n,
        "physical_period": geometry.physical_period,
        "p": p,
        "samples": samples,
        "seed": seed,
        "matching_function_cluster_side": cluster_mean,
        "cluster_side_se": cluster_se,
        "matching_function_wrapping_side": wrap_mean,
        "wrapping_side_se": wrap_se,
        "paired_difference": paired_mean,
        "paired_difference_se": paired_se,
    }


def json_ready(value: object) -> object:
    if isinstance(value, mp.mpf):
        return mp.nstr(value, 50)
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", choices=("axis", "diamond"), required=True)
    parser.add_argument("--L", type=int, required=True)
    parser.add_argument("--p", required=True, help="site occupation probability")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--exact", action="store_true")
    mode.add_argument("--samples", type=int)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args()

    if args.dps < 40:
        raise SystemExit("--dps must be at least 40")
    mp.mp.dps = args.dps

    geometry = axis_geometry(args.L) if args.geometry == "axis" else diamond_geometry(args.L)

    try:
        if args.exact:
            p = mp.mpf(args.p)
            if not 0 <= p <= 1:
                raise ValueError("p must lie in [0,1]")
            result = exact_check(geometry, p)
        else:
            result = monte_carlo_check(geometry, float(args.p), args.samples, args.seed)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    for key, value in result.items():
        shown = mp.nstr(value, 30) if isinstance(value, mp.mpf) else value
        print(f"{key}: {shown}")

    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(
            json.dumps(
                {key: json_ready(value) for key, value in result.items()},
                indent=2,
            ),
            encoding="utf-8",
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
