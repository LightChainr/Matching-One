#!/usr/bin/env python3
"""Finite controls for a staircase of NN crossings around an integer period.

The probability theorem is in docs/manuscripts/geometric-balance/manuscript.md.
This script does not estimate p_c or prove the imported RSW/sharpness inputs.
All geometry and graph calculations below use integers / Fraction.
"""
from __future__ import annotations

import argparse
import json
from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path
from typing import Iterable

Point = tuple[int, int]
NN: tuple[Point, ...] = ((1, 0), (-1, 0), (0, 1), (0, -1))


def det(a: Point, b: Point) -> int:
    return a[0] * b[1] - a[1] * b[0]


def add(a: Point, b: Point) -> Point:
    return a[0] + b[0], a[1] + b[1]


def sub(a: Point, b: Point) -> Point:
    return a[0] - b[0], a[1] - b[1]


@dataclass(frozen=True)
class Torus:
    u: Point
    v: Point

    def __post_init__(self) -> None:
        if det(self.u, self.v) <= 0:
            raise ValueError("The ordered integer period basis must have positive determinant")

    @property
    def n(self) -> int:
        return det(self.u, self.v)

    def reduce(self, z: Point) -> Point:
        # Half-open fundamental parallelogram; floor is correct also at negative z.
        i = det(z, self.v) // self.n
        j = det(self.u, z) // self.n
        return z[0] - i*self.u[0] - j*self.v[0], z[1] - i*self.u[1] - j*self.v[1]

    def vertices(self) -> tuple[Point, ...]:
        seen = {(0, 0)}
        queue = deque([(0, 0)])
        while queue:
            z = queue.popleft()
            for d in ((1, 0), (0, 1)):
                nxt = self.reduce(add(z, d))
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
        if len(seen) != self.n:
            raise AssertionError("Quotient cardinality mismatch")
        return tuple(sorted(seen))

    def period_coordinates(self, z: Point) -> Point:
        a, b = det(z, self.v), det(self.u, z)
        if a % self.n or b % self.n:
            raise AssertionError("Graph cycle displacement is not a period")
        return a // self.n, b // self.n


@dataclass(frozen=True)
class Rectangle:
    x0: int
    x1: int
    y0: int
    y1: int
    direction: str

    def __post_init__(self) -> None:
        if self.x0 >= self.x1 or self.y0 >= self.y1 or self.direction not in ("h", "v"):
            raise ValueError("Nondegenerate h/v crossing rectangle required")

    def points(self) -> tuple[Point, ...]:
        return tuple((x, y) for x in range(self.x0, self.x1+1)
                     for y in range(self.y0, self.y1+1))

    def translated(self, z: Point) -> Rectangle:
        return Rectangle(self.x0+z[0], self.x1+z[0], self.y0+z[1], self.y1+z[1], self.direction)


def staircase(u: Point, scale: int) -> tuple[Point, ...]:
    """Exact-endpoint integer staircase; zero steps omitted, signs allowed."""
    if scale < 1 or u == (0, 0):
        raise ValueError("Positive scale and nonzero period required")
    k = (max(abs(u[0]), abs(u[1])) + scale - 1) // scale
    centers: list[Point] = [(0, 0)]
    for j in range(k):
        x0, y0 = (j*u[0])//k, (j*u[1])//k
        x1, y1 = ((j+1)*u[0])//k, ((j+1)*u[1])//k
        for z in ((x1, y0), (x1, y1)):
            if z != centers[-1]:
                centers.append(z)
    assert centers[-1] == u
    assert all((a[0] == b[0]) ^ (a[1] == b[1]) for a, b in zip(centers, centers[1:]))
    assert all(max(abs(a[0]-b[0]), abs(a[1]-b[1])) <= scale
               for a, b in zip(centers, centers[1:]))
    return tuple(centers)


def corridor(u: Point, scale: int) -> tuple[Rectangle, ...]:
    """Two crossing directions in each hub; one connector per cyclic step."""
    z = staircase(u, scale)
    events: list[Rectangle] = []
    for x, y in z[:-1]:
        for d in ("h", "v"):
            events.append(Rectangle(x-scale, x+scale, y-scale, y+scale, d))
    for a, b in zip(z, z[1:]):
        if a[1] == b[1]:
            events.append(Rectangle(min(a[0], b[0])-scale, max(a[0], b[0])+scale,
                                    a[1]-scale, a[1]+scale, "h"))
        else:
            events.append(Rectangle(a[0]-scale, a[0]+scale,
                                    min(a[1], b[1])-scale, max(a[1], b[1])+scale, "v"))
    return tuple(events)


def support(events: Iterable[Rectangle]) -> set[Point]:
    return {z for r in events for z in r.points()}


def rectangle_crossing(torus: Torus, rect: Rectangle, occupied: set[Point]) -> bool:
    """Planar-lift BFS. Periodic occupancy, but no extra edges across box sides."""
    allowed = {z for z in rect.points() if torus.reduce(z) in occupied}
    starts = [z for z in allowed if (z[0] == rect.x0 if rect.direction == "h" else z[1] == rect.y0)]
    seen = set(starts)
    queue = deque(starts)
    while queue:
        z = queue.popleft()
        if z[0] == rect.x1 if rect.direction == "h" else z[1] == rect.y1:
            return True
        for d in NN:
            nxt = add(z, d)
            if nxt in allowed and nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return False


def compiled_crossing(torus: Torus, rect: Rectangle, ids: dict[Point, int]):
    """Compile box incidence, retaining a distinct planar vertex per box point."""
    pts = rect.points()
    pos = {z: i for i, z in enumerate(pts)}
    weights = [1 << ids[torus.reduce(z)] for z in pts]
    starts = [i for i, z in enumerate(pts)
              if (z[0] == rect.x0 if rect.direction == "h" else z[1] == rect.y0)]
    ends = {i for i, z in enumerate(pts)
            if (z[0] == rect.x1 if rect.direction == "h" else z[1] == rect.y1)}
    adjacency = [[pos[add(z, d)] for d in NN if add(z, d) in pos] for z in pts]

    def check(mask: int) -> bool:
        queue = [i for i in starts if mask & weights[i]]
        seen = set(queue)
        for i in queue:
            if i in ends:
                return True
            for j in adjacency[i]:
                if j not in seen and mask & weights[j]:
                    seen.add(j)
                    queue.append(j)
        return False
    return check


def winding_vectors(torus: Torus, occupied: set[Point]) -> list[Point]:
    """Independent spanning-forest detector on the physical quotient NN graph."""
    potential: dict[Point, Point] = {}
    gains: list[Point] = []
    for start in sorted(occupied):
        if start in potential:
            continue
        potential[start] = (0, 0)
        queue = deque([start])
        while queue:
            z = queue.popleft()
            for d in NN:
                nxt = torus.reduce(add(z, d))
                if nxt not in occupied:
                    continue
                proposed = add(potential[z], d)
                if nxt not in potential:
                    potential[nxt] = proposed
                    queue.append(nxt)
                else:
                    gain = sub(proposed, potential[nxt])
                    if gain != (0, 0):
                        gains.append(torus.period_coordinates(gain))
    return gains


def rank(torus: Torus, occupied: set[Point]) -> int:
    gains = winding_vectors(torus, occupied)
    if not gains:
        return 0
    return 2 if any(det(gains[0], z) for z in gains[1:]) else 1


def squared_distance_to_segment(z: Point, u: Point) -> Fraction:
    length2 = u[0]**2 + u[1]**2
    dot = z[0]*u[0] + z[1]*u[1]
    if dot < 0:
        return Fraction(z[0]**2+z[1]**2)
    if dot > length2:
        d = sub(z, u)
        return Fraction(d[0]**2+d[1]**2)
    return Fraction(det(u, z)**2, length2)


def circle_packing(modulus: int, pattern: set[int]) -> tuple[list[int], set[int]]:
    """Greedy disjoint translates, used only to check the finite-group lemma."""
    if modulus < 1 or not pattern:
        raise ValueError("Positive group order and nonempty pattern required")
    pattern = {x % modulus for x in pattern}
    differences = {(x-y) % modulus for x in pattern for y in pattern}
    available = set(range(modulus))
    centers: list[int] = []
    while available:
        c = min(available)
        centers.append(c)
        available.difference_update((c+d) % modulus for d in differences)
    return centers, differences


def exhaustive_case(u: Point, v: Point, scale: int = 1) -> dict:
    torus = Torus(u, v)
    if torus.n > 16:
        raise ValueError("Exhaustive control limited to 16 sites; no implicit large census")
    vertices = torus.vertices()
    ids = {z: i for i, z in enumerate(vertices)}
    events = corridor(u, scale)
    injective = all(len({torus.reduce(z) for z in r.points()}) == len(r.points()) for r in events)
    if not injective:
        raise ValueError("These test rectangles do not inject; choose a different control geometry")
    checks = [compiled_crossing(torus, r, ids) for r in events]
    ring_counts = [0]*(torus.n+1)
    failures: list[int] = []
    for mask in range(1 << torus.n):
        if all(f(mask) for f in checks):
            occ = {vertices[i] for i in range(torus.n) if mask >> i & 1}
            ring_counts[mask.bit_count()] += 1
            if rank(torus, occ) == 0:
                failures.append(mask)
    return {"u": u, "v": v, "scale": scale, "N": torus.n,
            "rectangles": len(events), "all_rectangles_inject": injective,
            "configurations": 1 << torus.n, "ring_configurations": sum(ring_counts),
            "ring_counts_by_occupation": ring_counts, "implication_failures": failures}


def execute_controls() -> dict:
    # Tiny controls verify the geometric gluing, not the asymptotic ell>=64s constants.
    small = [exhaustive_case((4, 0), (1, 3)),
             exhaustive_case((4, 1), (0, 4)),
             exhaustive_case((4, 2), (0, 4))]
    geometry = []
    for u, v, s in [((64,0),(17,83),1), ((63,16),(-33,127),1),
                    ((32,57),(-130,85),1), ((48,48),(-80,112),1),
                    ((128,33),(-50,263),2), ((-65,17),(-100,-171),1)]:
        # Require a reduced basis with shortest u; this criterion is exact in dimension two.
        uu = u[0]**2+u[1]**2
        vv = v[0]**2+v[1]**2
        dot = u[0]*v[0]+u[1]*v[1]
        if uu > vv or 2*abs(dot) > uu:
            q = (2*dot+uu)//(2*uu)
            v = v[0]-q*u[0], v[1]-q*u[1]
            vv = v[0]**2+v[1]**2
            dot = u[0]*v[0]+u[1]*v[1]
        assert uu <= vv and 2*abs(dot) <= uu
        torus = Torus(u,v)
        ev = corridor(u,s)
        pts = support(ev)
        maxdist = max(squared_distance_to_segment(z,u) for z in pts)
        assert maxdist <= 16*s*s
        assert uu >= (64*s)**2
        assert all(len({torus.reduce(z) for z in r.points()}) == len(r.points()) for r in ev)
        occ = {torus.reduce(z) for z in pts}
        assert all(rectangle_crossing(torus,r,occ) for r in ev)
        assert rank(torus,occ) > 0
        # No floating square roots enter these diameter/area inequalities.
        assert len(ev)**2 * s*s <= 64*uu  # event count <= 8 ell/s
        geometry.append({"u":u,"v":v,"scale":s,"ell_squared":uu,"N":torus.n,
                         "ambient_gcd_u":gcd(abs(u[0]),abs(u[1])),
                         "events":len(ev),"lifted_support_size":len(pts),
                         "max_distance_squared_to_segment":str(maxdist),
                         "rank_of_full_corridor_support":rank(torus,occ)})
    packing = []
    for n, pattern in [(23,{0,1,4}), (64,{0,2,5,11}), (101,{0,1,2,8,9})]:
        centers, diff = circle_packing(n,pattern)
        packed = [{(x+c)%n for x in pattern} for c in centers]
        assert all(not(a & b) for i,a in enumerate(packed) for b in packed[i+1:])
        assert len(centers)*len(diff) >= n
        packing.append({"group_order":n,"pattern":sorted(pattern),"difference_size":len(diff),
                        "centers":centers,"cover_inequality":len(centers)*len(diff)>=n})
    return {"schema":"matching-one.oblique-corridor.v1", "date":"2026-09-13",
            "scope":"deterministic NN staircase gluing and finite-group packing; not an RSW computation",
            "tiny_exhaustive":small,"geometry":geometry,"packing":packing,
            "total_exhaustive_configurations":sum(c["configurations"] for c in small),
            "all_implications_hold":not any(c["implication_failures"] for c in small)}


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        raise SystemExit("Refusing to replace an existing result; use a new --output path")
    result = execute_controls()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"configurations":result["total_exhaustive_configurations"],
                      "all_implications_hold":result["all_implications_hold"]}))


if __name__ == "__main__":
    main()
