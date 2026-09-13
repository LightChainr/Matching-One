#!/usr/bin/env python3
"""Finite geometry controls for the oblique RSW necklace construction.

Only integer arithmetic is used. This tests geometric supports, the periodic
seam, winding of explicit witnesses and disjoint translates. It does not
compute a near-critical probability, prove RSW, or independently prove the
all-size theorem. No old Matching-One implementation is imported.
"""
from __future__ import annotations

import argparse
from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from math import gcd, isqrt
import json
from pathlib import Path
import unittest

Point = tuple[int, int]


def det(a: Point, b: Point) -> int:
    return a[0] * b[1] - a[1] * b[0]


def dot(a: Point, b: Point) -> int:
    return a[0] * b[0] + a[1] * b[1]


def add(a: Point, b: Point) -> Point:
    return a[0] + b[0], a[1] + b[1]


def nearest(num: int, den: int) -> int:
    """Nearest integer, with ties upwards; commutes with integer translation."""
    if den <= 0:
        raise ValueError("denominator must be positive")
    return (2 * num + den) // (2 * den)


def ceil_scaled_sqrt(square: int, multiplier: int, divisor: int) -> int:
    """ceil(multiplier * sqrt(square) / divisor), exactly."""
    if square < 0 or multiplier <= 0 or divisor <= 0:
        raise ValueError("invalid square-root arguments")
    n = isqrt(multiplier * multiplier * square) // divisor
    return n + (n * n * divisor * divisor < multiplier * multiplier * square)


def reduced_basis(u: Point, v: Point) -> tuple[Point, Point]:
    """Two-dimensional Lagrange reduction; u is a shortest lattice vector."""
    if det(u, v) == 0:
        raise ValueError("periods must be independent")
    for _ in range(256):
        if dot(v, v) < dot(u, u):
            u, v = v, u
        q = nearest(dot(u, v), dot(u, u))
        if q == 0:
            break
        v = v[0] - q * u[0], v[1] - q * u[1]
    else:
        raise RuntimeError("reduction did not terminate")
    if u[0] < 0 or (u[0] == 0 and u[1] < 0):
        u = -u[0], -u[1]
    if det(u, v) < 0:
        v = -v[0], -v[1]
    return u, v


def bezout(a: int, b: int) -> tuple[int, int, int]:
    """Return g>=0, x, y with ax+by=g=gcd(a,b)."""
    r0, r1, x0, x1, y0, y1 = abs(a), abs(b), 1, 0, 0, 1
    while r1:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return r0, x0 * (1 if a >= 0 else -1), y0 * (1 if b >= 0 else -1)


@dataclass(frozen=True)
class Rectangle:
    x0: int
    x1: int
    y0: int
    y1: int
    direction: str

    def corners(self) -> tuple[Point, ...]:
        return ((self.x0, self.y0), (self.x0, self.y1),
                (self.x1, self.y0), (self.x1, self.y1))

    def vertices(self):
        for y in range(self.y0, self.y1 + 1):
            for x in range(self.x0, self.x1 + 1):
                yield x, y


def necklace(u: Point, s: int) -> tuple[list[Point], list[Rectangle]]:
    """Return n centres plus their translated endpoint, and exactly 5n boxes."""
    if s < 8 or s % 2:
        raise ValueError("s must be even and at least 8")
    S = dot(u, u)
    if S < (64 * s) ** 2:
        raise ValueError("this control uses the theorem regime ell >= 64s")
    n = ceil_scaled_sqrt(S, 4, s)
    centres = [(nearest(i * u[0], n), nearest(i * u[1], n))
               for i in range(n + 1)]
    boxes: list[Rectangle] = []
    for i, (x, y) in enumerate(centres[:-1]):
        boxes.extend((
            Rectangle(x - 2*s, x + 2*s, y + s, y + 2*s, 'H'),
            Rectangle(x - 2*s, x + 2*s, y - 2*s, y - s, 'H'),
            Rectangle(x - 2*s, x - s, y - 2*s, y + 2*s, 'V'),
            Rectangle(x + s, x + 2*s, y - 2*s, y + 2*s, 'V'),
        ))
        X, Y = centres[i + 1]
        boxes.append(Rectangle(min(x, X) - 3*s, max(x, X) + 3*s,
                               max(y, Y) - s//2, min(y, Y) + s//2, 'H'))
    return centres, boxes


def coset_key(x: Point, u: Point, v: Point) -> Point:
    N = det(u, v)
    if N <= 0:
        raise ValueError("basis must be positively oriented")
    return det(x, v) % N, det(u, x) % N


def crossing_path(box: Rectangle, salt: int) -> set[Point]:
    """A deterministic, sometimes backtracking, crossing inside one box."""
    if box.direction == 'H':
        y = box.y0 + (salt % (box.y1 - box.y0 + 1))
        pts = {(x, y) for x in range(box.x0, box.x1 + 1)}
        if box.y0 < y < box.y1:
            # Add a spur; it does not impose any stochastic model.
            xm = (box.x0 + box.x1) // 2
            pts.update((xm, yy) for yy in range(y, min(y + 3, box.y1) + 1))
        return pts
    x = box.x0 + (salt % (box.x1 - box.x0 + 1))
    return {(x, y) for y in range(box.y0, box.y1 + 1)}


def has_crossing(box: Rectangle, occupied, key) -> bool:
    """Free-boundary NN crossing BFS, restricted to this one lifted rectangle."""
    if box.direction == 'H':
        starts = [(box.x0, y) for y in range(box.y0, box.y1 + 1)]
        finish = lambda z: z[0] == box.x1
    else:
        starts = [(x, box.y0) for x in range(box.x0, box.x1 + 1)]
        finish = lambda z: z[1] == box.y1
    seen = {z for z in starts if key(z) in occupied}
    todo = deque(seen)
    while todo:
        z = todo.popleft()
        if finish(z):
            return True
        for e in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            t = add(z, e)
            if (box.x0 <= t[0] <= box.x1 and box.y0 <= t[1] <= box.y1
                    and t not in seen and key(t) in occupied):
                seen.add(t)
                todo.append(t)
    return False


def winding_gcd(occupied: dict[Point, Point], u: Point, v: Point) -> tuple[int, int]:
    """Independent graph-potential traversal. Return gcd of u,v windings."""
    N = det(u, v)
    visited: set[Point] = set()
    all_a = all_b = 0
    for initial in occupied:
        if initial in visited:
            continue
        pot = {initial: (0, 0)}
        visited.add(initial)
        todo = deque([initial])
        while todo:
            k = todo.popleft()
            r = occupied[k]
            for e in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                nb = coset_key(add(r, e), u, v)
                if nb not in occupied:
                    continue
                guess = add(pot[k], e)
                if nb not in pot:
                    pot[nb] = guess
                    visited.add(nb)
                    todo.append(nb)
                else:
                    c = guess[0] - pot[nb][0], guess[1] - pot[nb][1]
                    aN, bN = det(c, v), det(u, c)
                    if aN % N or bN % N:
                        raise AssertionError("non-period cycle gain")
                    all_a = gcd(all_a, abs(aN // N))
                    all_b = gcd(all_b, abs(bN // N))
    return all_a, all_b


def control(u0: Point, v0: Point, s: int = 8, seed: int = 0) -> dict:
    u, v = reduced_basis(u0, v0)
    S, N = dot(u, u), det(u, v)
    centres, boxes = necklace(u, s)
    n = len(centres) - 1
    assert centres[0] == (0, 0) and centres[-1] == u
    assert 4*N*N >= 3*S*S  # h >= sqrt(3)*ell/2
    assert (n*s)**2 >= 16*S and (n-1)**2*s*s < 16*S
    assert n*n*s*s <= 25*S
    for i in range(n):
        assert max(abs(centres[i+1][j]-centres[i][j]) for j in (0, 1)) <= s//2
    for box in boxes:
        width, height = box.x1-box.x0, box.y1-box.y0
        assert width > 0 and height > 0 and width*width + height*height < S
        long, short = (width, height) if box.direction == 'H' else (height, width)
        assert long <= 14*short
        assert all(det(u, z)**2 <= 36*s*s*S for z in box.corners())

    occupied: dict[Point, Point] = {}
    key = lambda x: coset_key(x, u, v)
    for j, box in enumerate(boxes):
        for z in crossing_path(box, seed + 17*j):
            occupied.setdefault(key(z), z)
    assert all(has_crossing(box, occupied, key) for box in boxes)
    ga, gb = winding_gcd(occupied, u, v)
    assert (ga, gb) == (1, 0)

    # Check WHOLE event supports for a few translates, not just occupied paths.
    support = {key(z): z for box in boxes for z in box.vertices()}
    g, a, b = bezout(*u)
    transverse_unit = -b, a
    assert det(u, transverse_unit) == g and N % g == 0
    D = ceil_scaled_sqrt(S, 12*s+2, g)
    gap, bands = D*g, N//(D*g)
    assert gap*gap > 144*s*s*S
    assert gap < (12*s+3) * (isqrt(S)+1)
    assert (30*s*bands)**2*S >= N*N
    shifts = sorted(set([0, 1, max(0, bands-1)]))
    supports = []
    for j in shifts:
        shift = j*D*transverse_unit[0], j*D*transverse_unit[1]
        here = {key(add(z, shift)) for z in support.values()}
        assert len(here) == len(support)
        for previous in supports:
            assert not (here & previous)
        supports.append(here)
    return dict(input_basis=[list(u0), list(v0)], reduced_basis=[list(u), list(v)],
                N=N, squared_shortest_period=S, ambient_gcd=g, block_scale=s,
                blocks=n, crossing_events=len(boxes), support_vertices=len(support),
                occupied_witness_vertices=len(occupied), winding_gcd=[ga, gb],
                packed_bands=bands, checked_translates=shifts,
                conditions_pass=True)


def finite_harris_control() -> dict:
    """Independent tiny overlapping-increasing-event Harris arithmetic."""
    # Three overlapping clauses, not a model for the RSW constant.
    events = [lambda m: bool(m & 3), lambda m: bool(m & 6),
              lambda m: bool(m & 5)]
    rows = []
    for p in (Fraction(1, 3), Fraction(1, 2), Fraction(2, 3)):
        probs = [Fraction(0)]*3
        joint = Fraction(0)
        for m in range(8):
            k = bin(m).count('1')
            w = p**k * (1-p)**(3-k)
            vals = [e(m) for e in events]
            for i, ok in enumerate(vals):
                probs[i] += w*ok
            joint += w*all(vals)
        product = probs[0]*probs[1]*probs[2]
        assert joint >= product
        rows.append(dict(p=str(p), joint=str(joint), product=str(product)))
    return dict(scope="toy Harris inequality, not a percolation/RSW estimate", rows=rows)


CASES = [
    ((512, 0), (37, 900)),
    ((0, 512), (-1200, 51)),
    ((384, 512), (-1536, 1152)),
    ((511, 129), (-387, 1533)),
    ((357, -407), (1221, 1071)),
    ((900, 0), (271, 711)),
    ((1100, 0), (473, 853)),
    ((1024, 0), (517, 515)),
]


class GeometryTests(unittest.TestCase):
    def test_integer_rounding_and_endpoints(self):
        for den in range(1, 16):
            for num in range(-40, 41):
                k = nearest(num, den)
                self.assertLessEqual(abs(k*den-num)*2, den)
                self.assertEqual(nearest(num+7*den, den), k+7)

    def test_lagrange_reduction(self):
        for u, v in CASES:
            a, b = reduced_basis(u, v)
            self.assertEqual(det(a, b), abs(det(u, v)))
            self.assertLessEqual(dot(a, a), dot(b, b))
            self.assertLessEqual(2*abs(dot(a, b)), dot(a, a))
            for i in range(-5, 6):
                for j in range(-5, 6):
                    if i or j:
                        x = i*a[0]+j*b[0], i*a[1]+j*b[1]
                        self.assertGreaterEqual(dot(x, x), dot(a, a))

    def test_seam_and_nonprimitive_period(self):
        r = control((384, 512), (-1536, 1152), seed=3)
        self.assertEqual(r['ambient_gcd'], 128)
        self.assertEqual(r['winding_gcd'], [1, 0])

    def test_oblique_primitive_and_packing(self):
        r = control((511, 129), (-387, 1533), seed=11)
        self.assertEqual(r['ambient_gcd'], 1)
        self.assertGreaterEqual(r['packed_bands'], 3)

    def test_reject_too_small_support_regime(self):
        with self.assertRaises(ValueError):
            necklace((4, 4), 8)

    def test_harris_control(self):
        self.assertEqual(len(finite_harris_control()['rows']), 3)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path)
    ap.add_argument('--tests', action='store_true')
    args = ap.parse_args()
    if args.tests:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(GeometryTests)
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        raise SystemExit(not result.wasSuccessful())
    if not args.output:
        ap.error('--output PATH or --tests is required')
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    rows = [control(u, v, seed=i*19) for i, (u, v) in enumerate(CASES)]
    out = dict(schema='matching-one.oblique-necklace-controls.v1',
               nature='deterministic geometric controls, not production evidence',
               cases=rows, cases_count=len(rows),
               crossing_event_checks=sum(r['crossing_events'] for r in rows),
               harris_control=finite_harris_control(),
               no_probability_or_novelty_estimate=True,
               full_repository_ci_run=False)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(cases=len(rows), crossing_checks=out['crossing_event_checks'],
                         output=str(args.output))))


if __name__ == '__main__':
    main()
