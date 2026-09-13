#!/usr/bin/env python3
"""Verify two-cycle blocker certificates; no discovery search or Monte Carlo."""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

SCHEMA = "matching-one/p429-dual-cycle-blocker/v1"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "results/p429-dual-cycle-blocker/certificate.json"
NN = ((1, 0), (-1, 0), (0, 1), (0, -1))
MATCHING = NN + ((1, 1), (1, -1), (-1, 1), (-1, -1))

class Geometry:
    def __init__(self, matrix: Sequence[Sequence[int]]) -> None:
        if len(matrix) != 2 or any(len(row) != 2 for row in matrix):
            raise ValueError("a 2-by-2 HNF matrix is required")
        a, b, c, d = (matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1])
        if any(type(x) is not int for x in (a, b, c, d)):
            raise ValueError("period entries must be integers")
        if a <= 0 or d <= 0 or c != 0 or not 0 <= b < a:
            raise ValueError("positive column-HNF periods required")
        self.a, self.b, self.d, self.n = a, b, d, a*d
        self.full = (1 << self.n) - 1
        self.adjacency = [[(self.target(v, step), step) for step in MATCHING]
                          for v in range(self.n)]
        self.edges = {
            False: [(v, self.target(v, step), step) for v in range(self.n)
                    for step in ((1, 0), (0, 1))],
            True: [(v, self.target(v, step), step) for v in range(self.n)
                   for step in ((1, 0), (0, 1), (1, 1), (-1, 1))],
        }

    def target(self, v: int, step: Tuple[int, int]) -> int:
        x, y = v % self.a + step[0], v // self.a + step[1]
        shift, ry = divmod(y, self.d)
        return (x - self.b*shift) % self.a + self.a*ry

    def winding(self, displacement: Tuple[int, int]) -> Tuple[int, int]:
        x, y = displacement
        if y % self.d:
            raise ValueError("unclosed vertical displacement")
        wy = y // self.d
        if (x-self.b*wy) % self.a:
            raise ValueError("unclosed horizontal displacement")
        return ((x-self.b*wy)//self.a, wy)

def ambient_rank(mask: int, geometry: Geometry, matching: bool) -> int:
    """Independent potential-union-find; discovery used lifted BFS instead."""
    if type(mask) is not int or mask < 0 or mask > geometry.full:
        raise ValueError("mask is outside the quotient")
    parent = list(range(geometry.n))
    size = [1]*geometry.n
    px = [0]*geometry.n
    py = [0]*geometry.n

    def find(v: int) -> Tuple[int, int, int]:
        if parent[v] != v:
            old = parent[v]
            root, x, y = find(old)
            px[v] += x
            py[v] += y
            parent[v] = root
        return parent[v], px[v], py[v]

    first = None
    for u, v, (dx, dy) in geometry.edges[matching]:
        if not (mask >> u & 1 and mask >> v & 1):
            continue
        ru, ux, uy = find(u)
        rv, vx, vy = find(v)
        ex, ey = ux+dx-vx, uy+dy-vy
        if ru != rv:
            if size[ru] < size[rv]:
                parent[ru] = rv
                px[ru], py[ru] = -ex, -ey
                size[rv] += size[ru]
            else:
                parent[rv] = ru
                px[rv], py[rv] = ex, ey
                size[ru] += size[rv]
        elif ex or ey:
            w = geometry.winding((ex, ey))
            if first is None:
                first = w
            elif first[0]*w[1] != first[1]*w[0]:
                return 2
    return int(first is not None)

def check_cycle(cycle: Dict[str, Any], white: int, geometry: Geometry) -> set:
    vertices = cycle["vertices"]
    if len(vertices) < 3 or len(set(vertices)) != len(vertices):
        raise ValueError("cycle must be simple with at least three vertices")
    if any(type(v) is not int or not 0 <= v < geometry.n or not (white >> v & 1)
           for v in vertices):
        raise ValueError("cycle contains a nonwhite or invalid vertex")
    dx = dy = 0
    for u, v in zip(vertices, vertices[1:]+vertices[:1]):
        steps = [step for target, step in geometry.adjacency[u] if target == v]
        if len(steps) != 1:
            raise ValueError("cycle edge is absent or has ambiguous lift")
        dx += steps[0][0]
        dy += steps[0][1]
    winding = geometry.winding((dx, dy))
    if winding == (0, 0) or list(winding) != cycle["winding"]:
        raise ValueError("cycle does not have the declared nonzero winding")
    return set(vertices)

def verify_case(case: Dict[str, Any], full_pairs: bool = False) -> Dict[str, Any]:
    geometry = Geometry(case["periods"])
    black = int(case["occupied_mask_hex"], 16)
    if black < 0 or black > geometry.full:
        raise ValueError("invalid occupied mask")
    white = geometry.full ^ black
    if bin(black).count("1") != case["occupied_count"]:
        raise ValueError("occupied count mismatch")
    if (ambient_rank(black, geometry, False), ambient_rank(white, geometry, True)) != (1, 1):
        raise ValueError("both complementary carriers must have ambient rank one")
    vacancies = [v for v in range(geometry.n) if white >> v & 1]
    singles = set()
    for v in vacancies:
        rb = ambient_rank(black | (1 << v), geometry, False)
        rw = ambient_rank(white ^ (1 << v), geometry, True)
        if rb+rw != 2:
            raise ValueError("single-insertion duality failed")
        if rb == 2:
            singles.add(v)
    if sorted(singles) != case["singleton_triggers"]:
        raise ValueError("singleton trigger set mismatch")
    if len(case["cycles"]) != 2:
        raise ValueError("exactly two cycle witnesses required")
    c0, c1 = [check_cycle(c, white, geometry) for c in case["cycles"]]
    if not (c0 & c1) <= singles:
        raise ValueError("cycles overlap on an individually safe vertex")
    answer = {"id": case["id"], "cycle_lengths": [len(c0), len(c1)],
              "singleton_count": len(singles), "two_cycle_certificate": True}
    if not full_pairs:
        return answer
    safe = [v for v in vacancies if v not in singles]
    degrees = Counter()
    count = checks = 0
    for i, u in enumerate(safe):
        for v in safe[i+1:]:
            pair = (1 << u) | (1 << v)
            rb = ambient_rank(black | pair, geometry, False)
            rw = ambient_rank(white ^ pair, geometry, True)
            if rb+rw != 2:
                raise ValueError("pair duality failed")
            checks += 1
            if rb == 2:
                if not ((u in c0 and v in c1) or (u in c1 and v in c0)):
                    raise ValueError("trigger edge does not cross cycle sides")
                count += 1
                degrees[u] += 1
                degrees[v] += 1
    measured = {"trigger_edges": count,
                "two_stars": sum(d*(d-1)//2 for d in degrees.values()),
                "side_sizes": [sum(v in degrees for v in c0-singles),
                               sum(v in degrees for v in c1-singles)],
                "isolated_safe_vertices": sum(v not in degrees for v in safe)}
    if measured != case["expected_pair_graph"]:
        raise ValueError("recomputed trigger graph summary differs from certificate")
    answer.update(measured)
    answer["all_safe_pair_checks"] = checks
    return answer

def verify(payload: Dict[str, Any], full_pairs: bool = False) -> List[Dict[str, Any]]:
    if payload.get("schema") != SCHEMA or not payload.get("cases"):
        raise ValueError("unknown or empty certificate")
    ids = [case["id"] for case in payload["cases"]]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate case id")
    return [verify_case(case, full_pairs) for case in payload["cases"]]

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path, nargs="?", default=DEFAULT)
    parser.add_argument("--full-pairs", action="store_true")
    args = parser.parse_args()
    print(json.dumps(verify(json.loads(args.certificate.read_text()), args.full_pairs), indent=2))

if __name__ == "__main__":
    main()
