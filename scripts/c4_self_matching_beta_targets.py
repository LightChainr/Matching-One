#!/usr/bin/env python3
"""Freeze exact Beta targets for the C4 self-matching quotient.

This module is intentionally independent of the exhaustive target enumerator.
It obtains the geometry hypothesis from the shortest nontrivial cycle in the
lifted checkerboard graph, then degree-elevates ``2 I_p(s,s)-1`` to the
unnormalized N-site Bernstein basis used by the exact enumeration.

The second frozen target is the majority law on the N/2 antipodal orbits:
``2 I_p((N+2)/4,(N+2)/4)-1`` for N=2 mod 4.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from c4_self_matching_exact import c4_self_matching_torus


@dataclass(frozen=True)
class CycleCertificate:
    support: int
    start_vertex: int
    vertices: tuple[int, ...]
    lifted_displacements: tuple[tuple[int, int], ...]
    winding: tuple[int, int]


def shortest_nontrivial_cycle(a: int, b: int) -> CycleCertificate:
    """Return a geometry-only shortest nonzero-homology cycle certificate."""

    geometry = c4_self_matching_torus(a, b)
    adjacency: list[list[tuple[int, int, int]]] = [
        [] for _ in range(geometry.n)
    ]
    for edge in geometry.primal_edges:
        adjacency[edge.i].append((edge.j, edge.dx, edge.dy))
        adjacency[edge.j].append((edge.i, -edge.dx, -edge.dy))

    best: CycleCertificate | None = None
    for start in range(geometry.n):
        initial = (start, 0, 0)
        queue = deque([initial])
        parent: dict[tuple[int, int, int], tuple[int, int, int] | None] = {
            initial: None
        }
        depth = {initial: 0}
        while queue:
            state = queue.popleft()
            vertex, x, y = state
            next_depth = depth[state] + 1
            if best is not None and next_depth > best.support:
                continue
            for target, dx, dy in adjacency[vertex]:
                candidate = (target, x + dx, y + dy)
                if candidate in parent:
                    continue
                parent[candidate] = state
                depth[candidate] = next_depth
                if target == start:
                    winding = geometry.periods.winding((x + dx, y + dy))
                    if winding != (0, 0):
                        path = [candidate]
                        cursor = candidate
                        while parent[cursor] is not None:
                            cursor = parent[cursor]  # type: ignore[assignment]
                            path.append(cursor)
                        path.reverse()
                        vertices = tuple(item[0] for item in path)
                        # A shortest nontrivial cycle must not reuse a quotient
                        # vertex before its closing return; certify that the
                        # edge length is also its occupied support size.
                        if len(set(vertices[:-1])) != next_depth:
                            raise AssertionError(
                                "shortest cycle repeats a quotient vertex"
                            )
                        certificate = CycleCertificate(
                            support=next_depth,
                            start_vertex=start,
                            vertices=vertices,
                            lifted_displacements=tuple(
                                (item[1], item[2]) for item in path
                            ),
                            winding=winding,
                        )
                        if best is None or certificate.support < best.support:
                            best = certificate
                        queue.clear()
                        break
                queue.append(candidate)
    if best is None:
        raise AssertionError("no nontrivial cycle found")
    return best


def beta_matching_bernstein_counts(n: int, parameter: int) -> list[int]:
    """Return unnormalized degree-N Bernstein counts for 2 I_p(r,r)-1."""

    from math import comb

    m = 2 * parameter - 1
    if m > n:
        raise ValueError("Beta polynomial degree exceeds target degree")
    base = [
        (-1 if k < parameter else 1) * comb(m, k)
        for k in range(m + 1)
    ]
    elevated = [0] * (n + 1)
    for k in range(n + 1):
        elevated[k] = sum(
            base[j] * comb(n - m, k - j)
            for j in range(max(0, k - (n - m)), min(m, k) + 1)
        )
    return elevated


def vector_sha256(values: Sequence[int]) -> str:
    canonical = json.dumps(list(values), separators=(",", ":"))
    return hashlib.sha256(canonical.encode("ascii")).hexdigest()


def freeze_payload(a: int = 5, b: int = 1) -> dict[str, object]:
    geometry = c4_self_matching_torus(a, b)
    cycle = shortest_nontrivial_cycle(a, b)
    if geometry.n % 4 != 2:
        raise ValueError("antipodal majority rule requires N=2 mod 4")
    majority_parameter = (geometry.n + 2) // 4

    hypotheses = []
    for name, parameter, provenance in (
        (
            "geometry_shortest_support",
            cycle.support,
            "parameter fixed from the geometry-only nontrivial-cycle support",
        ),
        (
            "antipodal_orbit_majority",
            majority_parameter,
            "parameter fixed as the strict majority of N/2 two-site orbits",
        ),
    ):
        vector = beta_matching_bernstein_counts(geometry.n, parameter)
        hypotheses.append(
            {
                "name": name,
                "beta_parameters": [parameter, parameter],
                "formula": f"M(p)=2*I_p({parameter},{parameter})-1",
                "bernstein_integer_coefficients": vector,
                "bernstein_vector_sha256": vector_sha256(vector),
                "provenance": provenance,
            }
        )

    return {
        "schema": "matching-one/c4-self-matching-beta-targets/v1",
        "status": "FROZEN_BEFORE_N26_ENUMERATION",
        "geometry": {
            "a": a,
            "b": b,
            "N": geometry.n,
            "period_matrix": geometry.periods.matrix,
            "wrapping_channel": "either",
        },
        "geometry_shortest_cycle_certificate": {
            "support": cycle.support,
            "start_vertex": cycle.start_vertex,
            "vertices_including_closing_return": cycle.vertices,
            "lifted_displacements": cycle.lifted_displacements,
            "winding": cycle.winding,
        },
        "antipodal_orbits": geometry.n // 2,
        "scoring_order": [
            "geometry_shortest_support",
            "antipodal_orbit_majority",
        ],
        "hypotheses": hypotheses,
        "post_failure_rule": "STOP_WITHOUT_GENERALIZED_BETA_FIT",
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a", type=int, default=5)
    parser.add_argument("--b", type=int, default=1)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args(argv)
    rendered = json.dumps(
        freeze_payload(args.a, args.b), indent=2, sort_keys=True
    ) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
