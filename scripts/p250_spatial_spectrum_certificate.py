#!/usr/bin/env python3
"""Exact N505 witnesses for the finite-group spatial spectrum of P250.

No Monte Carlo, repository imports, matrix eigensolver, or floating arithmetic.
The accompanying note supplies the Fourier/cyclotomic proof. This executable
constructs its finite topology witnesses, not statistical spectral weights.
"""
from __future__ import annotations

import argparse
from collections import deque
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

Point = Tuple[int, int]
NN = ((1, 0), (-1, 0), (0, 1), (0, -1))
MATCHING = NN + ((1, 1), (1, -1), (-1, 1), (-1, -1))
SCHEMA = "matching-one/p250-spatial-spectrum-certificate/v1"
SOURCE_COMMIT = "33c557b9aebed1bc9c07019b9cd5cee6c04be947"
HANDS = (("plus", 19, 12), ("minus", 21, -8))


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def digest(value: object) -> str:
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def parent_label(x: int, y: int) -> int:
    """The exact isomorphism Z^2/(10+i) -> Z/101."""
    return (x - 10 * y) % 101


class GaussianQuotient:
    """Adjugate-pair quotient keys; lifted BFS supplies cycle displacements."""

    def __init__(self, a: int, b: int) -> None:
        if not isinstance(a, int) or not isinstance(b, int) or a * a + b * b == 0:
            raise ValueError("a nonzero integer Gaussian period is required")
        self.a, self.b = a, b
        self.n = a * a + b * b
        self.representatives: Dict[Point, Point] = {(0, 0): (0, 0)}
        queue = deque([(0, 0)])
        while queue:
            vertex = queue.popleft()
            x, y = self.representatives[vertex]
            for dx, dy in NN:
                point = (x + dx, y + dy)
                target = self.key(*point)
                if target not in self.representatives:
                    self.representatives[target] = point
                    queue.append(target)
        if len(self.representatives) != self.n:
            raise AssertionError("adjugate quotient order mismatch")

    def key(self, x: int, y: int) -> Point:
        return ((self.a * x + self.b * y) % self.n,
                (-self.b * x + self.a * y) % self.n)

    def component_ranks(self, enabled: Iterable[Point], steps: Sequence[Point]):
        enabled = set(enabled)
        if not enabled <= self.representatives.keys():
            raise ValueError("enabled set contains a nonvertex")
        seen = set()
        ranks: Dict[Point, int] = {}
        components = []
        for start in self.representatives:
            if start not in enabled or start in seen:
                continue
            lifts = {start: self.representatives[start]}
            queue = deque([start])
            cycles = []
            seen.add(start)
            while queue:
                vertex = queue.popleft()
                x, y = lifts[vertex]
                for dx, dy in steps:
                    point = (x + dx, y + dy)
                    target = self.key(*point)
                    if target not in enabled:
                        continue
                    if target not in lifts:
                        lifts[target] = point
                        seen.add(target)
                        queue.append(target)
                    else:
                        old = lifts[target]
                        cycle = (point[0] - old[0], point[1] - old[1])
                        if cycle != (0, 0):
                            if self.key(*cycle) != (0, 0):
                                raise AssertionError("cycle displacement is not a period")
                            cycles.append(cycle)
            rank = 0 if not cycles else 1
            if cycles and any(cycles[0][0] * v[1] != cycles[0][1] * v[0] for v in cycles):
                rank = 2
            ranks.update((vertex, rank) for vertex in lifts)
            components.append({"size": len(lifts), "ambient_rank": rank})
        return sorted(components, key=lambda row: (row["size"], row["ambient_rank"])), ranks


def staircase(a: int, b: int) -> List[Point]:
    if a <= 0 or b == 0:
        raise ValueError("this witness requires a>0 and b!=0")
    x, y = 0, 0
    points = []
    for dx, dy in [(1, 0)] * a + [(0, 1 if b > 0 else -1)] * abs(b):
        points.append((x, y))
        x, y = x + dx, y + dy
    if (x, y) != (a, b):
        raise AssertionError("staircase endpoint mismatch")
    return points


def charge_coefficients(values: Sequence[int], charge: int) -> List[int]:
    """Coefficients of 5*F_r in basis 1,zeta5,zeta5^2,zeta5^3."""
    if len(values) != 5 or charge not in (1, 2, 3, 4):
        raise ValueError("five real integer fiber values and nonzero charge required")
    if any(type(value) is not int for value in values):
        raise ValueError("fiber values must be integers")
    coefficients = [0] * 5
    for fiber, value in enumerate(values):
        coefficients[(-charge * fiber) % 5] += value
    return [coefficients[i] - coefficients[4] for i in range(4)]


def phase_multiply(coefficients: Sequence[int], power: int) -> List[int]:
    """Exact fifth-root gauge action on a reduced coefficient vector."""
    if len(coefficients) != 4:
        raise ValueError("four cyclotomic coefficients required")
    expanded = [0] * 5
    for i, value in enumerate(coefficients):
        expanded[(i + power) % 5] += value
    return [expanded[i] - expanded[4] for i in range(4)]


def witness(hand: str, a: int, b: int) -> dict:
    geometry = GaussianQuotient(a, b)
    path = staircase(a, b)
    black = {geometry.key(*point) for point in path}
    white = set(geometry.representatives) - black
    bc, br = geometry.component_ranks(black, NN)
    wc, wr = geometry.component_ranks(white, MATCHING)
    if len(black) != a + abs(b) or geometry.key(a, b) != (0, 0):
        raise AssertionError("staircase is not a simple closed witness")
    scalar = {
        vertex: (1 if br[vertex] == 1 else 0) if vertex in black
        else (-1 if wr[vertex] == 1 else 0)
        for vertex in geometry.representatives
    }
    labels = [geometry.key(j + 10 * fiber, fiber)
              for j in range(101) for fiber in range(5)]
    if len(set(labels)) != 505 or len(scalar) != 505:
        raise AssertionError("parent/fiber section is not a bijection")
    charges = {}
    for charge in (1, 2, 3, 4):
        coefficients = [
            charge_coefficients([scalar[geometry.key(j + 10 * f, f)] for f in range(5)], charge)
            for j in range(101)
        ]
        zeros = [j for j, row in enumerate(coefficients) if not any(row)]
        nonzeros = [j for j, row in enumerate(coefficients) if any(row)]
        if not zeros or not nonzeros:
            raise AssertionError("a gauge-robust nonconstant field was not found")
        gauge_support = all(bool(any(phase_multiply(row, t))) == bool(any(row))
                            for row in coefficients for t in range(5))
        if not gauge_support:
            raise AssertionError("a fifth-root gauge changed zero support")
        charges[str(charge)] = {
            "zero_parents": len(zeros), "nonzero_parents": len(nonzeros),
            "zero_parent": zeros[0], "nonzero_parent": nonzeros[0],
            "zero_coefficients": coefficients[zeros[0]],
            "nonzero_coefficients": coefficients[nonzeros[0]],
            "all_coefficients_sha256": digest(coefficients),
            "all_fifth_root_gauges_preserve_zero_support": gauge_support,
        }
    return {
        "hand": hand, "gaussian_period": [a, b], "order": geometry.n,
        "period_matrix_columns": [[a, b], [-b, a]],
        "occupied_lifts": [list(point) for point in path],
        "occupied_count": len(black), "vacant_count": len(white),
        "black_NN_components": bc, "white_matching_components": wc,
        "parent_fiber_label_count": len(set(labels)),
        "charged_field_basis": "5*F_r=sum_{j=0}^3 coefficients[j]*zeta5^j",
        "charges": charges,
    }


def alias_audit(radius: int) -> dict:
    if type(radius) is not int or radius < 0:
        raise ValueError("radius must be a nonnegative integer")
    classes: Dict[int, List[List[int]]] = {}
    for a in range(-radius, radius + 1):
        for b in range(-radius, radius + 1):
            if abs(a) + abs(b) <= radius:
                classes.setdefault(parent_label(a, b), []).append([a, b])
    return {"radius": radius, "displacement_labels": sum(map(len, classes.values())),
            "distinct_parent_vertices": len(classes),
            "repeated_classes": [{"parent_label": label, "displacements": rows}
                                 for label, rows in sorted(classes.items()) if len(rows) > 1]}


def totient(n: int) -> int:
    if n <= 0:
        raise ValueError("positive integer required")
    return sum(gcd(k, n) == 1 for k in range(1, n + 1))


def matrix_product(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), Fraction(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def matrix_difference(a, b):
    return [[x - y for x, y in zip(arow, brow)] for arow, brow in zip(a, b)]


def projection_control() -> dict:
    """Commuting involutions can have noncommuting orthogonal compressions."""
    identity = [[Fraction(int(i == j)) for j in range(3)] for i in range(3)]
    q = [[Fraction(1, 3)] * 3 for _ in range(3)]
    p = matrix_difference(identity, q)
    u = [[Fraction((1, -1, 1)[i] if i == j else 0) for j in range(3)] for i in range(3)]
    v = [[Fraction((1, 1, -1)[i] if i == j else 0) for j in range(3)] for i in range(3)]
    mm = matrix_product
    pu, pv = mm(mm(p, u), p), mm(mm(p, v), p)
    commutator = matrix_difference(mm(pu, pv), mm(pv, pu))
    leakage = matrix_difference(mm(mm(mm(mm(p, v), q), u), p),
                                mm(mm(mm(mm(p, u), q), v), p))
    checks = {"microscopic_commutation": mm(u, v) == mm(v, u),
              "orthogonal_projector": mm(p, p) == p and list(map(list, zip(*p))) == p,
              "compression_commutator_nonzero": any(any(row) for row in commutator),
              "leakage_identity": commutator == leakage}
    if not all(checks.values()):
        raise AssertionError("projection leakage control failed")
    return {"checks": checks, "commutator": [[str(x) for x in row] for row in commutator],
            "boundary": "Algebraic control only; not a fitted P250 ordered response."}


def build_certificate() -> dict:
    return {
        "schema": SCHEMA, "issue": 406, "data_class": "exact finite witness; no Monte Carlo",
        "source_semantics_commit": SOURCE_COMMIT,
        "parent": {"gaussian_period": [10, 1], "order": 101,
                   "label": "(x-10*y) mod 101", "fiber_step": [10, 1]},
        "cyclotomic_degree": {"phi_5": totient(5), "phi_101": totient(101),
                              "phi_505": totient(505), "extension_degree": 100},
        "witnesses": [witness(*hand) for hand in HANDS],
        "endpoint_aliases": [alias_audit(radius) for radius in (4, 5, 6)],
        "projection_leakage_control": projection_control(),
        "conclusion": {"positive_nonzero_spatial_frequencies_per_hand_charge": 100,
                       "complete_spatial_rank_lower_bound": 100,
                       "ensemble": "0<p<1 independent Bernoulli and uniform parent-anchor target",
                       "boundary": "Not a statistical effective-rank or physical-field count; no production spectrum estimated."},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(build_certificate(), sort_keys=True, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
