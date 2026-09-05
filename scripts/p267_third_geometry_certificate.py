#!/usr/bin/env python3
"""Finite N50/N100/N125 geometry certificate; no stochastic calculation.

The fixed map is O=(1/5)[[4,-3],[3,4]].  Every Smith-(5,N/5)
target has P'=5B and det(B)=N/25, so bounded HNF enumeration is exhaustive.
All geometry arithmetic is rational.  E4 values are explicitly numerical
conditional-model coordinates, not measured responses or field identities.
"""
from __future__ import annotations

import argparse
import cmath
from fractions import Fraction as F
import json
from math import floor, gcd
from pathlib import Path

from integer_period_torus import IntegerPeriods


O = ((F(4, 5), F(-3, 5)), (F(3, 5), F(4, 5)))
M = ((4, 3), (-3, 4))  # 5 O^{-1}
SOURCE_COMMIT = "ce01e4d10abb03abf1b278192510937ee96db29d"
CELLS = (
    (50, "i", ((7, -1), (1, 7)), (F(0), F(1))),
    (50, "2i", ((3, -8), (4, 6)), (F(0), F(2))),
    (75, "3i", ((4, 9), (-3, 12)), (F(0), F(3))),
    (75, "1/2+3i/2", ((7, 2), (1, 11)), (F(1, 2), F(3, 2))),
    (100, "2i", ((7, -2), (1, 14)), (F(0), F(2))),
    (100, "4i", ((4, 12), (-3, 16)), (F(0), F(4))),
    (100, "1/2+i", ((8, 10), (-6, 5)), (F(1, 2), F(1))),
    (125, "i", ((11, 2), (-2, 11)), (F(0), F(1))),
    (125, "5i", ((4, 15), (-3, 20)), (F(0), F(5))),
    (125, "1/2+5i/2", ((7, 1), (1, 18)), (F(1, 2), F(5, 2))),
)


def multiply(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def determinant(p):
    return p[0][0] * p[1][1] - p[0][1] * p[1][0]


def smith(p):
    first = gcd(*(int(x) for row in p for x in row))
    return first, determinant(p) // first


def gram(p):
    return tuple(tuple(sum(p[k][i] * p[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def tau(p):
    g = gram(p)
    return F(g[0][1], g[0][0]), F(determinant(p), g[0][0])


def reduced_shape(z):
    """PSL2Z fundamental domain with complex conjugation/reflection identified."""
    x, y = z
    for _ in range(100):
        x -= floor(x + F(1, 2))
        radius = x * x + y * y
        if radius < 1 or (radius == 1 and x < 0):
            x, y = -x / radius, y / radius
        else:
            return abs(x), y
    raise ArithmeticError("nonterminating reduction")


def chi4(p):
    x, y = p[0][0], p[1][0]
    return F(x**4 - 6*x*x*y*y + y**4, (x*x+y*y)**2)


def complete_family(n):
    """Enumerate target sublattices, not arbitrary finite boxes of matrices."""
    if n % 25:
        return []
    m, rows = n // 25, []
    for a in range(1, m + 1):
        if m % a:
            continue
        d = m // a
        for b in range(a):
            target_basis = ((a, b), (0, d))
            source = multiply(M, target_basis)
            if gcd(a, b, d) != 1 or smith(source) != (1, n):
                continue
            rows.append({"target_divided_by_5_hnf": target_basis,
                         "source_periods": source,
                         "shape": reduced_shape(tau(source))})
    return rows


def sigma3(n):
    return sum(d**3 for d in range(1, n+1) if n % d == 0)


def e4(z):
    # At all selected cells Im(tau)>=1; omitted-series error is <1e-101.
    return 1 + 240 * sum(sigma3(n) * cmath.exp(2j*cmath.pi*n*z)
                        for n in range(1, 41))


def shape_coordinate(z):
    """Safe only for these reflection-boundary cells, where E4 is real."""
    x, y = z
    if x not in (F(0), F(1, 2)):
        raise ValueError("general shape needs complex spin-four projection")
    return float(y)**2 * e4(complex(float(x), float(y))).real / e4(1j).real


def serialize(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): serialize(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(v) for v in value]
    return value


def certificate():
    families = {}
    for n in (25, 50, 75, 100, 125):
        rows = complete_family(n)
        shapes = sorted({row["shape"] for row in rows})
        families[n] = {"source_smith": [1, n], "target_smith": [5, n//5],
                       "admissible_hnf_count": len(rows),
                       "distinct_shape_count": len(shapes), "shapes": shapes,
                       "complete_hnf_rows": rows}
    assert [families[n]["distinct_shape_count"] for n in families] == [1, 2, 2, 3, 3]
    cells = []
    for n, name, p, expected_tau in CELLS:
        rotated = multiply(O, p)
        assert all(x.denominator == 1 for row in rotated for x in row)
        q = tuple(tuple(int(x) for x in row) for row in rotated)
        assert determinant(p) == determinant(q) == n
        assert smith(p) == (1, n) and smith(q) == (5, n//5)
        assert gram(p) == gram(q)
        assert tau(p) == tau(q) == expected_tau
        assert reduced_shape(expected_tau) in families[n]["shapes"]
        delta = chi4(p) - chi4(q)
        assert delta != 0
        for matrix in (p, q):
            periods = IntegerPeriods(matrix)
            assert periods.order == n
            # Immediate NN and matching neighbours are distinct, with no loops.
            shifts = ((1, 0), (-1, 0), (0, 1), (0, -1),
                      (1, 1), (1, -1), (-1, 1), (-1, -1))
            assert len({periods.quotient_key(v) for v in ((0, 0), *shifts)}) == 9
        g = shape_coordinate(expected_tau)
        t = (g - 1) / float(F(7, 4))
        cells.append({"N": n, "tau": name, "tau_exact": expected_tau,
                      "first_period_matrix": p, "second_period_matrix": q,
                      "first_smith": smith(p), "second_smith": smith(q),
                      "gram": gram(p), "shortest_period_squared": gram(p)[0][0],
                      "first_chi4": chi4(p), "second_chi4": chi4(q),
                      "delta_chi4": delta, "spin4_projector_degenerate": False,
                      "immediate_graph_alias": False,
                      "conditional_E4_shape_g": g,
                      "affine_E4_old_i_weight": 1-t,
                      "affine_E4_old_2i_weight": t,
                      "runner_matrix_arguments": ["--first-matrix", *sum(p, ()),
                                                  "--second-matrix", *sum(q, ())]})
    g = lambda x, y: shape_coordinate((F(x), F(y)))
    return serialize({
        "schema": "matching-one/p267-third-geometry-feasibility/v1",
        "status": "exact finite-geometry certificate plus conditional model coordinates; no MC",
        "factorial_source_commit": SOURCE_COMMIT,
        "map": O,
        "proof": "P'=5B, det(B)=N/25; primitive B plus gcd((5O^-1)B)=1 are necessary and sufficient. Column HNF enumerates all target sublattices.",
        "N50_third_shape_impossible_even_if_map_changes": True,
        "minimum_N_with_three_shapes_in_fixed_map_smith_family": 100,
        "families": families, "cells": cells,
        "conditional_covector": {
            "basis_invariant_complex": "Phi4(P)=Im(tau)^2 exp(-4i arg(omega1)) E4(tau)",
            "normalized_difference": "Re(Phi4(P)-Phi4(OP))/(chi4(P)-chi4(OP))/E4(i)",
            "boundary_simplification": "g(tau)=Im(tau)^2 E4(tau)/E4(i) when Re(tau)=0 or 1/2",
            "third_cell_predictor": "v(N,tau)=s_N*((1-t)*v(50,i)+t*v(50,2i)), t=(g(tau)-1)/(7/4)",
            "scale_boundary": "s_N is an explicit additional hypothesis or requires a new-area bridge cell; no scalar or exponent is inferred here",
        },
        "exact_conditional_E4_relations": [
            {"identity": "g(2i)=11/4", "rational_rhs": F(11, 4),
             "display_residual": g(0, 2)-float(F(11, 4))},
            {"identity": "g(4i)+g(1/2+i)=91/8", "rational_rhs": F(91, 8),
             "display_residual": g(0, 4)+g(F(1, 2), 1)-float(F(91, 8))},
            {"identity": "g(5i)-g(1/2+5i/2)=322/25", "rational_rhs": F(322, 25),
             "display_residual": g(0, 5)-g(F(1, 2), F(5, 2))-float(F(322, 25))},
        ],
        "scientific_boundary": "Feasibility and normalization are exact. E4 predictions are conditional, not established amplitudes for A/E/C/W, and do not identify a continuum field. Finite widths differ. No samples or additional model fitting.",
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    result = json.dumps(certificate(), indent=2) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(result)
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
