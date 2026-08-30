#!/usr/bin/env python3
"""Exact C4 character obstruction for a cyclic length-five Laurent quotient.

Conditional model theorem, not a statement that P250 has such a quotient.
Use R U R^-1=V, R V R^-1=U^-1, R b=b. A reduced quotient with five distinct
joint spectral points must have C4 character (5,1,1,1). The explicit
nonreduced control below instead has character (5,-1,1,-1).
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from typing import Sequence

Matrix = tuple[tuple[F, ...], ...]


def mat(rows: Sequence[Sequence[object]]) -> Matrix:
    if len(rows) != 5 or any(len(row) != 5 for row in rows):
        raise ValueError("expected a five by five exact matrix")
    if any(isinstance(x, (float, bool)) for row in rows for x in row):
        raise TypeError("use rational strings or integers")
    return tuple(tuple(F(x) for x in row) for row in rows)


def eye() -> Matrix:
    return tuple(tuple(F(i == j) for j in range(5)) for i in range(5))


def mul(a: Matrix, b: Matrix) -> Matrix:
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(5))
                       for j in range(5)) for i in range(5))


def add(*terms: tuple[F, Matrix]) -> Matrix:
    return tuple(tuple(sum(scale * value[i][j] for scale, value in terms)
                       for j in range(5)) for i in range(5))


def power(a: Matrix, exponent: int) -> Matrix:
    result = eye()
    for _ in range(exponent):
        result = mul(result, a)
    return result


def trace(a: Matrix) -> F:
    return sum(a[i][i] for i in range(5))


def apply(a: Matrix, vector: Sequence[F]) -> tuple[F, ...]:
    return tuple(sum(a[i][j] * vector[j] for j in range(5)) for i in range(5))


def rank_columns(columns: Sequence[Sequence[F]]) -> int:
    rows = [list(row) for row in zip(*columns)]
    rank = 0
    for col in range(len(columns)):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        divisor = rows[rank][col]
        rows[rank] = [x / divisor for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank:
                factor = rows[i][col]
                rows[i] = [x - factor * y for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def serial(a: Matrix) -> list[list[str]]:
    return [[str(x) for x in row] for row in a]


def build_certificate() -> dict:
    # Algebra: C[x,y] / ((x,y)^3, x^2+y^2), basis 1,x,y,x^2,xy.
    x = [[0] * 5 for _ in range(5)]
    y = [[0] * 5 for _ in range(5)]
    x[1][0], x[3][1], x[4][2] = 1, 1, 1
    y[2][0], y[4][1], y[3][2] = 1, 1, -1
    x, y = mat(x), mat(y)
    r = mat(((1, 0, 0, 0, 0), (0, 0, -1, 0, 0),
             (0, 1, 0, 0, 0), (0, 0, 0, -1, 0), (0, 0, 0, 0, -1)))
    u = add((F(1), eye()), (F(1), x), (F(1, 2), power(x, 2)))
    v = add((F(1), eye()), (F(1), y), (F(1, 2), power(y, 2)))
    inverse_u = add((F(1), eye()), (F(-1), x), (F(1, 2), power(x, 2)))
    b = (F(1), F(0), F(0), F(0), F(0))
    zero = tuple(tuple(F(0) for _ in range(5)) for _ in range(5))
    checks = {
        "log_generators_commute": mul(x, y) == mul(y, x),
        "quadratic_relation": add((F(1), power(x, 2)), (F(1), power(y, 2))) == zero,
        "all_cubic_monomials_vanish": all(mul(power(x, a), power(y, 3-a)) == zero for a in range(4)),
        "nonzero_nilpotent": x != zero and power(x, 3) == zero,
        "rotation_order_four": power(r, 4) == eye(),
        "unit_source_invariant": apply(r, b) == b,
        "translation_inverse": mul(u, inverse_u) == eye(),
        "translations_commute": mul(u, v) == mul(v, u),
        "rotation_U_to_V": mul(mul(r, u), power(r, 3)) == v,
        "rotation_V_to_U_inverse": mul(mul(r, v), power(r, 3)) == inverse_u,
        "cyclic_Laurent_quotient": rank_columns([b, apply(u, b), apply(v, b),
                                                apply(power(u, 2), b), apply(mul(u, v), b)]) == 5,
    }
    character = tuple(trace(power(r, k)) for k in range(4))
    checks["nonreduced_character"] = character == (5, -1, 1, -1)

    # Independent reduced control: one fixed point and a four-point orbit.
    points = [(F(1), F(1)), (F(2), F(3)), (F(3), F(1, 2)),
              (F(1, 2), F(1, 3)), (F(1, 3), F(2))]
    ordinary_u = tuple(tuple(points[i][0] if i == j else F(0) for j in range(5)) for i in range(5))
    ordinary_v = tuple(tuple(points[i][1] if i == j else F(0) for j in range(5)) for i in range(5))
    rr = [[0] * 5 for _ in range(5)]
    for source, target in enumerate((0, 4, 1, 2, 3)):
        rr[target][source] = 1
    rr = mat(rr)
    ordinary_inverse_u = tuple(tuple(1 / points[i][0] if i == j else F(0) for j in range(5)) for i in range(5))
    ones = (F(1),) * 5
    reduced_character = tuple(trace(power(rr, k)) for k in range(4))
    checks.update({
        "reduced_points_distinct": len(set(points)) == 5,
        "reduced_rotation_order_four": power(rr, 4) == eye(),
        "reduced_unit_invariant": apply(rr, ones) == ones,
        "reduced_rotation_U_to_V": mul(mul(rr, ordinary_u), power(rr, 3)) == ordinary_v,
        "reduced_rotation_V_to_U_inverse": mul(mul(rr, ordinary_v), power(rr, 3)) == ordinary_inverse_u,
        "reduced_cyclic": rank_columns([apply(power(ordinary_u, k), ones) for k in range(5)]) == 5,
        "reduced_character": reduced_character == (5, 1, 1, 1),
    })
    # Fixed spectral points: (+1,+1),(-1,-1). One two-cycle: (+1,-1),(-1,+1).
    # All remaining C4 orbits have length four. A reduced support uses a point once.
    decompositions = [(fixed, pairs, quartets)
                      for fixed in range(3) for pairs in range(2) for quartets in range(2)
                      if fixed + 2 * pairs + 4 * quartets == 5]
    checks["only_reduced_orbit_decomposition"] = decompositions == [(1, 0, 1)]
    if not all(checks.values()):
        raise AssertionError([name for name, passed in checks.items() if not passed])
    return {
        "schema": "matching-one/c4-five-state-obstruction/v1",
        "kind": "exact_conditional_model_control_not_P250_data",
        "basis": ["1", "x", "y", "x^2", "xy"],
        "hypotheses": ["flat cyclic length-five quotient of C[U^+-1,V^+-1]",
                       "commuting invertible U,V", "R^4=I", "R U R^-1=V",
                       "R V R^-1=U^-1", "R b=b with b cyclic"],
        "nonreduced_control": {"log_U": serial(x), "log_V": serial(y),
                               "U": serial(u), "V": serial(v), "R": serial(r),
                               "character_R_powers": [int(value) for value in character],
                               "multiplicities_1_i_minus1_minusi": [1, 1, 2, 1]},
        "reduced_control": {"joint_points": [[str(a), str(b)] for a, b in points],
                            "R": serial(rr),
                            "character_R_powers": [int(value) for value in reduced_character],
                            "multiplicities_1_i_minus1_minusi": [2, 1, 1, 1]},
        "reduced_support_orbit_counts_1_2_4": decompositions,
        "exact_checks": checks,
        "passed_exact_checks": len(checks),
        "boundary": "Not valid without flatness, cyclicity, the declared global Laurent action, and the source character. Matching the allowed character is necessary, not sufficient, for reducedness.",
    }


def verify_certificate(payload: dict) -> dict:
    """Verify the stored matrices without invoking the certificate builder."""
    def require(condition: bool, message: str) -> None:
        if not condition:
            raise ValueError(message)

    require(payload["schema"] == "matching-one/c4-five-state-obstruction/v1", "schema")
    control = payload["nonreduced_control"]
    x, y = mat(control["log_U"]), mat(control["log_V"])
    u, v, r = mat(control["U"]), mat(control["V"]), mat(control["R"])
    zero = add((F(0), eye()))
    b = (F(1), F(0), F(0), F(0), F(0))
    require(mul(x, y) == mul(y, x), "log commutation")
    require(add((F(1), power(x, 2)), (F(1), power(y, 2))) == zero, "quadratic relation")
    require(all(mul(power(x, a), power(y, 3-a)) == zero for a in range(4)), "cubic relations")
    require(power(x, 2) != zero, "nilpotency really has order three")
    for translation, log_translation in ((u, x), (v, y)):
        require(translation == add((F(1), eye()), (F(1), log_translation),
                                   (F(1, 2), power(log_translation, 2))), "exponential")
    inverse_u = add((F(1), eye()), (F(-1), x), (F(1, 2), power(x, 2)))
    require(mul(u, inverse_u) == eye(), "translation inverse")
    require(mul(u, v) == mul(v, u), "translation commutation")
    require(power(r, 4) == eye(), "rotation order")
    require(apply(r, b) == b, "unit source")
    require(mul(mul(r, u), power(r, 3)) == v, "R U R^-1=V")
    require(mul(mul(r, v), power(r, 3)) == inverse_u, "R V R^-1=U^-1")
    require(rank_columns([b, apply(u, b), apply(v, b), apply(power(u, 2), b),
                          apply(mul(u, v), b)]) == 5, "cyclic source")
    require([trace(power(r, k)) for k in range(4)] == [5, -1, 1, -1], "nonreduced traces")
    require(control["character_R_powers"] == [5, -1, 1, -1], "stored nonreduced traces")
    require(control["multiplicities_1_i_minus1_minusi"] == [1, 1, 2, 1], "nonreduced multiplicities")

    ordinary = payload["reduced_control"]
    points = [tuple(F(x) for x in point) for point in ordinary["joint_points"]]
    require(len(points) == 5 and all(len(point) == 2 for point in points), "five spectral points")
    require(len(set(points)) == 5 and all(a and b for a, b in points), "distinct torus points")
    diagonal_u = tuple(tuple(points[i][0] if i == j else F(0) for j in range(5)) for i in range(5))
    diagonal_v = tuple(tuple(points[i][1] if i == j else F(0) for j in range(5)) for i in range(5))
    inverse = tuple(tuple(1 / points[i][0] if i == j else F(0) for j in range(5)) for i in range(5))
    rr = mat(ordinary["R"])
    ones = (F(1),) * 5
    require(power(rr, 4) == eye() and apply(rr, ones) == ones, "reduced rotation and source")
    require(mul(mul(rr, diagonal_u), power(rr, 3)) == diagonal_v, "reduced U rotation")
    require(mul(mul(rr, diagonal_v), power(rr, 3)) == inverse, "reduced V rotation")
    require(rank_columns([apply(power(diagonal_u, k), ones) for k in range(5)]) == 5, "reduced cyclicity")
    require([trace(power(rr, k)) for k in range(4)] == [5, 1, 1, 1], "reduced traces")
    require(ordinary["character_R_powers"] == [5, 1, 1, 1], "stored reduced traces")
    require(ordinary["multiplicities_1_i_minus1_minusi"] == [2, 1, 1, 1], "reduced multiplicities")
    expected = [[fixed, pairs, quartets]
                for fixed in range(3) for pairs in range(2) for quartets in range(2)
                if fixed + 2*pairs + 4*quartets == 5]
    require([list(row) for row in payload["reduced_support_orbit_counts_1_2_4"]] == expected == [[1, 0, 1]],
            "reduced support orbit count")
    require(payload["passed_exact_checks"] == 20 and len(payload["exact_checks"]) == 20
            and all(value is True for value in payload["exact_checks"].values()), "stored checks")
    return {"status": "verified_exactly", "model_controls": 2,
            "matrix_dimension": 5, "floating_point_operations": 0}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", help="path to a stored JSON certificate")
    args = parser.parse_args()
    if args.verify:
        with open(args.verify, encoding="utf-8") as source:
            result = verify_certificate(json.load(source))
    else:
        result = build_certificate()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
