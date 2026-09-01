#!/usr/bin/env python3
"""Exact N=9 single-geometry landing control for Issue #537.

The calculation is deliberately bounded.  It reuses the physical torus and
canonical Kreg pair kernel from the thermal-gate audit, sums one clean landing
state over its C4 orbit, and forms the two-channel (kernel/readout) covariance
derivative matrix for q=r-1 and E=q^2.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
THERMAL_GATE = (
    ROOT / "experiments" / "p337-thermal-gate-audit-20260901" / "thermal_gate.py"
)


def _load_thermal_gate():
    spec = importlib.util.spec_from_file_location("p337_thermal_gate", THERMAL_GATE)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load the pinned thermal-gate implementation")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


gate = _load_thermal_gate()


# Polynomials are coefficient tuples in ascending degree.  This tiny exact
# implementation avoids making SymPy a repository dependency.
def trim(a):
    a = tuple(F(x) for x in a)
    while len(a) > 1 and a[-1] == 0:
        a = a[:-1]
    return a


def add(a, b):
    n = max(len(a), len(b))
    return trim(tuple((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)))


def neg(a):
    return tuple(-x for x in a)


def sub(a, b):
    return add(a, neg(b))


def scale(a, c):
    return trim(tuple(F(c) * x for x in a))


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def power(a, n):
    out = (F(1),)
    for _ in range(n):
        out = mul(out, a)
    return out


def evaluate(a, x):
    out = F(0)
    for coefficient in reversed(a):
        out = out * x + coefficient
    return out


def derivative(a):
    return trim(tuple(i * a[i] for i in range(1, len(a))) or (F(0),))


def divmod_poly(a, b):
    a = list(trim(a)); b = trim(b)
    if b == (0,):
        raise ZeroDivisionError
    q = [F(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        degree = len(a) - len(b)
        coefficient = a[-1] / b[-1]
        q[degree] = coefficient
        for j, value in enumerate(b):
            a[degree + j] -= coefficient * value
        a = list(trim(a))
    return trim(q), trim(a)


def gcd_poly(a, b):
    a, b = trim(a), trim(b)
    while b != (0,):
        _, r = divmod_poly(a, b)
        a, b = b, r
    return scale(a, 1 / a[-1])


P = (F(0), F(1))
ONE_MINUS_P = (F(1), F(-1))


def bernoulli_weight(n, k):
    return mul(power(P, k), power(ONE_MINUS_P, n - k))


def mean_polynomial(values, n):
    out = (F(0),)
    for mask, value in enumerate(values):
        out = add(out, scale(bernoulli_weight(n, mask.bit_count()), value))
    return out


def rotate_orbit(torus, mask, z, pair):
    n = torus.n
    key = lambda x, y: (
        (torus.a * x + torus.b * y) % n,
        (-torus.b * x + torus.a * y) % n,
    )
    index = {key(x, y): i for i, (x, y) in enumerate(torus.reps)}

    def rotate_vertex(v):
        x, y = torus.reps[v]
        return index[key(-y, x)]

    def rotate_mask(value):
        output = 0
        for v in range(n):
            if (value >> v) & 1:
                output |= 1 << rotate_vertex(v)
        return output

    orbit = []
    for _ in range(4):
        pair = tuple(sorted(pair))
        orbit.append((mask, z, pair))
        mask, z, pair = (
            rotate_mask(mask),
            rotate_vertex(z),
            tuple(rotate_vertex(v) for v in pair),
        )
    if len(set(orbit)) != 4:
        raise AssertionError("the selected landing does not have a four-point C4 orbit")
    return orbit


def _matrix_entry(torus, orbit, observable, mean_observable):
    n = torus.n
    kernel = (F(0),)
    readout = (F(0),)
    for mask, z, pair in orbit:
        pair_index = torus.pairs.index(pair)
        source = [F(value, 16) for value in (
            torus.pair_values(state)[pair_index] for state in range(1 << n)
        )]
        mean_source = mean_polynomial(source, n)
        filled = mask | (1 << z)
        conditional_weight = bernoulli_weight(n - 1, mask.bit_count())
        source_mid = F(source[mask] + source[filled], 2)
        observable_mid = F(observable[mask] + observable[filled], 2)
        kernel = add(
            kernel,
            mul(
                conditional_weight,
                scale(
                    sub((observable_mid,), mean_observable),
                    source[filled] - source[mask],
                ),
            ),
        )
        readout = add(
            readout,
            mul(
                conditional_weight,
                scale(
                    sub((source_mid,), mean_source),
                    observable[filled] - observable[mask],
                ),
            ),
        )
    return kernel, readout


def exact_landing_minor():
    torus = gate.Torus(3, 0)
    n = torus.n
    # Two occupied sites are a branch-free path.  Filling z closes the unique
    # length-three primitive horizontal cycle.  The selected vacant pair is a
    # physical canonical-kernel port pair whose reconnection is nonzero.
    orbit = rotate_orbit(torus, mask=3, z=3, pair=(2, 6))
    ranks = [torus.rank(mask) for mask in range(1 << n)]
    q = [rank - 1 for rank in ranks]
    energy = [value * value for value in q]
    mean_q = mean_polynomial(q, n)
    mean_energy = mean_polynomial(energy, n)

    landing_rows = []
    for mask, z, pair in orbit:
        filled = mask | (1 << z)
        occupancy = tuple((mask >> v) & 1 for v in torus.nb[z])
        pair_index = torus.pairs.index(pair)
        before = F(torus.pair_values(mask)[pair_index], 16)
        after = F(torus.pair_values(filled)[pair_index], 16)
        degrees_before = sorted(
            sum((mask >> w) & 1 for w in torus.nb[v])
            for v in range(n) if (mask >> v) & 1
        )
        degrees_after = sorted(
            sum((filled >> w) & 1 for w in torus.nb[v])
            for v in range(n) if (filled >> v) & 1
        )
        if occupancy not in ((1, 0, 1, 0), (0, 1, 0, 1)):
            raise AssertionError("landing ports are not alternating")
        if (mask.bit_count(), filled.bit_count()) != (2, 3):
            raise AssertionError("the witness contains an extra occupied branch")
        if degrees_before != [1, 1] or degrees_after != [2, 2, 2]:
            raise AssertionError("the branch-free path/cycle geometry drifted")
        if (ranks[mask], ranks[filled]) != (0, 1):
            raise AssertionError("the selected flip is not the ordinary first-rank pivotal")
        if (before, after) != (F(0), F(1, 4)):
            raise AssertionError("canonical kernel reconnection drifted")
        landing_rows.append({
            "mask": mask, "site": z, "pair": list(pair),
            "neighbor_occupancy_NESW": list(occupancy),
            "rank_before_after": [ranks[mask], ranks[filled]],
            "kernel_before_after": [str(before), str(after)],
        })

    matrix = [
        _matrix_entry(torus, orbit, q, mean_q),
        _matrix_entry(torus, orbit, energy, mean_energy),
    ]
    determinant = sub(mul(matrix[0][0], matrix[1][1]), mul(matrix[0][1], matrix[1][0]))

    factor_a = (F(9), F(0), F(-18), F(9), F(1))
    factor_b = (
        F(1), F(0), F(-4), F(-4), F(16), F(96), F(-332), F(404), F(-224), F(48)
    )
    expected = mul(mul(power(P, 9), power(ONE_MINUS_P, 12)), mul(factor_a, factor_b))
    if determinant != expected:
        raise AssertionError("landing determinant factorization drifted")

    root_polynomial = trim(mean_q)
    common_factor = gcd_poly(root_polynomial, determinant)
    if common_factor != (F(1),):
        raise AssertionError("the landing minor vanishes on a matching-mean root")
    root_lo, root_hi = F(1173, 2000), F(2933, 5000)
    if not (evaluate(root_polynomial, root_lo) < 0 < evaluate(root_polynomial, root_hi)):
        raise AssertionError("the exact finite root bracket drifted")

    # The thermal Schur row operation E -> E-Rq has determinant one.  Report R
    # numerically only as a locator; nonvanishing is certified by polynomial gcd.
    for _ in range(80):
        midpoint = (root_lo + root_hi) / 2
        if evaluate(root_polynomial, midpoint) < 0:
            root_lo = midpoint
        else:
            root_hi = midpoint
    root_mid = (root_lo + root_hi) / 2
    thermal_ratio = evaluate(derivative(mean_energy), root_mid) / evaluate(derivative(mean_q), root_mid)

    return {
        "schema": "matching-one/p537-single-geometry-kernel-readout-minor/v2",
        "geometry": {"period": [3, 0], "N": n},
        "landing_definition": {
            "occupied_geometry": "two-site branch-free path closed by the flip into one primitive length-three cycle",
            "port_condition": "alternating occupied/vacant NESW incident ports",
            "angular_sum": "four C4 rotations with spin-4 character +1",
            "source": "canonical Kreg pair kernel g16/16",
        },
        "orbit": landing_rows,
        "matrix_rows": ["q", "E=q^2"],
        "matrix_columns": ["kernel_reconnection", "readout_pivotal"],
        "matrix_polynomials_ascending": [
            [[str(x) for x in entry] for entry in row] for row in matrix
        ],
        "root_polynomial_ascending": [str(x) for x in root_polynomial],
        "root_polynomial_display": "-4*p^9+18*p^8-18*p^7+6*p^3-1",
        "root_bracket": [str(F(1173, 2000)), str(F(2933, 5000))],
        "root_decimal": float(root_mid),
        "thermal_ratio_R_decimal": float(thermal_ratio),
        "canonical_determinant_factorization": (
            "p^9*(1-p)^12*(p^4+9*p^3-18*p^2+9)*"
            "(48*p^9-224*p^8+404*p^7-332*p^6+96*p^5+16*p^4-4*p^3-4*p^2+1)"
        ),
        "determinant_at_half": str(evaluate(determinant, F(1, 2))),
        "determinant_at_root_decimal": float(evaluate(determinant, root_mid)),
        "gcd_root_and_determinant": [str(x) for x in common_factor],
        "schur_projection": {
            "row_operation": "E -> E-R*q, R=(d_p mean E)/(d_p mean q) at the finite root",
            "determinant_invariant": True,
        },
        "decision": "nonzero_single_geometry_C4_kernel_readout_minor",
        "scope": (
            "exact supporting control only: quarter-turn C4 is not the original "
            "axis-versus-tilted P4 projector, and kernel/readout are two terms of "
            "one covariance derivative rather than independent source/thermal columns"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = exact_landing_minor()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "decision": result["decision"],
        "root_decimal": result["root_decimal"],
        "determinant_at_root_decimal": result["determinant_at_root_decimal"],
    }, indent=2))


if __name__ == "__main__":
    main()
