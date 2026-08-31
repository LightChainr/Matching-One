#!/usr/bin/env python3
"""P398: exact positive-weight width-four A/landing propagation, not another jet."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

from noncrossing_connectivity_codec import noncrossing_states
from p321_homology_trace_certificate import action_matrix, identity, join_adjacent, rotate_state
from p333_generic_q_detach_intertwiner import detach_jet
from p333_gram_source_intertwiner import multiply, rref_solve, transpose
from p333_source_landing_doublet_width4 import landing_gram_jet
from p398_rooted_gr1_completion import selected_completion_families

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/p398_positive_cylinder_protocol.json"
RESULT = ROOT / "results/p398-positive-cylinder/latest.json"


def encode(value):
    if isinstance(value, F):
        return int(value) if value.denominator == 1 else str(value)
    if isinstance(value, (list, tuple)):
        return [encode(x) for x in value]
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    return value


def zadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def zmul(a, b):
    return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])


def zconj(a):
    return (a[0], -a[1])


def zdet(a):
    x, y = zmul(a[0][0], a[1][1]), zmul(a[0][1], a[1][0])
    return (x[0]-y[0], x[1]-y[1])


def zmatmul(a, b):
    return [[zadd(zmul(a[i][0], b[0][j]), zmul(a[i][1], b[1][j]))
             for j in range(2)] for i in range(2)]


def bond_factor(generator):
    size = len(generator)
    return [[(F(generator[r][c]) + int(r == c))/2 for c in range(size)] for r in range(size)]


def build_model():
    states = noncrossing_states(4)
    size = len(states)
    horizontal, vertical = identity(size), identity(size)
    for site in range(4):
        join = action_matrix(4, lambda state, site=site: join_adjacent(state, site))
        detach, _ = detach_jet(4, site)
        horizontal = multiply(bond_factor(join), horizontal)
        vertical = multiply(bond_factor(detach), vertical)
    transfer = multiply(horizontal, vertical)
    equations = [[transfer[r][c]-int(r == c) for c in range(size)] for r in range(size)]
    stationary_solution = rref_solve(equations + [[1]*size], [0]*size+[1], size)
    assert stationary_solution["consistent"] and stationary_solution["dimension"] == 0
    stationary = stationary_solution["particular"]
    assert all(p > 0 for p in stationary)
    assert all(sum(transfer[r][c] for r in range(size)) == 1 for c in range(size))
    assert all(value >= 0 for row in transfer for value in row)
    ap = selected_completion_families()["rooted_charge1"]["columns"]
    gram0, _ = landing_gram_jet()
    readouts = [ap[r] + gram0[r][-2:] for r in range(size)]
    # All four real columns have exact zero stationary mean by C4 covariance.
    assert all(sum(stationary[r]*readouts[r][c] for r in range(size)) == 0 for c in range(4))
    rotation = action_matrix(4, lambda state: rotate_state(state, 1))
    assert multiply(transfer, rotation) == multiply(rotation, transfer)
    rotation_squared = multiply(rotation, rotation)
    odd_space = rref_solve([[rotation_squared[r][c]+int(r == c) for c in range(size)]
                            for r in range(size)], [0]*size, size)
    assert odd_space["dimension"] == 4  # complete +/-i sector, not a fitted subspace
    return states, transfer, stationary, readouts


def sector_matrix(transfer, readouts):
    propagated = multiply(transpose(transfer), readouts)
    columns = []
    for col in range(4):
        solution = rref_solve(readouts, [row[col] for row in propagated], 4)
        assert solution["consistent"] and solution["dimension"] == 0
        columns.append(solution["particular"])
    real = transpose(columns)
    complex_matrix = [[(real[2*a][2*b], real[2*a][2*b+1]) for b in range(2)] for a in range(2)]
    for a in range(2):
        for b in range(2):
            re, im = complex_matrix[a][b]
            assert real[2*a+1][2*b] == -im and real[2*a+1][2*b+1] == re
    return real, complex_matrix


def correlation(stationary, readouts, propagated):
    real = [[sum(stationary[r]*readouts[r][a]*propagated[r][b]
                 for r in range(len(stationary))) for b in range(4)] for a in range(4)]
    return [[(real[2*a][2*b]+real[2*a+1][2*b+1],
              real[2*a+1][2*b]-real[2*a][2*b+1]) for b in range(2)] for a in range(2)]


def result():
    states, transfer, stationary, readouts = build_model()
    sector_real, sector = sector_matrix(transfer, readouts)
    trace = zadd(sector[0][0], sector[1][1])
    determinant = zdet(sector)
    assert trace[1] == determinant[1] == 0
    discriminant = trace[0]**2 - 4*determinant[0]
    assert discriminant > 0
    eigenvalues = [(float(trace[0])+s*math.sqrt(float(discriminant)))/2 for s in (1,-1)]
    c0 = correlation(stationary, readouts, readouts)
    rows, propagated, previous = [], readouts, None
    for distance in range(9):
        value = correlation(stationary, readouts, propagated)
        if previous is not None:
            assert value == zmatmul(previous, [[zconj(z) for z in row] for row in sector])
        det = zdet(value)
        expected_det = zmul(zdet(c0), (determinant[0]**distance, F(0)))
        assert det == expected_det
        rows.append({"d": distance, "matrix_complex_re_im": value,
                     "matrix_decimal_re_im": [[[float(v) for v in z] for z in row] for row in value],
                     "determinant_re_im": det,
                     "determinant_decimal": float(det[0]),
                     "rank_over_C": 2 if det != (0,0) else 1})
        previous = value
        propagated = multiply(transpose(transfer), propagated)
    cross_asymmetry = zadd(rows[1]["matrix_complex_re_im"][0][1],
                          tuple(-x for x in zconj(rows[1]["matrix_complex_re_im"][1][0])))
    input_paths = [PROTOCOL, Path(__file__),
                   *(ROOT / "scripts" / name for name in (
                       "noncrossing_connectivity_codec.py", "p321_homology_trace_certificate.py",
                       "p333_generic_q_detach_intertwiner.py", "p333_gram_source_intertwiner.py",
                       "p333_source_landing_doublet_width4.py", "p398_rooted_gr1_completion.py"))]
    return encode({
        "schema": "matching-one/p398-positive-cylinder-result/v1",
        "protocol_sha256": hashlib.sha256(PROTOCOL.read_bytes()).hexdigest(),
        "input_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in input_paths},
        "arithmetic": "fractions.Fraction; decimals are display only",
        "state_order": states,
        "stationary_probability": stationary,
        "transfer_column_source": transfer,
        "readouts_Ar_Ai_Lr_Li": readouts,
        "backward_sector_real": sector_real,
        "backward_sector_complex_re_im": sector,
        "complete_charge_one_complex_dimension": 2,
        "characteristic_polynomial_high_to_low": [F(1), -trace[0], determinant[0]],
        "discriminant": discriminant,
        "eigenvalues_decimal": eigenvalues,
        "eigenvalues_exact": ["(3+sqrt(5))/64", "(3-sqrt(5))/64"],
        "eigenfunctions_exact": ["A+(1-sqrt(5))*(1-i)*L/4", "A+(1+sqrt(5))*(1-i)*L/4"],
        "decay_lengths_in_rows": [-1/math.log(x) for x in eigenvalues],
        "subleading_to_leading_decay_ratio": "(7-3*sqrt(5))/2",
        "matrix_recurrence": "C(d+2)=(3/32)*C(d+1)-(1/1024)*C(d)",
        "determinant_law": "det C(d)=(73216/1940449)*(1/1024)^d",
        "d1_cross_asymmetry_C_AL_minus_conjugate_C_LA": cross_asymmetry,
        "correlations": rows,
        "decision": "two_distinct_nonzero_charge_one_transfer_eigenvalues; no_exact_common_rank_one_response_ray",
        "boundary": "width-four positive square-bond frontier response only; not site Matching, continuum or physical Jordan identification"
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=RESULT)
    args = parser.parse_args()
    value = result()
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(value, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: value[k] for k in ("backward_sector_complex_re_im", "characteristic_polynomial_high_to_low", "eigenvalues_decimal", "decay_lengths_in_rows", "decision")}, indent=2))
    print(json.dumps(value["correlations"][1], indent=2))


if __name__ == "__main__":
    main()
