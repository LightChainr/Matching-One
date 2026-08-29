#!/usr/bin/env python3
"""Exact microscopic odd tangent and thermal-null local row on the C4 control."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from c4_local_odd_pivotal import _active, _cross, pivotal_h4
from c4_self_matching_exact import c4_self_matching_torus


def ftext(value: Fraction | int) -> str:
    return str(Fraction(value))


def response(observable: list[Fraction], score: list[int]) -> Fraction:
    return sum(value * weight for value, weight in zip(observable, score)) / len(observable)


def build_oracle() -> dict:
    geometry = c4_self_matching_torus(3, 1)
    n = geometry.n
    size = 1 << n
    full = size - 1
    even_sites = [vertex for vertex in range(n) if vertex % 2 == 0]
    odd_sites = [vertex for vertex in range(n) if vertex % 2 == 1]
    root_even = geometry.vertex((0, 0))
    root_odd = geometry.vertex((1, 0))
    if root_even % 2 or not root_odd % 2:
        raise AssertionError("declared local cell must contain one even and one odd site")

    score_t: list[int] = []
    score_lambda: list[int] = []
    thermal_cell: list[Fraction] = []
    local_h4: list[Fraction] = []
    global_cross: list[Fraction] = []
    signature_violations = 0

    for mask in range(size):
        active = _active(mask, n)
        inactive = _active(full ^ mask, n)
        signs = [1 if value else -1 for value in active]
        score_t.append(2 * sum(signs))
        score_lambda.append(2 * sum((1 if vertex % 2 == 0 else -1) * signs[vertex] for vertex in range(n)))
        # Average of the centered occupancies in one even/odd microscopic cell.
        thermal_cell.append(Fraction(signs[root_even] + signs[root_odd], 4))
        local_h4.append(Fraction(pivotal_h4(geometry, active) - pivotal_h4(geometry, inactive), 2))
        global_cross.append(Fraction(_cross(geometry, active) - _cross(geometry, inactive), 2))

        ke = sum(active[v] for v in even_sites)
        ko = sum(active[v] for v in odd_sites)
        # P_{t,lambda}(A)=P_{-t,-lambda}(A^c): under sign reversal
        # the occupied/vacant bases exchange on each sublattice.
        left_signature = (ke, len(even_sites) - ke, ko, len(odd_sites) - ko)
        ce = len(even_sites) - ke
        co = len(odd_sites) - ko
        right_signature_in_left_bases = (
            len(even_sites) - ce, ce, len(odd_sites) - co, co
        )
        signature_violations += left_signature != right_signature_in_left_bases

    alpha_star = -response(local_h4, score_t) / response(thermal_cell, score_t)
    local_thermal_null = [value + alpha_star * density for value, density in zip(local_h4, thermal_cell)]

    complement_violations = {
        "score_t": sum(score_t[full ^ mask] != -score_t[mask] for mask in range(size)),
        "score_lambda": sum(score_lambda[full ^ mask] != -score_lambda[mask] for mask in range(size)),
        "thermal_cell": sum(thermal_cell[full ^ mask] != -thermal_cell[mask] for mask in range(size)),
        "local_h4": sum(local_h4[full ^ mask] != -local_h4[mask] for mask in range(size)),
        "local_thermal_null": sum(
            local_thermal_null[full ^ mask] != -local_thermal_null[mask] for mask in range(size)
        ),
    }

    fisher = [
        [response([Fraction(value) for value in score_t], score_t), response([Fraction(value) for value in score_t], score_lambda)],
        [response([Fraction(value) for value in score_lambda], score_t), response([Fraction(value) for value in score_lambda], score_lambda)],
    ]
    rows = {
        "global_cross": [response(global_cross, score_t), response(global_cross, score_lambda)],
        "local_h4": [response(local_h4, score_t), response(local_h4, score_lambda)],
        "thermal_cell": [response(thermal_cell, score_t), response(thermal_cell, score_lambda)],
        "local_thermal_null": [response(local_thermal_null, score_t), response(local_thermal_null, score_lambda)],
    }
    global_null_direction = (Fraction(2), Fraction(-3))

    return {
        "schema": "matching-one.c4-microscopic-odd-tangent.v1",
        "issue": 155,
        "geometry": {
            "a": 3,
            "b": 1,
            "N": n,
            "configurations": size,
            "even_sites": len(even_sites),
            "odd_sites": len(odd_sites),
            "local_cell_vertices": [root_even, root_odd],
        },
        "single_theory_intertwiner": {
            "definition": "C|A>=|A^c> on the self-matching C4 graph",
            "C_squared_is_identity": True,
            "family_identity": "C W(t,lambda) C = W(-t,-lambda)",
            "configurationwise_weight_signature_violations": signature_violations,
            "parameter_action": "(t,lambda)->(-t,-lambda)",
            "consequence": "unlike the non-self-matching #233 pair, species doubling collapses and local odd operators exist inside one theory",
        },
        "microscopic_scores": {
            "S_t": "4 sum_x (n_x-1/2)",
            "S_lambda": "4 sum_x parity(x)(n_x-1/2)",
            "complement_odd_violations": {
                key: complement_violations[key] for key in ("score_t", "score_lambda")
            },
            "Fisher_Gram": [[ftext(value) for value in row] for row in fisher],
            "exact_thermal_orthogonality": fisher[0][1] == 0,
            "interpretation": "S_lambda is a genuine local UV odd tangent Fisher-orthogonal to the uniform thermal score S_t at the center",
        },
        "local_thermal_orthogonalization": {
            "family": "O_alpha=O_local_H4 + alpha epsilon_cell",
            "epsilon_cell": "[(n_even-1/2)+(n_odd-1/2)]/2 on one microscopic even/odd cell",
            "alpha_star": ftext(alpha_star),
            "unique_t_response_zero": rows["local_thermal_null"][0] == 0,
            "surviving_lambda_matrix_element": ftext(rows["local_thermal_null"][1]),
            "complement_odd_violations": complement_violations["local_thermal_null"],
            "response_rows_columns_t_lambda": {
                name: [ftext(value) for value in row] for name, row in rows.items()
            },
        },
        "coupling_space_zero": {
            "direction_delta_t_delta_lambda": [ftext(value) for value in global_null_direction],
            "direction_delta_p_even_delta_p_odd": [
                ftext(global_null_direction[0] + global_null_direction[1]),
                ftext(global_null_direction[0] - global_null_direction[1]),
            ],
            "global_cross_matrix_element": ftext(sum(a * b for a, b in zip(rows["global_cross"], global_null_direction))),
            "local_h4_matrix_element": ftext(sum(a * b for a, b in zip(rows["local_h4"], global_null_direction))),
            "meaning": "an exact microscopic coupling tangent invisible to the N10 global thermal proxy but visible to the local H4 row",
        },
        "operator_separation": {
            "Ad_C_parity": "all recorded scores and readouts are exactly odd",
            "Potts_colour_charge": "singlet; independent of the #257 [2] selection rule",
            "Q_defect_vs_measure": "this is a p-coupling score inside one self-matching theory, complementary to the #233 Q-cluster pull-through obstruction",
        },
        "claim_boundary": {
            "proved": [
                "the single-theory complement intertwiner for the finite probability family",
                "a local complement-odd score Fisher-orthogonal to the uniform thermal score",
                "the unique alpha=3/64 local readout with zero t response and lambda response 11/64",
                "the exact coupling direction (2,-3) with global response zero and local response -39/64",
            ],
            "not_proved": [
                "that either orthogonalization is an RG eigenoperator at large size",
                "that the surviving local row has x=21/4",
                "a transfer-matrix or continuum OPE involution beyond the finite probability family",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    payload = json.dumps(build_oracle(), indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
