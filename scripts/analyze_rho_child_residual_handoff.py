#!/usr/bin/env python3
"""Post-reveal decomposition of the frozen E4+E6 rho-child residual."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from derive_hexagonal_degree2_hecke import normalized_eisenstein
from score_rho_child_primitive_h4 import CHILD_DESIGNS, tau_from_matrix


N30_CHILDREN = (
    ("2omega", ((6, 6), (0, 10))),
    ("omega_over_2", ((12, 3), (0, 5))),
    ("omega_plus_1_over_2", ((12, 9), (0, 5))),
)


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def cross(left, right):
    return [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]


def modular_columns(designs, *, dps):
    e4, e6 = [], []
    for _, matrix in designs:
        tau = tau_from_matrix(matrix)
        e4.append(normalized_eisenstein(4, tau, dps=dps)[0])
        e6.append(normalized_eisenstein(6, tau, dps=dps)[0])
    return e4, e6, [value * value for value in e4]


def decompose(score, *, dps=60):
    mp.mp.dps = dps
    model = score["two_character_mixture_opponents"]["E4_r1+E6_r0"]
    residual = [
        mp.mpc(model["residual"][2 * j], model["residual"][2 * j + 1])
        for j in range(3)
    ]
    covariance = mp.matrix(
        [[mp.mpf(value) for value in row] for row in score["full_covariance_6x6"]]
    )
    child_contributions = []
    for j, name in enumerate(score["child_order"]):
        block = mp.matrix([[covariance[2*j+i, 2*j+k] for k in range(2)] for i in range(2)])
        vector = mp.matrix([mp.re(residual[j]), mp.im(residual[j])])
        chi2 = (vector.T * block**-1 * vector)[0]
        child_contributions.append({
            "child": name,
            "residual_re_im": [mp.nstr(mp.re(residual[j]), 20), mp.nstr(mp.im(residual[j]), 20)],
            "residual_phase_rad": mp.nstr(mp.arg(residual[j]), 20),
            "chi_square_contribution": mp.nstr(chi2, 20),
        })
    zeta = mp.exp(2 * mp.pi * mp.j / 3)
    dft = [sum(residual[j] * zeta**(-r*j) for j in range(3)) / 3 for r in range(3)]

    annihilators = {}
    current_coefficient = None
    delta = [mp.mpc(score["delta_H4_re_im"][2*j], score["delta_H4_re_im"][2*j+1]) for j in range(3)]
    for label, designs in (("N112_Dplus1", CHILD_DESIGNS), ("N60_Dminus2", N30_CHILDREN)):
        e4, e6, e4sq = modular_columns(designs, dps=dps)
        null = cross(e4, e6)
        denominator = dot(null, e4sq)
        row = {
            "matrices": [matrix for _, matrix in designs],
            "w_dot_E4_abs": mp.nstr(abs(dot(null, e4)), 12),
            "w_dot_E6_abs": mp.nstr(abs(dot(null, e6)), 12),
            "w_dot_E4sq_re_im": [mp.nstr(mp.re(denominator), 20), mp.nstr(mp.im(denominator), 20)],
            "w_dot_E4sq_phase_rad": mp.nstr(mp.arg(denominator), 20),
        }
        if label == "N112_Dplus1":
            current_coefficient = dot(null, delta) / denominator
            row["inferred_E4sq_completion_re_im"] = [
                mp.nstr(mp.re(current_coefficient), 20), mp.nstr(mp.im(current_coefficient), 20)
            ]
            row["inferred_E4sq_completion_phase_rad"] = mp.nstr(mp.arg(current_coefficient), 20)
        annihilators[label] = row

    return {
        "schema": "matching-one/rho-child-residual-handoff/v1",
        "status": "post_reveal_discovery_not_a_rescore",
        "source_score_status": score["status"],
        "frozen_E4_plus_E6": {
            "chi_square": model["chi_square"],
            "dof": model["dof"],
            "survival_p": model["survival_p"],
            "child_contributions": child_contributions,
            "ideal_DFT_of_residual": [
                {
                    "r": r,
                    "value_re_im": [mp.nstr(mp.re(value), 20), mp.nstr(mp.im(value), 20)],
                    "abs": mp.nstr(abs(value), 20),
                    "phase_rad": mp.nstr(mp.arg(value), 20),
                }
                for r, value in enumerate(dft)
            ],
            "r2_to_largest_other_DFT_magnitude": mp.nstr(abs(dft[2]) / max(abs(dft[0]), abs(dft[1])), 20),
        },
        "next_design": {
            "parent": "opposite-side Pell Dminus2 N30",
            "children": "three degree-2 N60 children",
            "reason": "the current failed complex relation is a coherent r2 residual; the opposite Pell side tests its phase while reversing signed Pell leakage",
            "primary_normalization_free_quantity": "phase of (w dot delta_H4)/(w dot E4_squared), where w dot E4=w dot E6=0 within each triple",
            "annihilator_geometry": annihilators,
            "production_authorized": False,
        },
    }


def render(result):
    rows = result["frozen_E4_plus_E6"]["child_contributions"]
    dft = result["frozen_E4_plus_E6"]["ideal_DFT_of_residual"]
    return "\n".join([
        "# Post-reveal rho-child residual handoff", "",
        "The frozen E4+E6 relation fails coherently, not at one bad child.", "",
        "| child | residual (re, im) | chi-square contribution |", "|---|---:|---:|",
        *[f"| {row['child']} | {row['residual_re_im']} | {row['chi_square_contribution']} |" for row in rows],
        "", "The residual DFT is", "",
        *[f"- r={row['r']}: `{row['value_re_im']}`, magnitude `{row['abs']}`" for row in dft],
        "",
        f"The r2 magnitude is `{result['frozen_E4_plus_E6']['r2_to_largest_other_DFT_magnitude']}` times the largest other residual row.",
        "Thus the frozen two-character function misses one coherent r2-shaped complex relation; this is distinct from the raw-data statement that standalone r2 remains unresolved.",
        "", "## Next recognition experiment", "",
        "Use the three degree-2 children of the opposite-side Pell N30 parent:", "",
        "- `[[6,6],[0,10]]`", "- `[[12,3],[0,5]]`", "- `[[12,9],[0,5]]`", "",
        "Within each triple construct the exact complex annihilator `w` of E4 and E6, then compare the phase of `(w dot delta_H4)/(w dot E4_squared)` between N60 and N112. This removes overall normalization. E4-squared completion predicts phase preservation; signed Pell leakage predicts reversal or loss of coherence. Production is not authorized by this post-reveal note.", "",
    ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("score", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = decompose(json.loads(args.score.read_text()))
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
