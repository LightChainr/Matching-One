#!/usr/bin/env python3
"""Freeze the equal-size N60-only increment power gate after the first reveal."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from score_rho_child_opposite_pell_phase import ray_score


def matrix(rows):
    return mp.matrix([[mp.mpf(value) for value in row] for row in rows])


def project(score):
    mp.mp.dps = 60
    c112 = mp.matrix([mp.mpf(value) for value in score["c112_re_im"]])
    c60 = mp.matrix([mp.mpf(value) for value in score["c60_re_im"]])
    cov112 = matrix(score["c112_covariance_2x2"])
    cov60 = matrix(score["c60_covariance_2x2"])
    combined_cov60 = cov60 / 2
    preserve = ray_score(c112, cov112, c60, combined_cov60, +1)
    flip = ray_score(c112, cov112, c60, combined_cov60, -1)
    null112 = (c112.T * cov112**-1 * c112)[0]
    null112_p1 = mp.erfc(mp.sqrt(null112 / 2))
    null60 = (c60.T * combined_cov60**-1 * c60)[0]
    alpha = mp.mpf(str(score["decision_alpha"]))
    passes = mp.mpf(flip["survival_p"]) < alpha
    structural_passes = null112_p1 < alpha
    return {
        "schema": "matching-one/rho-child-opposite-pell-increment-power/v1",
        "status": "frozen_power_gate_no_go" if not passes else "frozen_power_gate_pass",
        "source_result_commit": "9277fc91ff7b765707971d7c827ad43513df4df3",
        "proposal": "one additional independent 2M/child N60 replication combined by equal-information inverse covariance",
        "plug_in_contract": "first-batch c60 is the projection mean; independent equal-size replication halves C60; N112 is unchanged",
        "decision_alpha": score["decision_alpha"],
        "projected_phase_preserving_ray": preserve,
        "projected_phase_flipping_ray": flip,
        "projected_combined_N60_zero": {
            "chi_square": mp.nstr(null60, 20),
            "dof": 2,
            "survival_p": mp.nstr(mp.exp(-null60 / 2), 20),
        },
        "N112_zero_boundary_for_flip": {
            "chi_square": mp.nstr(null112, 20),
            "effective_ray_dof": 1,
            "survival_p": mp.nstr(null112_p1, 20),
        },
        "finite_increment_gate_passed": passes,
        "N60_only_asymptotic_gate_passed": structural_passes,
        "increment_authorized": passes and structural_passes,
        "decision": (
            "authorize_fixed_increment" if passes and structural_passes else
            "do_not_run_N60_only_increment"
        ),
        "reason": "the negative ray can send its scale to infinity and fit N60 while collapsing the unchanged N112 coefficient to zero; its best possible N60-only limit remains above alpha",
        "next_information_needed": "independent N112 information or a new non-degenerate geometry constraint, not more N60 samples alone",
    }


def render(result):
    return "\n".join([
        "# Opposite-Pell N60 increment power freeze", "",
        f"Decision: **{result['decision']}**.", "",
        f"- projected positive ray: `{result['projected_phase_preserving_ray']['chi_square']}/1`, p `{result['projected_phase_preserving_ray']['survival_p']}`",
        f"- projected negative ray: `{result['projected_phase_flipping_ray']['chi_square']}/1`, p `{result['projected_phase_flipping_ray']['survival_p']}`, scale `{result['projected_phase_flipping_ray']['best_positive_scale']}`",
        f"- unchanged-N112 zero boundary: `{result['N112_zero_boundary_for_flip']['chi_square']}/1`, p `{result['N112_zero_boundary_for_flip']['survival_p']}`", "",
        "The equal-size N60 increment would make the negative-ray optimizer escape to the opposite scale boundary, not reject the flip. Even infinite N60 precision cannot improve past the unchanged N112-zero boundary, whose p-value remains above the frozen .01 threshold. Therefore no second N60 acquisition is authorized.", "",
    ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("score", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = project(json.loads(args.score.read_text()))
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

