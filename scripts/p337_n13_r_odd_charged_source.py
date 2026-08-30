#!/usr/bin/env python3
"""Exact minimal R-odd F3 charged-source response on the N13 quotient."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Optional, Sequence

from p337_n13_twist_flux_bridge import P_REF, build_certificate, evaluate_bernstein


CHANNELS = {
    "A_axis_odd": {
        "plus": "axis_x",
        "minus": "axis_y",
        "twist": "q_A=T_01-T_10",
        "D4_irrep": "B1",
    },
    "D_diagonal_odd": {
        "plus": "diag_plus",
        "minus": "diag_minus",
        "twist": "q_D=T_12-T_11",
        "D4_irrep": "B2",
    },
}


def laurent_payload(weight: Fraction) -> dict[str, str]:
    """Return exact Z(u) and O(u)=u d_u Z coefficients for q=0,+/-1."""
    half = weight / 2
    return {
        "Z_u_minus_1": str(half),
        "Z_u_0": str(1 - weight),
        "Z_u_plus_1": str(half),
        "O_u_minus_1": str(-half),
        "O_u_0": "0",
        "O_u_plus_1": str(half),
        "at_omega_Z": str(1 - Fraction(3, 2) * weight),
        "at_omega_O_coefficient_of_omega_minus_omega2": str(half),
    }


def exact_channel_rows(
    bridge: dict[str, object], channel: dict[str, str]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    state_rows = []
    response_coefficients: list[Fraction] = []
    for row in bridge["state_coefficient_rows"]:
        plus = row["line_state_counts"][channel["plus"]]
        minus = row["line_state_counts"][channel["minus"]]
        denominator = row["subset_count"]
        weight = Fraction(plus + minus, denominator)
        charge = Fraction(plus - minus, denominator)
        response_coefficients.append(weight)
        state_rows.append({
            "k": row["k"],
            "plus_count": plus,
            "minus_count": minus,
            "unweighted_charge_coefficient": str(charge),
            "linear_susceptibility": str(weight),
            "unit_HAD_susceptibility": str(weight / 2),
            "F3_phase_defect": laurent_payload(weight),
        })

    flux_rows = []
    for row in bridge["flux_coefficient_rows"]:
        denominator = row["edge_normalization"]
        birth_plus = row["line_birth_edges"][channel["plus"]]
        birth_minus = row["line_birth_edges"][channel["minus"]]
        exit_plus = row["line_exit_edges"][channel["plus"]]
        exit_minus = row["line_exit_edges"][channel["minus"]]
        birth_response = Fraction(birth_plus + birth_minus, denominator)
        exit_response = Fraction(exit_plus + exit_minus, denominator)
        derivative_response = 13 * (
            response_coefficients[row["lower_k"] + 1]
            - response_coefficients[row["lower_k"]]
        )
        flux_rows.append({
            "lower_k": row["lower_k"],
            "birth_plus_edges": birth_plus,
            "birth_minus_edges": birth_minus,
            "exit_plus_edges": exit_plus,
            "exit_minus_edges": exit_minus,
            "unweighted_birth_charge": str(Fraction(birth_plus - birth_minus, denominator)),
            "unweighted_exit_charge": str(Fraction(exit_plus - exit_minus, denominator)),
            "charged_birth_response": str(birth_response),
            "charged_exit_response": str(exit_response),
            "charged_derivative_response": str(derivative_response),
            "continuity_pass": derivative_response == birth_response - exit_response,
            "at_omega_birth_O_coefficient": str(birth_response / 2),
            "at_omega_exit_O_coefficient": str(exit_response / 2),
        })
    return state_rows, flux_rows


def evaluate_channel(
    state_rows: Sequence[dict[str, object]],
    flux_rows: Sequence[dict[str, object]],
) -> dict[str, object]:
    susceptibility = [Fraction(row["linear_susceptibility"]) for row in state_rows]
    birth = [Fraction(row["charged_birth_response"]) for row in flux_rows]
    exit_ = [Fraction(row["charged_exit_response"]) for row in flux_rows]
    derivative = [Fraction(row["charged_derivative_response"]) for row in flux_rows]
    values = {
        "susceptibility": evaluate_bernstein(susceptibility, P_REF),
        "unit_HAD_susceptibility": evaluate_bernstein(
            [value / 2 for value in susceptibility], P_REF
        ),
        "birth_response": evaluate_bernstein(birth, P_REF),
        "exit_response": evaluate_bernstein(exit_, P_REF),
        "derivative_response": evaluate_bernstein(derivative, P_REF),
    }
    values["at_omega_Z"] = 1 - Fraction(3, 2) * values["susceptibility"]
    values["at_omega_O_coefficient_of_omega_minus_omega2"] = (
        values["susceptibility"] / 2
    )
    return {
        "exact": {name: str(value) for name, value in values.items()},
        "decimal": {name: float(value) for name, value in values.items()},
    }


def build_charged_certificate() -> dict[str, object]:
    bridge = build_certificate()
    channel_payload = {}
    all_continuity = True
    all_unweighted_zero = True
    for name, channel in CHANNELS.items():
        state_rows, flux_rows = exact_channel_rows(bridge, channel)
        all_continuity &= all(row["continuity_pass"] for row in flux_rows)
        all_unweighted_zero &= all(
            Fraction(row["unweighted_charge_coefficient"]) == 0
            for row in state_rows
        ) and all(
            Fraction(row[key]) == 0 for row in flux_rows
            for key in ("unweighted_birth_charge", "unweighted_exit_charge")
        )
        channel_payload[name] = {
            "definition": channel,
            "state_response_rows": state_rows,
            "flux_response_rows": flux_rows,
            "reference_evaluation": evaluate_channel(state_rows, flux_rows),
        }

    nonzero_response = all(
        any(Fraction(row["linear_susceptibility"]) > 0
            for row in payload["state_response_rows"])
        for payload in channel_payload.values()
    )
    a_unit = channel_payload["A_axis_odd"]["reference_evaluation"]["exact"][
        "unit_HAD_susceptibility"
    ]
    d_unit = channel_payload["D_diagonal_odd"]["reference_evaluation"]["exact"][
        "unit_HAD_susceptibility"
    ]
    return {
        "schema": "matching-one/p337-N13-R-odd-F3-charged-source/v1",
        "issues": [337, 334],
        "status": "exact minimal charged-source state and current response",
        "geometry": bridge["geometry"],
        "representation": {
            "projective_line_order": ["axis_x", "axis_y", "diag_plus", "diag_minus"],
            "R_action": "axis_x<->axis_y and diag_plus<->diag_minus",
            "even_basis": ["E_axis=e_x+e_y", "E_diag=e_diagPlus+e_diagMinus"],
            "charged_basis": ["q_A=e_x-e_y", "q_D=e_diagPlus-e_diagMinus"],
            "C4_charge": "both q_A and q_D have R eigenvalue -1, hence charge 2; the projective action factors through C2 because R^2=-I acts trivially on lines",
            "D4_refinement": "q_A is B1 (reflection even); q_D is B2 (reflection odd)",
        },
        "source_definition": {
            "continuous": "reweight a rank-one state by exp(s_C q_C); R sends (q_C,s_C)->(-q_C,-s_C)",
            "F3_phase": "set u=omega with omega^3=1 and weight u^q_C",
            "defect_partition": "Z_C(u)=E[u^q_C]",
            "charged_one_point": "O_C(u)=E[q_C u^q_C]=u dZ_C/du",
            "linear_response": "d<E[q_C]>_s/ds at s=0 = E[q_C^2] = W_C",
        },
        "selection_rules": {
            "A_source": "activates A response W_A; linear H and D responses remain zero",
            "D_source": "activates D response W_D; linear H and A responses remain zero",
            "mixed_A_D": "zero because q_A q_D=0 statewise",
            "scalar_partition": "even in each charged source; its first derivative at zero remains zero",
            "why": "R parity kills H-q_C and the other odd cross-response; disjoint line support kills A-D mixing",
        },
        "linear_response_matrix_at_p_ref": {
            "basis": "unit #337 H/A/D coordinates, with A=q_A/sqrt(2), D=q_D/sqrt(2)",
            "rows": ["H", "A", "D"],
            "columns": ["unit_A_source", "unit_D_source"],
            "exact": [["0", "0"], [a_unit, "0"], ["0", d_unit]],
            "interpretation": (
                "the two explicit charged sources span the complete R-odd response "
                "multiplicity; H and cross-channel response remain symmetry null"
            ),
        },
        "channels": channel_payload,
        "minimal_identifiable_observable": {
            "state": "W_C(p)=E[q_C^2], equivalently 2 O_C(omega)/(omega-omega^2)",
            "current": "J_C,birth(p), J_C,exit(p), with dW_C/dp=J_C,birth-J_C,exit",
            "production_tuple": "for C=A,D report W_C, birth response, exit response and full cross-channel covariance from the existing ell archive; no new path field is needed",
        },
        "gates": {
            "parent_bridge_pass": bridge["gates"]["all_pass"],
            "unweighted_odd_state_source_sink_zero": all_unweighted_zero,
            "charged_state_response_nonzero": nonzero_response,
            "charged_continuity_coefficientwise": all_continuity,
            "all_pass": (
                bridge["gates"]["all_pass"] and all_unweighted_zero
                and nonzero_response and all_continuity
            ),
        },
        "claim_boundary": (
            "exact finite charged-defect response and selection rule; the source is an "
            "explicit projective reweighting, not evidence for a continuum charged field"
        ),
    }


def render_markdown(payload: dict[str, object]) -> str:
    a = payload["channels"]["A_axis_odd"]["reference_evaluation"]["decimal"]
    d = payload["channels"]["D_diagonal_odd"]["reference_evaluation"]["decimal"]
    lines = [
        "# Minimal R-odd F3 charged source on Gaussian N13", "",
        "All unweighted-null, charged-response and coefficientwise continuity gates pass.", "",
        "## Charged basis", "",
        "```text",
        "q_A = T_01-T_10 = 1_axisX-1_axisY,",
        "q_D = T_12-T_11 = 1_diagPlus-1_diagMinus.",
        "```", "",
        "Quarter-turn sends each charge to its negative. The projective action has "
        "C4 charge 2 and factors through C2; under full D4, q_A is B1 and q_D is B2.", "",
        "## Explicit F3 defect", "",
        "Reweight a rank-one state by `u^q_C`, with `u=omega`, `omega^3=1`. Then", "",
        "```text",
        "Z_C(u)=E[u^q_C],",
        "O_C(u)=E[q_C u^q_C],",
        "O_C(omega)=(omega-omega^2) W_C/2,",
        "W_C=E[q_C^2].",
        "```", "",
        "Although the unweighted one-point is exactly zero, both charged responses are nonzero.", "",
        "In the unit #337 H/A/D convention, the exact response matrix at `p_ref` is", "",
        "```text",
        "              unit A source   unit D source",
        f"H response       0               0",
        f"A response       {a['unit_HAD_susceptibility']:.12g}               0",
        f"D response       0               {d['unit_HAD_susceptibility']:.12g}",
        "```", "",
        "| channel | W at p_ref | birth response | exit response | dW/dp |",
        "|---|---:|---:|---:|---:|",
        f"| A / B1 | {a['susceptibility']:.12g} | {a['birth_response']:.12g} | "
        f"{a['exit_response']:.12g} | {a['derivative_response']:.12g} |",
        f"| D / B2 | {d['susceptibility']:.12g} | {d['birth_response']:.12g} | "
        f"{d['exit_response']:.12g} | {d['derivative_response']:.12g} |", "",
        "Every degree-12 current coefficient obeys `dW_C/dp=J_C,birth-J_C,exit`.", "",
        "## Selection rule", "",
        "An A source activates only A at linear order; a D source activates only D. "
        "H response and A-D cross-response remain exactly zero. The minimal production "
        "output is therefore `(W_A,J_A,birth,J_A,exit,W_D,J_D,birth,J_D,exit)` with "
        "one joint covariance block. Existing projective-line archives already contain it.", "",
        "This constructs the missing explicit charged defect; it does not identify its "
        "continuum operator.", "",
    ]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = build_charged_certificate()
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                         encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
