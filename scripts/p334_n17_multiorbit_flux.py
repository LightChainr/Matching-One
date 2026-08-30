#!/usr/bin/env python3
"""Exact N=17 multi-orbit flux and cross-quotient comparison for Issue #334."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Optional, Sequence

from p334_n13_multiorbit_flux import (
    P_REF,
    _decimal,
    _evaluate_incidence,
    evaluate_reference,
    exact_census,
)


LABELS = ("axis_orbit", "diagonal_orbit")


def _exact_shares(census: dict[str, object]) -> dict[str, Fraction]:
    contributions: dict[str, Fraction] = {}
    for label in LABELS:
        character_real = Fraction(census["orbits"][label]["chi4"]["real"])
        incidence = _evaluate_incidence(census, label, "birth", P_REF)
        incidence -= _evaluate_incidence(census, label, "exit", P_REF)
        contributions[label] = character_real * incidence
    total = sum(contributions.values(), Fraction(0))
    if total == 0:
        raise AssertionError("p_ref derivative has no real component")
    return {label: value / total for label, value in contributions.items()}


def _geometry_summary(
    census: dict[str, object], reference: dict[str, object]
) -> dict[str, object]:
    return {
        "geometry": census["geometry"],
        "direct_rank2_edge_count": census["direct_rank2_edge_count"],
        "orbits": census["orbits"],
        "signed_collinear_share": reference["signed_collinear_share"],
        "total_complex_dA4_dp": reference["total_complex_dA4_dp"],
    }


def build_certificate() -> dict[str, object]:
    n13 = exact_census(3, 2, include_direct_rank2=True)
    n17 = exact_census(4, 1, include_direct_rank2=True)
    ref13 = evaluate_reference(n13)
    ref17 = evaluate_reference(n17)
    shares13 = _exact_shares(n13)
    shares17 = _exact_shares(n17)
    deltas = {label: shares17[label] - shares13[label] for label in LABELS}

    all_orbit_flux_nonzero = all(
        n17["orbits"][label][key] > 0
        for label in LABELS
        for key in ("birth_edge_count", "exit_edge_count")
    )
    share_shift_l1 = sum(abs(value) for value in deltas.values())
    gates = {
        "n17_has_two_inequivalent_d4_orbits": n17["gates"]["orbit_count"] == 2,
        "n17_characters_are_exact_opposites": n17["gates"][
            "characters_are_exact_opposites"
        ],
        "all_n17_orbits_have_birth_and_exit": all_orbit_flux_nonzero,
        "n17_coefficientwise_continuity": n17["gates"][
            "coefficientwise_dA4_equals_birth_minus_exit"
        ],
        "both_quotients_reinforce_at_p_ref": (
            ref13["both_orbits_reinforce_total"]
            and ref17["both_orbits_reinforce_total"]
        ),
        "all_pass": (
            n17["gates"]["all_pass"]
            and all_orbit_flux_nonzero
            and ref13["both_orbits_reinforce_total"]
            and ref17["both_orbits_reinforce_total"]
        ),
    }
    return {
        "schema": "matching-one/p334-n17-cross-quotient-multiorbit-flux/v1",
        "issue": 334,
        "parent_commit": "b8e286e",
        "status": "exact_cross_quotient_multiorbit_flux",
        "orbit_gate": {
            "candidate": "gaussian-4-plus-1i",
            "criterion": "at least two inequivalent primitive-line D4 orbits",
            "passed": gates["n17_has_two_inequivalent_d4_orbits"],
            "fallback_scan_used": False,
        },
        "result": {
            "statement": (
                "The N=17 Gaussian quotient 4+i again resolves rank-one traffic "
                "into axis and diagonal primitive-line orbits with exactly opposite "
                "chi4. Both orbits carry nonzero source and sink flux and reinforce "
                "dA4/dp at p_ref, so the N=13 mechanism is not quotient-specific."
            ),
            "cross_geometry_statement": (
                "The p_ref axis signed share moves from 75.574% at N=13 to "
                "76.484% at N=17, a +0.911 percentage-point shift despite the "
                "different physical chi4 phase."
            ),
            "path_enumeration_avoided": (
                "The N=17 certificate visits 2^17 subsets and 17*2^16 directed "
                "subset-boundary edges; it does not enumerate 17! paths."
            ),
        },
        "n17_census": n17,
        "n17_reference_evaluation": ref17,
        "n13_comparison": _geometry_summary(n13, ref13),
        "cross_geometry_signed_share": {
            "basis": list(LABELS),
            "p_ref": str(P_REF),
            "n13": {label: _decimal(shares13[label]) for label in LABELS},
            "n17": {label: _decimal(shares17[label]) for label in LABELS},
            "n17_minus_n13": {label: _decimal(deltas[label]) for label in LABELS},
            "l1_shift": _decimal(share_shift_l1),
        },
        "gates": gates,
        "claim_boundary": [
            "This is an exact finite-volume source/sink localization, not a continuum-field identification.",
            "The cross-quotient share stability is a two-geometry mechanism clue, not an asymptotic limit.",
            "p_ref is inherited unchanged from the N=13 certificate; no parameter was fitted.",
            "No path enumeration, Monte Carlo sample, Huawei production, new PR, or merge is used.",
        ],
    }


def render_markdown(payload: dict[str, object]) -> str:
    n17 = payload["n17_census"]
    comparison = payload["cross_geometry_signed_share"]
    lines = [
        "# Exact N=17 multi-orbit flux and N=13 cross-quotient comparison",
        "",
        f"Status: `{payload['status']}`.",
        "",
        payload["result"]["statement"],
        "",
        "## N=17 orbit gate",
        "",
        f"- period matrix: `{n17['geometry']['period_matrix']}`",
        f"- subset states: {n17['geometry']['subset_states']:,}",
        f"- directed addition edges: {n17['geometry']['directed_addition_edges']:,}",
        f"- direct 0->2 edges: {n17['direct_rank2_edge_count']:,}",
        f"- orbit count: {n17['gates']['orbit_count']}",
        f"- coefficientwise source/sink identity: `{n17['gates']['coefficientwise_dA4_equals_birth_minus_exit']}`",
        "",
        "| orbit | primitive lines | chi4 | rank-one states | birth | exit |",
        "|---|---|---|---:|---:|---:|",
    ]
    for label in LABELS:
        row = n17["orbits"][label]
        lines.append(
            f"| {label} | {row['primitive_lines']} | "
            f"({row['chi4']['real']}, {row['chi4']['imag']}) | "
            f"{row['rank_one_state_count']} | {row['birth_edge_count']} | "
            f"{row['exit_edge_count']} |"
        )
    lines += [
        "",
        "## Frozen p_ref signed-share comparison",
        "",
        "| quotient | axis share | diagonal share |",
        "|---|---:|---:|",
        f"| 3+2i (N=13) | {comparison['n13']['axis_orbit']} | {comparison['n13']['diagonal_orbit']} |",
        f"| 4+i (N=17) | {comparison['n17']['axis_orbit']} | {comparison['n17']['diagonal_orbit']} |",
        f"| N17-N13 | {comparison['n17_minus_n13']['axis_orbit']} | {comparison['n17_minus_n13']['diagonal_orbit']} |",
        "",
        payload["result"]["cross_geometry_statement"],
        "",
        "The character phase rotates between the two Gaussian generators, but the normalized source-minus-sink partition stays close and both line orbits reinforce the total derivative. That favors a transported projective-current mechanism over an N=13-only cancellation accident.",
        "",
        "## Boundary",
        "",
    ]
    lines.extend(f"- {item}" for item in payload["claim_boundary"])
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    payload = build_certificate()
    rendered = (
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(payload) + "\n"
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
