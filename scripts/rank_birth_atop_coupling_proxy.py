#!/usr/bin/env python3
"""Exact tiny connected-coupling proxy between A_top and rank-birth H4 fields."""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Optional, Sequence

from homology_rank_birth_insertion import (
    _fraction_text,
    _mask_active,
    geometry_specs,
)
from marked_pivotal_h4_reference import landing_mark
from rank_birth_parity_channels import (
    _gate_record,
    _rank_states,
    spin4_character,
)


SOURCE_NAMES = (
    "S",
    "D",
    "line_cos4_S",
    "line_sin4_S",
    "line_cos4_D",
    "line_sin4_D",
    "landing_h4_S",
    "landing_h4_D",
)

EVEN_SOURCES = ("S", "line_cos4_S", "line_sin4_S", "landing_h4_S")
ODD_SOURCES = ("D", "line_cos4_D", "line_sin4_D", "landing_h4_D")

ROOT = Path(__file__).resolve().parents[1]


def _source_at_root(
    geometry,
    states,
    mask: int,
    vertex: int,
    *,
    matching: bool,
    local_radius: Optional[int],
) -> dict[str, Fraction]:
    environment = mask & ~(1 << vertex)
    record = _gate_record(states[environment], states[environment | (1 << vertex)])
    even = Fraction(record["even"])
    odd = Fraction(record["odd"])
    values = {name: Fraction(0) for name in SOURCE_NAMES}
    values["S"] = even
    values["D"] = odd
    if record["ell"] is not None:
        physical = geometry.periods.period_vector(record["ell"])
        cos4, sin4 = spin4_character(physical)
        values["line_cos4_S"] = even * cos4
        values["line_sin4_S"] = even * sin4
        values["line_cos4_D"] = odd * cos4
        values["line_sin4_D"] = odd * sin4
    if local_radius is not None:
        active = _mask_active(environment, geometry.n)
        h4 = landing_mark(
            geometry,
            active,
            local_radius,
            open_matching=matching,
        )["h4"]
        values["landing_h4_S"] = even * h4
        values["landing_h4_D"] = odd * h4
    return values


def _side_moments(
    geometry,
    *,
    roots: Sequence[int],
    root_multiplicity: int,
    local_radius: Optional[int],
    matching: bool,
) -> dict[str, Any]:
    states = _rank_states(geometry, matching=matching)
    configuration_count = 1 << geometry.n
    sum_q = Fraction(0)
    sum_q2 = Fraction(0)
    sums = {name: Fraction(0) for name in SOURCE_NAMES}
    sum_q_source = {name: Fraction(0) for name in SOURCE_NAMES}

    for mask in range(configuration_count):
        q = Fraction(states[mask][0] - 1)
        source = {name: Fraction(0) for name in SOURCE_NAMES}
        for vertex in roots:
            root_values = _source_at_root(
                geometry,
                states,
                mask,
                vertex,
                matching=matching,
                local_radius=local_radius,
            )
            for name, value in root_values.items():
                source[name] += root_multiplicity * value
        sum_q += q
        sum_q2 += q * q
        for name, value in source.items():
            sums[name] += value
            sum_q_source[name] += q * value

    denominator = Fraction(configuration_count)
    mean_q = sum_q / denominator
    variance_q = sum_q2 / denominator - mean_q * mean_q
    means = {name: value / denominator for name, value in sums.items()}
    covariances = {
        name: sum_q_source[name] / denominator - mean_q * means[name]
        for name in SOURCE_NAMES
    }
    birth_mass = means["S"]
    if birth_mass == 0:
        raise AssertionError("unmarked rank-birth mass vanished")
    normalized = {
        name: covariances[name] / birth_mass
        for name in ("line_cos4_D", "line_sin4_D", "landing_h4_D")
    }
    return {
        "matching_backend": matching,
        "configuration_count": configuration_count,
        "root_reduction": {
            "roots": list(roots),
            "multiplicity": root_multiplicity,
            "statement": (
                "translation reduction is exact for means and covariance with the translation-invariant q"
                if root_multiplicity > 1
                else "every root is summed explicitly"
            ),
        },
        "mean_A_top": _fraction_text(mean_q),
        "variance_A_top": _fraction_text(variance_q),
        "source_means": {
            name: _fraction_text(value) for name, value in means.items()
        },
        "connected_covariance_with_A_top": {
            name: _fraction_text(value) for name, value in covariances.items()
        },
        "normalized_D_H4_proxy_over_birth_mass": {
            name: _fraction_text(value) for name, value in normalized.items()
        },
    }


def _geometry_summary(spec: dict[str, Any]) -> dict[str, Any]:
    geometry = spec["geometry"]
    common = {
        "roots": spec["roots"],
        "root_multiplicity": spec["root_multiplicity"],
        "local_radius": spec["local_radius"],
    }
    primal = _side_moments(geometry, matching=False, **common)
    matching = _side_moments(geometry, matching=True, **common)
    primal_cov = {
        name: Fraction(value)
        for name, value in primal["connected_covariance_with_A_top"].items()
    }
    matching_cov = {
        name: Fraction(value)
        for name, value in matching["connected_covariance_with_A_top"].items()
    }
    complement_residuals = {}
    for name in EVEN_SOURCES:
        # A_top is odd and an S-source is even, so the covariance is odd.
        complement_residuals[name] = primal_cov[name] + matching_cov[name]
    for name in ODD_SOURCES:
        # A_top and a D-source are both odd, so their covariance is even.
        complement_residuals[name] = primal_cov[name] - matching_cov[name]
    if any(complement_residuals.values()):
        raise AssertionError(f"connected complement parity failed: {complement_residuals}")
    return {
        "id": spec["name"],
        "N": geometry.n,
        "period_matrix": [list(row) for row in geometry.periods.matrix],
        "primal": primal,
        "matching": matching,
        "complement_covariance_residuals": {
            name: _fraction_text(value)
            for name, value in complement_residuals.items()
        },
        "finite_lattice_nonzero": {
            "line_D_complex": (
                primal_cov["line_cos4_D"] != 0
                or primal_cov["line_sin4_D"] != 0
            ),
            "landing_D": primal_cov["landing_h4_D"] != 0,
        },
    }


def scaling_dictionary() -> dict[str, Any]:
    x = Fraction(21, 4)
    y_thermal = Fraction(3, 4)
    integrated_n = Fraction(2 - x, 2)
    local_n = -x / 2
    raw_pivotal_sum = y_thermal / 2 + integrated_n
    raw_pivotal_site = raw_pivotal_sum - 1
    return {
        "candidate_dimension_x": _fraction_text(x),
        "direct_CFT_density_normalization": {
            "one_local_insertion_N_power": _fraction_text(local_n),
            "sum_over_N_sites_N_power": _fraction_text(integrated_n),
            "derivation": "N^(-x/2) locally and N^(1-x/2) after summing N sites",
        },
        "rank_birth_measure_normalization": {
            "assumption": "critical unmarked birth mass B_N=M_prime scales as L^(1/nu)=N^(3/8), nu=4/3",
            "normalized_proxy_C_over_B_N_power": _fraction_text(integrated_n),
            "raw_connected_sum_C_N_power": _fraction_text(raw_pivotal_sum),
            "raw_per_site_contribution_N_power": _fraction_text(raw_pivotal_site),
            "values": "C/B ~ N^-13/8; C ~ N^-5/4; C/N ~ N^-9/4",
        },
        "jordan_variant": "multiply the declared powers by an affine A+B log N factor where the top partner is read out",
    }


def archived_birth_snapshot() -> dict[str, Any]:
    """Consume the bedc94b archive view and expose its precise stopping point."""

    path = ROOT / "results/essential-birth-histogram/latest.json"
    archive = json.loads(path.read_text(encoding="utf-8"))
    root = archive["archived_root_evaluation"]
    return {
        "source": str(path.relative_to(ROOT)),
        "source_N": archive["source_N"],
        "source_sample_count": archive["source_sample_count"],
        "balance_root_p": root["p"],
        "available_unmarked_birth_mass_B_equals_M_prime": root["M_prime"],
        "available_unmarked_D_equals_f12_minus_f01": str(
            -Decimal(root["rank_one_derivative"])
        ),
        "connected_D_H4_proxy": None,
        "reason": "ell/local H4 and same-sample A_top times J cross moments were not retained",
    }


def build_artifact() -> dict[str, Any]:
    geometries = [_geometry_summary(spec) for spec in geometry_specs()]
    return {
        "schema": "matching-one/rank-birth-atop-coupling-proxy/v1",
        "issues": [215, 275],
        "base_commit": "dabe28e mechanism assets rebased over main c500640",
        "status": "tiny_exact_connected_proxy_and_conditional_CFT_bridge",
        "finite_lattice_definition": {
            "A_top_sample": "q(omega)=rank(omega)-1",
            "odd_spin4_source": "J_D4=sum_v chi4(ell_v) (I_12-I_01)_v, with a separate landing-H4 version",
            "tilted_measure": "P_h(omega) proportional to P_p(omega) exp(h Re[e^(-4 i phi) J_D4(omega)])",
            "response": "partial_h E_h[A_top]|_0 = Cov_p(A_top, Re[e^(-4 i phi) J_D4])",
            "normalized_proxy": "gamma_D4=Cov(A_top,J_D4)/B_N with B_N=E[sum_v(I_01+I_12)]=M_prime",
        },
        "symmetry": {
            "allowed": "A_top and J_D4 are both complement-odd, so their connected correlator is complement-even",
            "null_sign_control": "A_top is odd while J_S4 is even, so the primal and matching connected responses have opposite sign",
            "rotation": "chi4 supplies spin four and changes sign under multiplication of the physical line by 1+i",
        },
        "tiny_exact": geometries,
        "scaling_if_Q4_epsilon": scaling_dictionary(),
        "archive_reconstructibility": {
            "bedc94b_snapshot": archived_birth_snapshot(),
            "bedc94b_available": [
                "K1/K2 marginals and joint histogram",
                "unmarked B_N=M_prime",
                "unmarked D_N=-partial_p P(rank=1)",
            ],
            "missing": [
                "ell/iota at the two births",
                "line or landing H4 source J_D4",
                "same-sample A_top times J_D4 cross moment",
            ],
            "conclusion": "the committed rank and marked-pivotal archives cannot reconstruct the connected coupling proxy",
        },
        "next_sufficient_statistics": {
            "estimator": "one uniform random root per Bernoulli configuration, multiplied by N, is unbiased for the site sum",
            "per_sample": [
                "global q=rank-1",
                "root gate pair I_01/I_12 and direct-0-to-2 flag",
                "ell, iota, chi4(ell) when non-null",
                "landing H4 mark at the same root",
                "J_S4 and J_D4 real/imaginary coordinates",
            ],
            "per_aligned_batch_raw_moments": [
                "sum q and sum q^2",
                "sum J for every S/D line and landing coordinate",
                "sum q*J for every coordinate",
                "birth mass sum S",
            ],
            "common_field": "reuse the same configuration counter, random root, and delete-one batch across paired orientations",
        },
        "claim_layers": {
            "exact": "finite exponential-tilt response equals the connected covariance; complement parity and tiny nonzero values",
            "conditional_CFT_bridge": "gamma_D4 scaling and modulus phase identify Q4 epsilon only if J_D4 flows to that field",
            "conjecture": "the line-resolved D source has nonzero overlap with the thermal Q4 epsilon/Jordan module",
        },
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    lines = [
        "# A_top--rank-birth H4 connected coupling proxy",
        "",
        "For `J_D4=sum chi4(ell)(I12-I01)`, an exponential source gives exactly",
        "",
        "```text",
        "d/dh E_h[A_top]|0 = Cov(A_top,J_D4).",
        "```",
        "",
        "Both factors are complement-odd, so this response is allowed/even. The S-H4",
        "control is complement-odd and reverses sign between primal and matching backends.",
        "",
        "| geometry | Re Cov(A,J_D4) | Im Cov(A,J_D4) | Cov(A,J_Dlanding) | Re gamma=C/B | Im gamma=C/B |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in artifact["tiny_exact"]:
        cov = row["primal"]["connected_covariance_with_A_top"]
        gamma = row["primal"]["normalized_D_H4_proxy_over_birth_mass"]
        lines.append(
            f"| {row['id']} | {cov['line_cos4_D']} | {cov['line_sin4_D']} | "
            f"{cov['landing_h4_D']} | {gamma['line_cos4_D']} | {gamma['line_sin4_D']} |"
        )
    lines.extend(
        [
            "",
            "The tiny proxy is nonzero, proving that no finite-lattice symmetry forbids the",
            "coupling. It does not identify the continuum field. If it flows to x=21/4, the",
            "birth-mass-normalized proxy scales as N^-13/8 (with an optional affine log N Jordan",
            "factor); raw summed and per-site rank-birth correlators scale as N^-5/4 and N^-9/4",
            "assuming M'~N^3/8.",
            "",
        ]
    )
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    artifact = build_artifact()
    rendered = (
        json.dumps(artifact, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(artifact)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
