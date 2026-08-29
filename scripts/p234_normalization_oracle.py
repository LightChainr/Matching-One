#!/usr/bin/env python3
"""Exact bookkeeping oracle for the P234 natural-cutoff normalization."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def normalization_dictionary(
    C1: float,
    C2: float,
    CL: float,
    G2: float = 1.0,
    F3: float = 1.0,
) -> dict[str, float]:
    if min(C1, C2, CL, G2, F3) <= 0.0:
        raise ValueError("C1, C2, CL, G2, and F3 must be positive")
    spin_pair_amplitude = math.sqrt(C1)
    mixed = C2 * G2 / spin_pair_amplitude
    top_top_slope = -2.0 * CL * G2
    proxy = -top_top_slope / (2.0 * mixed)
    top_spin_spin_slope = -CL * F3
    invariant_CL = -2.0 * (top_spin_spin_slope / F3) ** 2 / (
        top_top_slope / G2
    )
    shear_gate = (top_top_slope / G2) / (
        2.0 * top_spin_spin_slope / F3
    )
    return {
        "spin_pair_amplitude_sqrt_C1": spin_pair_amplitude,
        "paper_kappa_C1_CL_over_C2": C1 * CL / C2,
        "LD_continuum": mixed,
        "DD_log_2delta_slope": top_top_slope,
        "kappa_proxy_sqrt_C1_CL_over_C2": proxy,
        "top_spin_spin_log_2delta_slope": top_spin_spin_slope,
        "CL_from_gauge_invariant_extra_observable": invariant_CL,
        "same_shear_gate": shear_gate,
    }


def gauge_exponents() -> dict[str, list[int]]:
    """Powers of independent (lambda_phi, mu_top) field rescalings."""
    return {
        "mixed_phi_top": [1, 1],
        "top_top_slope": [0, 2],
        "top_spin_spin_slope": [0, 1],
        "kappa_proxy": [-1, 1],
        "CL_invariant_minus_2_t3_squared_over_s2": [0, 0],
    }


def render(partial: float, standard_error: float) -> dict[str, object]:
    if standard_error <= 0.0:
        raise ValueError("standard error must be positive")
    target = 8.0 / 3.0
    return {
        "schema": "matching-one.p234-natural-cutoff-normalization-audit.v1",
        "issue": 234,
        "status": "theory_audit_partial_not_final",
        "paper_convention": {
            "spin_two_point_coefficient": "sqrt(C1)",
            "paper_kappa": "C1*CL/C2",
        },
        "scorer_dictionary_after_p_connection_normalization": {
            "LD_continuum": "(C2/sqrt(C1))*G2",
            "DD_log_2delta_slope": "-2*CL*G2",
            "kappa_proxy": "sqrt(C1)*CL/C2",
        },
        "partial_8_over_3_diagnostic": {
            "estimate": partial,
            "standard_error": standard_error,
            "target": target,
            "z_score_estimate_minus_target": (partial - target) / standard_error,
            "classification": "high_risk_amplitude_conjecture_not_exponent_theorem",
            "identity_if_true": "C2=(3/8)*sqrt(C1)*CL",
        },
        "gauge_exponents_lambda_phi_mu_top": gauge_exponents(),
        "minimal_extra_observable": {
            "definition": "t3=d<T_delta chi chi>/dlog(2delta), chi=psi/C1^(1/4)",
            "continuum": "t3=-CL*F3",
            "existing_top_top_slope": "s2=-2*CL*G2",
            "gauge_invariant_CL": "-2*(t3/F3)^2/(s2/G2)",
            "same_shear_gate_in_pconn_gauge": "(s2/G2)/(2*t3/F3)=1",
        },
        "numeric_self_check": normalization_dictionary(9.0, 5.0, 7.0, 11.0, 13.0),
        "scope": [
            "No frozen scorer or run artifact is modified.",
            "The supplied 2.653 +/- 0.448 value is partial and not a final score.",
            "Plane shape factors are explicit; individual torus amplitudes need a torus-to-plane geometry bridge.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partial", type=float, default=2.653)
    parser.add_argument("--standard-error", type=float, default=0.448)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis/p234_natural_cutoff_normalization_audit.json"),
    )
    args = parser.parse_args()
    payload = render(args.partial, args.standard_error)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
