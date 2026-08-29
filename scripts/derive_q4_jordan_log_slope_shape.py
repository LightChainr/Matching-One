#!/usr/bin/env python3
"""Freeze the module-specific logarithmic-slope shape of the Q4 Jordan pair.

The exact inputs are the rank-2 dilation law, the thermal Q4 Ward coefficient,
and the area-normalized degree-2 E4 Hecke ratio.  Direct high-precision q-series
evaluation is retained as an independent numerical check of the CM targets.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp

from derive_rectangular_thermal_q4_hecke import exact_ratios
from derive_thermal_level4_torus_e4 import coefficients


def divisor_power_sum(n: int, power: int) -> int:
    total = 0
    divisor = 1
    while divisor * divisor <= n:
        if n % divisor == 0:
            partner = n // divisor
            total += divisor**power
            if partner != divisor:
                total += partner**power
        divisor += 1
    return total


def e4_qseries(tau: mp.mpc, *, dps: int = 90) -> tuple[mp.mpc, int, mp.mpf]:
    """Evaluate E4(tau)=1+240 sum sigma_3(n) q^n adaptively."""

    with mp.workdps(dps):
        q = mp.exp(2 * mp.pi * mp.j * tau)
        tolerance = mp.power(10, -(dps + 8))
        value = mp.mpc(1)
        small = 0
        for n in range(1, 100_001):
            term = 240 * divisor_power_sum(n, 3) * q**n
            value += term
            if abs(term) < tolerance:
                small += 1
                if small >= 8:
                    return +value, n, +abs(term)
            else:
                small = 0
    raise RuntimeError("E4 q-series did not converge")


def e4hat(tau: mp.mpc, *, dps: int = 90) -> tuple[mp.mpc, int, mp.mpf]:
    value, terms, last_term = e4_qseries(tau, dps=dps)
    return mp.im(tau) ** 2 * value, terms, last_term


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def complex_text(value: mp.mpc, digits: int = 70) -> dict[str, str]:
    return {
        "real": mp.nstr(mp.re(value), digits),
        "imag": mp.nstr(mp.im(value), digits),
    }


def render(dps: int = 90) -> dict:
    ward = coefficients()["Q4_over_primary_g2"]
    logn_over_primary_g2 = -ward / 2
    rectangular_ratio = exact_ratios()["E4hat_2i_over_E4hat_i"]

    with mp.workdps(dps):
        rho = (1 + mp.sqrt(3) * mp.j) / 2
        points = {
            "rho_exp_i_pi_over_3": rho,
            "i": mp.j,
            "2i": 2 * mp.j,
        }
        values = {}
        for name, tau in points.items():
            value, terms, last_term = e4hat(tau, dps=dps)
            values[name] = {
                "tau": complex_text(tau),
                "E4hat": complex_text(value),
                "terms_used": terms,
                "last_term_abs": mp.nstr(last_term, 12),
            }

        def as_mpc(row: dict) -> mp.mpc:
            return mp.mpc(row["E4hat"]["real"], row["E4hat"]["imag"])

        rho_value = as_mpc(values["rho_exp_i_pi_over_3"])
        square_value = as_mpc(values["i"])
        rectangle_value = as_mpc(values["2i"])
        ratio = rectangle_value / square_value

        return {
            "schema": "matching-one.q4-jordan-log-slope-shape.v1",
            "frozen_at": "2026-08-29",
            "claim_level": "C0_exact_conditional_bridge",
            "related": {"issue": 216, "stacked_on_pr": 217},
            "module_specific_derivation": {
                "rank2_dilation": "q_tilde -> s^(-x) [q_tilde-log(s) q]",
                "dimension": "x=21/4",
                "logN_slope": "B_logN(tau)=-(lambda_top/2)*A_q(tau)",
                "bottom_Ward_ratio": "A_q(tau)/A_epsilon(tau)=(493/96)*g2(tau)",
                "B_logN_over_lambda_top_A_epsilon_g2": fraction_text(
                    logn_over_primary_g2
                ),
                "root_normalized_shape": "B_root(tau)=C_J*Re[g2(tau)] (or C_J*g2 in a chiral projection)",
            },
            "frozen_predictions": {
                "hexagonal_zero": "B_root(exp(i*pi/3))=0",
                "rectangular_ratio": "B_root(2i)/B_root(i)=11/4",
                "rectangular_covariant_residual": "B_root(2i)-(11/4)*B_root(i)=0",
                "cross_channel_ratio": "[B_Sprime(tau1)/A_bottom(tau1)]/[B_Sprime(tau2)/A_bottom(tau2)]=1",
                "covariance_rule": "recompute slopes, bottom amplitudes, root normalization, ratios, and residuals inside each common replicate",
            },
            "exact_inputs": {
                "Q4_Ward_coefficient": fraction_text(ward),
                "logN_coefficient_after_Ward": fraction_text(
                    logn_over_primary_g2
                ),
                "area_normalized_E4hat_2i_over_i": fraction_text(
                    rectangular_ratio
                ),
            },
            "numerical_qseries": {
                "precision_decimal_digits": dps,
                "normalization": "E4hat(tau)=Im(tau)^2*E4(tau)",
                "values": values,
                "checks": {
                    "hexagonal_zero_abs": mp.nstr(abs(rho_value), 12),
                    "rectangular_ratio": complex_text(ratio),
                    "ratio_error_abs": mp.nstr(
                        abs(ratio - mp.mpf(11) / 4), 12
                    ),
                },
            },
            "interpretation_boundary": [
                "the_norm4_scale_cocycle_alone_cannot_identify_the_Q4_Jordan_module",
                "the_prediction_is_module_specific_only_when_log_scale_and_bottom_Q4_shape_are_jointly_satisfied",
                "a_generic_q2_coefficient_can_accidentally_share_an_E4_factor",
                "the_full_top_intercept_shape_A_tilde_tau_has_not_been_derived",
                "do_not_guess_weight_or_c_derivatives_E2_or_log_eta_terms_for_the_top_intercept",
                "lambda_top_C_J_and_absolute_lattice_amplitudes_remain_nonuniversal",
            ],
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    text = json.dumps(render(args.dps), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
