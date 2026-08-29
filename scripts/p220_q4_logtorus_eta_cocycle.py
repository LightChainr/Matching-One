#!/usr/bin/env python3
"""Exact eta-cocycle fingerprint for the Q4 energy-block derivative."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp


def divisor_power_sum(n: int, power: int) -> int:
    total = 0
    divisor = 1
    while divisor * divisor <= n:
        if n % divisor == 0:
            other = n // divisor
            total += divisor**power
            if other != divisor:
                total += other**power
        divisor += 1
    return total


def eta(tau: mp.mpc, *, dps: int = 90) -> mp.mpc:
    """Dedekind eta from its q product, with q=exp(2 pi i tau)."""
    with mp.workdps(dps):
        q = mp.exp(2 * mp.pi * mp.j * tau)
        value = mp.exp(mp.pi * mp.j * tau / 12)
        tolerance = mp.power(10, -(dps + 8))
        small = 0
        for n in range(1, 100_001):
            value *= 1 - q**n
            if abs(q**n) < tolerance:
                small += 1
                if small >= 8:
                    return +value
            else:
                small = 0
    raise RuntimeError("eta product did not converge")


def e4(tau: mp.mpc, *, dps: int = 90) -> mp.mpc:
    with mp.workdps(dps):
        q = mp.exp(2 * mp.pi * mp.j * tau)
        value = mp.mpc(1)
        tolerance = mp.power(10, -(dps + 8))
        small = 0
        for n in range(1, 100_001):
            term = 240 * divisor_power_sum(n, 3) * q**n
            value += term
            if abs(term) < tolerance:
                small += 1
                if small >= 8:
                    return +value
            else:
                small = 0
    raise RuntimeError("E4 series did not converge")


def fixed_q4_coefficient(h: Fraction) -> Fraction:
    """<Q4 phi>/<phi>/g2 along the level-2 degenerate family.

    This keeps the repository's c=0 Q4 coefficients fixed.  Its h derivative
    is normalization-path dependent and is not used in the shape cocycle.
    """
    return h * (Fraction(111) + Fraction(120, 2 * h + 1)) / 20


def fixed_q4_coefficient_derivative(h: Fraction) -> Fraction:
    return Fraction(111, 20) + Fraction(6, (2 * h + 1) ** 2)


def mptext(value: mp.mpf | mp.mpc, digits: int = 70) -> str | dict[str, str]:
    if isinstance(value, mp.mpc):
        return {"real": mp.nstr(mp.re(value), digits), "imag": mp.nstr(mp.im(value), digits)}
    return mp.nstr(value, digits)


def analyze(dps: int = 90) -> dict:
    with mp.workdps(dps):
        tau_i = mp.j
        tau_2i = 2 * mp.j
        tau_shear = (1 + mp.j) / 2
        tau_rho = (1 + mp.sqrt(3) * mp.j) / 2
        values = {name: eta(tau, dps=dps) for name, tau in {
            "i": tau_i,
            "2i": tau_2i,
            "(1+i)/2": tau_shear,
            "rho": tau_rho,
        }.items()}
        log_shape_2i = 2 * mp.log(abs(values["2i"] / values["i"]))
        log_shape_shear = 2 * mp.log(abs(values["(1+i)/2"] / values["i"]))
        rational_target_2i = -mp.mpf(3) * mp.log(2) / 4
        rational_target_shear = mp.log(2) / 2
        ratio = log_shape_2i / log_shape_shear
        rho_e4 = e4(tau_rho, dps=dps)

        h = Fraction(5, 8)
        coefficient = fixed_q4_coefficient(h)
        coefficient_derivative = fixed_q4_coefficient_derivative(h)
        return {
            "schema": "matching-one/q4-logtorus-eta-cocycle/v1",
            "issue": 220,
            "status": "exact_energy_block_derivative_and_frozen_gauge_free_shape_ratio",
            "generic_identity": {
                "bottom": "H_Q4(u,tau)=C_Q4(u) g2(tau) |eta(tau)|^(4h(u))",
                "tangent_over_bottom": "R_u(tau)=(partial_u H_Q4)/H_Q4=partial_u log C_Q4+4 h'(u) log|eta(tau)|",
                "dimension_velocity": "x'(u)=2h'(u)",
                "gauge_free_difference": "[R_u(tau1)-R_u(tau0)]/x'(u)=2 log|eta(tau1)/eta(tau0)|",
                "top_shift_invariance": "H_top -> H_top+alpha H_bottom adds a constant to R and cancels in every modulus difference",
            },
            "fixed_Q4_path_audit": {
                "C_Q4_at_h_5_8": str(coefficient),
                "dC_Q4_dh_at_h_5_8": str(coefficient_derivative),
                "dlogC_Q4_dh_at_h_5_8": str(coefficient_derivative / coefficient),
                "boundary": "these constants depend on the generic descendant-normalization path; the frozen modulus difference does not",
            },
            "exact_CM_eta_identities": {
                "eta(2i)/eta(i)": "2^(-3/8)",
                "abs_eta((1+i)/2)/eta(i)": "2^(1/4)",
                "proof_inputs": [
                    "theta2(tau)=2 eta(2tau)^2/eta(tau)",
                    "lambda(i)=theta2(i)^4/theta3(i)^4=1/2",
                    "eta(i)^3=theta2(i)theta3(i)theta4(i)/2 with theta2(i)=theta4(i)",
                    "eta(gamma tau)=epsilon(gamma)(c tau+d)^(1/2)eta(tau), gamma=[[1,0],[1,1]]",
                ],
            },
            "frozen_predictions": {
                "Xi_2i_i": "-3 log(2)/4",
                "Xi_shear_i": "log(2)/2",
                "amplitude_and_velocity_free_ratio": "Xi(2i,i)/Xi((1+i)/2,i)=-3/2",
                "hexagonal_joint_zero": "H_bottom(rho)=partial_u H_bottom(rho)=0 because E4(rho)=0",
                "scoring_definition": "Xi(tau1,tau0)=[H_top/H_bottom(tau1)-H_top/H_bottom(tau0)]/x_collision_velocity",
            },
            "numerical_oracle": {
                "dps": dps,
                "eta_values": {name: mptext(value) for name, value in values.items()},
                "Xi_2i_i": mptext(log_shape_2i),
                "Xi_shear_i": mptext(log_shape_shear),
                "ratio": mptext(ratio),
                "errors": {
                    "Xi_2i_i": mptext(abs(log_shape_2i - rational_target_2i)),
                    "Xi_shear_i": mptext(abs(log_shape_shear - rational_target_shear)),
                    "ratio_minus_minus3over2": mptext(abs(ratio + mp.mpf(3) / 2)),
                    "E4_rho_abs": mptext(abs(rho_e4)),
                },
            },
            "claim_boundary": {
                "exact": "the derivative of the generic Potts energy eta block after a homogeneous level-4 torus Ward operation",
                "conditional_bridge": "the lattice top and bottom observables must realize the same energy/[2] collision normalization across moduli",
                "not_claimed": "that an arbitrary Q4 Jordan top contains no additional torus blocks, or that P4[S'] already equals this field",
            },
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.dumps(analyze(args.dps), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")


if __name__ == "__main__":
    main()
