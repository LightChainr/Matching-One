#!/usr/bin/env python3
"""High-precision degree-2 Hecke fingerprints at the hexagonal CM point.

This is a modular-form calculation only.  It evaluates E4 and E6 from their
q-series and checks the same values against modular transformations and the
degree-2 Hecke relation.  No lattice observable is identified with either
modular form here.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp


def divisor_power_sum(n: int, power: int) -> int:
    """Return sigma_power(n) using paired divisors."""

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


def eisenstein_qseries(
    weight: int,
    tau: mp.mpc,
    *,
    dps: int = 90,
) -> tuple[mp.mpc, int, mp.mpf]:
    """Evaluate E4 or E6 with adaptive direct q-series truncation."""

    if weight not in (4, 6):
        raise ValueError("weight must be 4 or 6")
    with mp.workdps(dps):
        q = mp.exp(2 * mp.pi * mp.j * tau)
        coefficient = 240 if weight == 4 else -504
        tolerance = mp.power(10, -(dps + 8))
        value = mp.mpc(1)
        consecutive_small = 0
        last_term = mp.inf
        for n in range(1, 100_001):
            term = coefficient * divisor_power_sum(n, weight - 1) * q**n
            value += term
            last_term = abs(term)
            if last_term < tolerance:
                consecutive_small += 1
                if consecutive_small >= 8:
                    return +value, n, +last_term
            else:
                consecutive_small = 0
    raise RuntimeError("q-series did not converge within 100000 terms")


def normalized_eisenstein(weight: int, tau: mp.mpc, *, dps: int) -> tuple:
    value, terms, last_term = eisenstein_qseries(weight, tau, dps=dps)
    return mp.im(tau) ** (mp.mpf(weight) / 2) * value, terms, last_term


def _complex_payload(value: mp.mpc, digits: int = 70) -> dict[str, str]:
    return {
        "real": mp.nstr(mp.re(value), digits),
        "imag": mp.nstr(mp.im(value), digits),
    }


def fingerprint(dps: int = 90) -> dict:
    """Return the frozen numerical and exact hexagonal Hecke fingerprint."""

    with mp.workdps(dps):
        omega = (1 + mp.sqrt(3) * mp.j) / 2
        zeta = mp.exp(2 * mp.pi * mp.j / 3)
        points = [2 * omega, omega / 2, (omega + 1) / 2]
        names = ["2omega", "omega_over_2", "omega_plus_1_over_2"]

        numerical: dict[str, dict] = {}
        checks: dict[str, dict] = {}
        for weight in (4, 6):
            parent, parent_terms, parent_last = normalized_eisenstein(
                weight, omega, dps=dps
            )
            raw_parent, _, _ = eisenstein_qseries(weight, omega, dps=dps)
            children = []
            raw_children = []
            metadata = []
            for name, point in zip(names, points):
                raw, terms, last_term = eisenstein_qseries(weight, point, dps=dps)
                normalized = mp.im(point) ** (mp.mpf(weight) / 2) * raw
                raw_children.append(raw)
                children.append(normalized)
                metadata.append(
                    {
                        "name": name,
                        "tau": _complex_payload(point),
                        "raw_Ek": _complex_payload(raw),
                        "Ehat_k": _complex_payload(normalized),
                        "terms_used": terms,
                        "last_term_abs": mp.nstr(last_term, 12),
                    }
                )

            # S maps omega/2 to 2omega modulo T^-2.  Real q-coefficients
            # make the third child the reflection conjugate of the second.
            modular_second = (2**weight) * omega ** (-weight) * raw_children[0]
            modular_third = mp.conj(raw_children[1])
            hecke_lhs = (
                2 ** (weight - 1) * raw_children[0]
                + (raw_children[1] + raw_children[2]) / 2
            )
            hecke_rhs = (1 + 2 ** (weight - 1)) * raw_parent

            numerical[f"weight_{weight}"] = {
                "parent": {
                    "tau": _complex_payload(omega),
                    "raw_Ek": _complex_payload(raw_parent),
                    "Ehat_k": _complex_payload(parent),
                    "terms_used": parent_terms,
                    "last_term_abs": mp.nstr(parent_last, 12),
                },
                "children": metadata,
            }
            checks[f"weight_{weight}"] = {
                "S_then_T_error_abs": mp.nstr(
                    abs(raw_children[1] - modular_second), 12
                ),
                "reflection_error_abs": mp.nstr(
                    abs(raw_children[2] - modular_third), 12
                ),
                "Hecke_residual_abs": mp.nstr(abs(hecke_lhs - hecke_rhs), 12),
            }

        e4_children = [
            mp.mpc(row["Ehat_k"]["real"], row["Ehat_k"]["imag"])
            for row in numerical["weight_4"]["children"]
        ]
        e6_parent = mp.mpc(
            numerical["weight_6"]["parent"]["Ehat_k"]["real"],
            numerical["weight_6"]["parent"]["Ehat_k"]["imag"],
        )
        e6_children = [
            mp.mpc(row["Ehat_k"]["real"], row["Ehat_k"]["imag"])
            for row in numerical["weight_6"]["children"]
        ]
        checks["fingerprint_errors"] = {
            "E4_parent_zero_abs": mp.nstr(
                abs(
                    mp.mpc(
                        numerical["weight_4"]["parent"]["Ehat_k"]["real"],
                        numerical["weight_4"]["parent"]["Ehat_k"]["imag"],
                    )
                ),
                12,
            ),
            "E4_child_phase_max_abs": mp.nstr(
                max(
                    abs(e4_children[1] / e4_children[0] - zeta),
                    abs(e4_children[2] / e4_children[0] - zeta**2),
                ),
                12,
            ),
            "E6_child_ratio_max_abs": mp.nstr(
                max(abs(child / e6_parent - mp.mpf(11) / 4) for child in e6_children),
                12,
            ),
        }

        return {
            "schema": "matching-one.hexagonal-degree2-hecke.v1",
            "frozen_at": "2026-08-29",
            "claim_level": "C0",
            "related_issues": [161, 164],
            "precision_decimal_digits": dps,
            "convention": {
                "omega": "exp(i*pi/3)=(1+i*sqrt(3))/2",
                "zeta": "exp(2*pi*i/3)=omega^2",
                "period_basis": "(1,tau)",
                "child_order": names,
                "normalization": "Ehat_k(tau)=Im(tau)^(k/2)*E_k(tau)",
            },
            "exact_fingerprint": {
                "E4_parent": "Ehat_4(omega)=0",
                "E4_children": "A*(1,zeta,zeta^2), A=Ehat_4(2omega)",
                "E4_reflection_even_real_pattern": "A*(1,-1/2,-1/2)",
                "E4_reflection_odd_imag_pattern": "A*(0,sqrt(3)/2,-sqrt(3)/2)",
                "E6_children_over_parent": ["11/4", "11/4", "11/4"],
                "E6_raw_children_over_parent": ["11/32", "22", "22"],
            },
            "numerical_qseries": numerical,
            "independent_checks": checks,
            "interpretation_boundary": {
                "shape_ratio": "11/4 is an area-normalized degree-2 shape ratio",
                "similarity_scaling": "under scalar similarity m, weight k scales as m^(-k)",
                "spin4": "child differences/projectors remove a common scalar child contribution",
                "spin6": "the constant child vector needs Issue #161 scalar normalization",
                "claim": "exact modular-form identity, not a percolation identification",
            },
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    payload = fingerprint(args.dps)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.json is None:
        print(rendered, end="")
    else:
        args.json.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
