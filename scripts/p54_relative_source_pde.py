#!/usr/bin/env python3
"""Exact source PDE and topology-sector inversion for the relative Q fugacity."""

from __future__ import annotations

import argparse
import json
import math
from functools import lru_cache
from fractions import Fraction
from pathlib import Path

from euler_motif_controls import configuration_identity
from integer_period_torus import gaussian_integer_torus


def fraction_record(value: Fraction) -> dict[str, int | str]:
    return {"numerator": value.numerator, "denominator": value.denominator, "text": str(value)}


def tiny_joint_distribution() -> list[tuple[int, int, Fraction]]:
    """N=5, p=1/2 distribution of (matching charge q, thermal score t)."""
    geometry = gaussian_integer_torus(2, 1)
    rows: dict[tuple[int, int], int] = {}
    for mask in range(1 << geometry.n):
        active = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
        record = configuration_identity(geometry, active, mask)
        occupied = record.motifs["V"]
        # This is the integer multiple of the Bernoulli score at p=1/2.
        thermal_score = 2 * occupied - geometry.n
        key = (record.q, thermal_score)
        rows[key] = rows.get(key, 0) + 1
    total = Fraction(1, 1 << geometry.n)
    return [(q, score, count * total) for (q, score), count in sorted(rows.items())]


def joint_moment(distribution: list[tuple[int, int, Fraction]], q_power: int, u_power: int) -> Fraction:
    return sum(weight * (q**q_power) * (score**u_power) for q, score, weight in distribution)


def cumulant_engine(distribution: list[tuple[int, int, Fraction]]):
    """Return exact joint cumulants kappa(q repeated a, score repeated b)."""

    @lru_cache(maxsize=None)
    def kappa(q_count: int, u_count: int) -> Fraction:
        if q_count + u_count == 0:
            return Fraction(0)
        moment = joint_moment(distribution, q_count, u_count)
        # Use a q as the distinguished variable when present.  The moment-cumulant
        # recursion sums over the block containing that distinguished variable.
        if q_count:
            subtotal = Fraction(0)
            for q_extra in range(q_count):
                for u_in_block in range(u_count + 1):
                    if q_extra == q_count - 1 and u_in_block == u_count:
                        continue
                    multiplicity = math.comb(q_count - 1, q_extra) * math.comb(u_count, u_in_block)
                    subtotal += (
                        multiplicity
                        * kappa(1 + q_extra, u_in_block)
                        * joint_moment(distribution, q_count - 1 - q_extra, u_count - u_in_block)
                    )
            return moment - subtotal

        # Pure score cumulants use a score as distinguished variable.
        subtotal = Fraction(0)
        for u_extra in range(u_count):
            if u_extra == u_count - 1:
                continue
            subtotal += (
                math.comb(u_count - 1, u_extra)
                * kappa(0, 1 + u_extra)
                * joint_moment(distribution, 0, u_count - 1 - u_extra)
            )
        return moment - subtotal

    return kappa


def build_oracle() -> dict:
    distribution = tiny_joint_distribution()
    kappa = cumulant_engine(distribution)

    sectors = {
        q: sum(weight for row_q, _, weight in distribution if row_q == q)
        for q in (-1, 0, 1)
    }
    mu = kappa(1, 0)
    variance = kappa(2, 0)
    raw_second = variance + mu**2
    reconstructed = {
        1: (raw_second + mu) / 2,
        -1: (raw_second - mu) / 2,
        0: 1 - raw_second,
    }

    raw_mixed = []
    for u_power in range(5):
        first = joint_moment(distribution, 1, u_power)
        third = joint_moment(distribution, 3, u_power)
        raw_mixed.append({
            "u_power": u_power,
            "d_s_d_u": fraction_record(first),
            "d_s3_d_u": fraction_record(third),
            "equal": first == third,
        })

    pure_closure_rhs = mu - 3 * mu * variance - mu**3
    mixed_first_u_rhs = (
        (1 - 3 * variance - 3 * mu**2) * kappa(1, 1)
        - 3 * mu * kappa(2, 1)
    )

    return {
        "schema": "matching-one.p54-relative-source-pde.v1",
        "issues": [54, 114],
        "support_identity": {
            "charge_support": [-1, 0, 1],
            "minimal_polynomial": "q(q-1)(q+1)=0",
            "configurationwise_identity": "q^3=q",
        },
        "exact_source_equations": {
            "relative_fugacity": "G(s,u)=E_u[exp(s q)]=P_0(u)+P_+(u) exp(s)+P_-(u) exp(-s), with G(0,u)=1",
            "linear_PDE": "partial_s^3 G=partial_s G",
            "mixed_linear_PDE": "partial_u^alpha partial_s^3 G=partial_u^alpha partial_s G for every commuting thermal/score multi-derivative alpha",
            "log_cumulant_closure": "F_sss=F_s-3 F_s F_ss-(F_s)^3, F=log G",
            "one_u_derivative": "F_sssu=(1-3 F_ss-3 F_s^2)F_su-3 F_s F_ssu",
            "derivative_recurrence": "all positive odd raw s derivatives equal G_s; all positive even raw s derivatives equal G_ss",
        },
        "sector_inversion": {
            "inputs": "mu=F_s and v=F_ss at s=0",
            "raw_second": "E[q^2]=v+mu^2",
            "P_plus": "(v+mu^2+mu)/2",
            "P_minus": "(v+mu^2-mu)/2",
            "P_zero": "1-v-mu^2",
            "meaning": "the first two connected Q-source derivatives reconstruct all three finite matching-topology sector weights",
        },
        "tiny_exact_oracle": {
            "geometry": "gaussian(2,1)",
            "N": 5,
            "p": "1/2",
            "thermal_score": "the unnormalized numerator source exp(u t) uses t=2*occupied-N, a fixed multiple of the Bernoulli score; division by its s-independent partition function gives normalized E_u and preserves the s-PDE",
            "joint_distribution": [
                {"q": q, "thermal_score": score, "probability": fraction_record(weight)}
                for q, score, weight in distribution
            ],
            "sector_weights": {str(q): fraction_record(value) for q, value in sectors.items()},
            "mu": fraction_record(mu),
            "variance": fraction_record(variance),
            "third_q_cumulant": fraction_record(kappa(3, 0)),
            "third_q_cumulant_from_closure": fraction_record(pure_closure_rhs),
            "sector_weights_from_mu_variance": {
                str(q): fraction_record(value) for q, value in reconstructed.items()
            },
            "mixed_raw_PDE_checks": raw_mixed,
            "first_thermal_log_closure": {
                "F_sssu_direct": fraction_record(kappa(3, 1)),
                "closure_rhs": fraction_record(mixed_first_u_rhs),
                "equal": kappa(3, 1) == mixed_first_u_rhs,
            },
        },
        "mixed_ledger_consequence": {
            "eliminated_independent_rows": "every row with at least three Q-source insertions is algebraically reducible to rows with at most two Q-source insertions",
            "retained_rows": "thermal/score derivative order is unrestricted and remains dynamical",
            "issue_54_use": "compute or estimate only the 0-, 1-, and 2-Q-source columns of the mixed Q/thermal ledger; generate higher-Q columns from the exact closure",
        },
        "claim_boundary": {
            "proved": [
                "the linear source PDE and every commuting mixed derivative",
                "the logarithmic cumulant closure",
                "the inversion from (F_s,F_ss) to the three q-sector weights",
                "the N5 finite oracle including one thermal-score insertion",
            ],
            "not_proved": [
                "the thermal third-derivative kappa3 sought by issue 54",
                "a relation reducing F_suuu to F_su",
                "a continuum value for any thermal cumulant",
            ],
            "critical_distinction": "F_sss closes three Q-source insertions; issue 54 requires three thermal insertions (schematically F_suuu), which this identity does not close",
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
