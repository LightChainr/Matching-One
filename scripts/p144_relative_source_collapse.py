#!/usr/bin/env python3
"""Exact terminal relative-source collapse of the P144 incidence spine."""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from p144_typed_incidence_spine import build_oracle as build_spine_oracle  # noqa: E402


def fraction_record(value: Fraction) -> dict[str, int | str]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "text": str(value),
    }


def axis_L3_relative_collapse(p: Fraction = Fraction(2, 5)) -> dict[str, object]:
    if not 0 < p < 1:
        raise ValueError("exact check requires 0<p<1")
    spine = build_spine_oracle()["honest_torus_exact_oracle"]
    N = int(spine["N"])
    sectors = {-1: Fraction(0), 0: Fraction(0), 1: Fraction(0)}
    terminal_support: set[tuple[int, int]] = set()
    derivative_specialization = Fraction(0)

    for row in spine["coefficient_rows"]:
        occupied = int(row["occupied"])
        rank_black = int(row["rank_black"])
        rank_white = int(row["rank_white"])
        count = int(row["count"])
        if rank_black + rank_white != 2:
            raise AssertionError("digital-Alexander rank sum failed")
        q_numerator = rank_black - rank_white
        if q_numerator % 2:
            raise AssertionError("rank difference is not an even charge")
        q = q_numerator // 2
        if q not in sectors:
            raise AssertionError("charge escaped {-1,0,+1}")
        weight = count * p**occupied * (1 - p) ** (N - occupied)
        sectors[q] += weight
        terminal_support.add((rank_black, rank_white))
        derivative_specialization += Fraction(rank_black - rank_white, 2) * weight

    if sum(sectors.values()) != 1:
        raise AssertionError("sector probabilities do not normalize")
    if terminal_support != {(0, 2), (1, 1), (2, 0)}:
        raise AssertionError("terminal rank support did not collapse to three monomials")

    mean = sectors[1] - sectors[-1]
    raw_second = sectors[1] + sectors[-1]
    variance = raw_second - mean**2
    covariance = ((variance, -variance), (-variance, variance))
    determinant = covariance[0][0] * covariance[1][1] - covariance[0][1] * covariance[1][0]
    third_cumulant_direct = sum(
        probability * (Fraction(q) - mean) ** 3
        for q, probability in sectors.items()
    )
    third_cumulant_closure = mean - 3 * mean * variance - mean**3
    reconstructed = {
        1: (variance + mean**2 + mean) / 2,
        -1: (variance + mean**2 - mean) / 2,
        0: 1 - variance - mean**2,
    }
    if reconstructed != sectors:
        raise AssertionError("mean/variance sector inversion failed")
    if derivative_specialization != mean:
        raise AssertionError("relative source and bivariate derivative disagree")
    if third_cumulant_direct != third_cumulant_closure:
        raise AssertionError("three-state cumulant closure failed")

    bernstein_sector_counts = {
        q: [0] * (N + 1) for q in (-1, 0, 1)
    }
    for row in spine["coefficient_rows"]:
        q = (int(row["rank_black"]) - int(row["rank_white"])) // 2
        bernstein_sector_counts[q][int(row["occupied"])] += int(row["count"])

    return {
        "schema": "matching-one.p144-relative-rank-source-collapse.v1",
        "issues": [144, 269, 54],
        "geometry": spine["geometry"],
        "N": N,
        "p": fraction_record(p),
        "exact_factorization": {
            "Phi": "x*y*Z_rel(p,x/y)",
            "Z_rel": "P_minus*Q^-1+P_zero+P_plus*Q",
            "terminal_monomials": ["y^2", "x*y", "x^2"],
            "terminal_rank_pairs": [list(pair) for pair in sorted(terminal_support)],
        },
        "sector_probabilities": {
            "P_minus": fraction_record(sectors[-1]),
            "P_zero": fraction_record(sectors[0]),
            "P_plus": fraction_record(sectors[1]),
        },
        "terminal_Bernstein_counts_by_q": {
            str(q): values for q, values in bernstein_sector_counts.items()
        },
        "relative_source": {
            "mean_q_equals_M": fraction_record(mean),
            "raw_second_q": fraction_record(raw_second),
            "variance_q": fraction_record(variance),
            "matching_from_bivariate_derivative": fraction_record(derivative_specialization),
        },
        "diagonal_source": {
            "rank_sum": 2,
            "operator_identity": "(x*d_x+y*d_y)^m Phi=2^m Phi",
            "stochastic_diagonal_cumulants_order_ge_2": 0,
        },
        "rank_covariance": {
            "matrix": [
                [fraction_record(value) for value in row] for row in covariance
            ],
            "determinant": fraction_record(determinant),
            "rank": 1 if variance > 0 else 0,
            "strict_rank_one_for_0_lt_p_lt_1": variance > 0,
        },
        "three_state_closure": {
            "support_identity": "q^3=q",
            "linear_source_PDE": "G_sss=G_s",
            "log_cumulant_identity": "F_sss=F_s-3*F_s*F_ss-F_s^3",
            "third_cumulant_direct": fraction_record(third_cumulant_direct),
            "third_cumulant_from_closure": fraction_record(third_cumulant_closure),
            "sector_inversion_from_mean_variance": {
                str(q): fraction_record(value) for q, value in reconstructed.items()
            },
        },
        "checks": {
            "rank_sum_two_for_every_coefficient_row": True,
            "terminal_support_exactly_three_monomials": True,
            "Phi_equals_xy_Zrel_symbolically_per_monomial": True,
            "relative_source_equals_matching_specialization": True,
            "covariance_is_strict_rank_one": variance > 0 and determinant == 0,
            "third_cumulant_closure_exact": True,
        },
        "boundary": "terminal output collapses to three sectors; the typed connectivity frontier does not",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--json",
        type=Path,
        default=Path("results/exact-typed-incidence-spine/relative-source-collapse.json"),
    )
    args = parser.parse_args()
    payload = axis_L3_relative_collapse()
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
