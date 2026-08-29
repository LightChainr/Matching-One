#!/usr/bin/env python3
"""Minimal lattice Q-score interface for the #263 boundary CLE tangent."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Iterable, Sequence

from p263_boundary_tangent_ode import frobenius_high_branch


LINK_PATTERNS = ("1234", "12|34", "14|23")
LAMBDAS = (Fraction(1, 4), Fraction(1, 3), Fraction(2, 3), Fraction(3, 4))
ANCHOR_INDEX = 1
SERIES_ORDER = 100
CHECK_ORDER = 140


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _find(parent: list[int], vertex: int) -> int:
    if parent[vertex] != vertex:
        parent[vertex] = _find(parent, parent[vertex])
    return parent[vertex]


def terminal_partition(mask: int) -> str:
    """Partition of four boundary vertices on the open square cycle."""

    edges = ((0, 1), (1, 2), (2, 3), (3, 0))
    parent = list(range(4))
    for index, (first, second) in enumerate(edges):
        if not ((mask >> index) & 1):
            continue
        root_first = _find(parent, first)
        root_second = _find(parent, second)
        if root_first != root_second:
            parent[root_second] = root_first
    blocks: dict[int, list[int]] = {}
    for vertex in range(4):
        blocks.setdefault(_find(parent, vertex), []).append(vertex + 1)
    return "|".join(
        "".join(str(vertex) for vertex in block)
        for block in sorted(blocks.values(), key=lambda value: value[0])
    )


def square_cycle_rows() -> list[dict]:
    rows = []
    for mask in range(16):
        partition = terminal_partition(mask)
        cluster_count = len(partition.split("|"))
        bonds = mask.bit_count()
        rows.append(
            {
                "mask": mask,
                "partition": partition,
                "b": bonds,
                "k": cluster_count,
                "J": 2 * cluster_count + bonds,
            }
        )
    return rows


def sufficient_statistics(rows: Sequence[dict]) -> dict:
    output = {
        "samples": len(rows),
        "sum_J": sum(row["J"] for row in rows),
        "sum_J2": sum(row["J"] ** 2 for row in rows),
        "link_patterns": {},
    }
    for pattern in LINK_PATTERNS:
        selected = [row for row in rows if row["partition"] == pattern]
        output["link_patterns"][pattern] = {
            "count": len(selected),
            "sum_J": sum(row["J"] for row in selected),
            "sum_J2": sum(row["J"] ** 2 for row in selected),
            "sum_explicit_projector_field_Q_derivative": "0",
        }
    return output


def measure_tangent_from_sufficient(statistics: dict, pattern: str) -> Fraction:
    """Cov(1_pattern,J/2) at Q=1 on v=sqrt(Q)."""

    samples = statistics["samples"]
    channel = statistics["link_patterns"][pattern]
    return Fraction(channel["sum_J"], 2 * samples) - Fraction(
        channel["count"] * statistics["sum_J"], 2 * samples * samples
    )


def direct_polynomial_tangent(rows: Sequence[dict], pattern: str) -> Fraction:
    """Differentiate sum I*x^J / sum x^J, then use d/dQ=(1/2)d/dx."""

    numerator = Counter(row["J"] for row in rows if row["partition"] == pattern)
    denominator = Counter(row["J"] for row in rows)
    n0 = sum(numerator.values())
    z0 = sum(denominator.values())
    nprime = sum(power * count for power, count in numerator.items())
    zprime = sum(power * count for power, count in denominator.items())
    return Fraction(nprime * z0 - n0 * zprime, 2 * z0 * z0)


def _high_branch_unit_shape(lam: float, order: int) -> float:
    """(pi/sqrt(3))*partial_Q log V_{3h+1} at Q=1, modulo amplitude."""

    _, coefficients = frobenius_high_branch(order)
    ordinary = sum(float(coefficient.value) * lam**n for n, coefficient in enumerate(coefficients))
    # d_Q c_n = (-3/2 d_kappa c_n) * sqrt(3)/pi.
    regular_tangent_unit = sum(
        float(-Fraction(3, 2) * coefficient.derivative) * lam**n
        for n, coefficient in enumerate(coefficients)
    )
    return math.log(lam) + regular_tangent_unit / ordinary


def high_branch_targets() -> dict:
    unit = [_high_branch_unit_shape(float(value), SERIES_ORDER) for value in LAMBDAS]
    check = [_high_branch_unit_shape(float(value), CHECK_ORDER) for value in LAMBDAS]
    anchor = unit[ANCHOR_INDEX]
    gauge_free_unit = [value - anchor for value in unit]
    scale = math.sqrt(3) / math.pi
    actual = [scale * value for value in gauge_free_unit]
    return {
        "channel": "14|23",
        "continuum_identification": "U^(14)(23)(lambda)=C1*V_(3h+1)(lambda)",
        "lambda_order": [fraction_text(value) for value in LAMBDAS],
        "anchor_lambda": fraction_text(LAMBDAS[ANCHOR_INDEX]),
        "series_order": SERIES_ORDER,
        "independent_check_order": CHECK_ORDER,
        "max_abs_order_difference": max(abs(a - b) for a, b in zip(unit, check)),
        "anchored_dQ_logU_in_sqrt3_over_pi_units": gauge_free_unit,
        "anchored_dQ_logU": actual,
        "reflected_channel": "U^(12)(34)(lambda)=C1*V_(3h+1)(1-lambda)",
    }


def boundary_geometry(lam: Fraction) -> dict:
    """Rational four-point representative with the requested cross-ratio."""

    # x=(0,s,1,2), lambda=s/(2-s), hence s=2lambda/(1+lambda).
    s = 2 * lam / (1 + lam)
    k_prefactor = (2 - s) / (2 * s * (1 - s))
    return {
        "lambda": fraction_text(lam),
        "normalized_boundary_points": ["0", fraction_text(s), "1", "2"],
        "K": fraction_text(k_prefactor),
        "K_decimal": float(k_prefactor),
    }


def render() -> dict:
    rows = square_cycle_rows()
    sufficient = sufficient_statistics(rows)
    derivatives = {}
    for pattern in LINK_PATTERNS:
        score = measure_tangent_from_sufficient(sufficient, pattern)
        direct = direct_polynomial_tangent(rows, pattern)
        if score != direct:
            raise AssertionError("measure score did not match direct FK differentiation")
        derivatives[pattern] = {
            "probability_Q1": fraction_text(
                Fraction(sufficient["link_patterns"][pattern]["count"], len(rows))
            ),
            "measure_score_covariance": fraction_text(score),
            "direct_critical_manifold_Q_derivative": fraction_text(direct),
            "explicit_projector_field_derivative": "0",
        }
    if sum(Fraction(row["measure_score_covariance"]) for row in derivatives.values()) != Fraction(-39, 256):
        raise AssertionError("unexpected three-link total tangent")

    return {
        "schema": "matching-one.p263-boundary-lattice-qscore.v1",
        "issue": 263,
        "status": "frozen_minimal_lattice_interface",
        "continuum_source": {
            "paper": "Gefei Cai, Boundary four-point connectivities of conformal loop ensembles",
            "arxiv": "2603.28161v2",
            "lattice_link_patterns_equation": "Eq. (1.1)",
            "conformal_prefactor_equation": "Eq. (1.4)",
            "branch_identification": "Theorem 1.4",
        },
        "boundary_observable": {
            "boundary_condition": "free",
            "marked_vertices": "four distinct boundary vertices x1<x2<x3<x4",
            "primary_vector_order": list(LINK_PATTERNS),
            "definitions": {
                "1234": "all four marked vertices lie in one FK cluster",
                "12|34": "x1,x2 share one cluster; x3,x4 share a distinct cluster",
                "14|23": "x1,x4 share one cluster; x2,x3 share a distinct cluster",
            },
            "planar_forbidden": "13|24 as two distinct clusters",
            "other_terminal_partitions": "recorded for normalization/debugging but zero in the primary three-vector",
        },
        "frozen_geometries": [boundary_geometry(value) for value in LAMBDAS],
        "q_score_decomposition": {
            "critical_square_FK_manifold": "v=sqrt(Q), p=v/(1+v)",
            "measure_score": "S_measure=(J-E[J])/2 with J=2k+b; d_Q E[O]|_1=Cov(O,J/2)",
            "explicit_projector_derivative": "zero for the three plain connectivity indicators",
            "explicit_field_derivative": "zero for the bare indicators; renormalized boundary fields add 4*h_prime*log(L)*G plus an unknown lambda-independent normalization tangent",
            "h_prime_at_Q1": "sqrt(3)/(3*pi)",
            "conformal_prefactor_derivative": "2*h_prime*log(K)*G, kept separate and subtracted when forming d_Q log U",
            "no_double_count_rule": "use normalized continuum x_i and lattice spacing delta=1/L; log(L) belongs to field renormalization and log(K) to Eq. (1.4)",
        },
        "batch_sufficient_statistics_schema": {
            "required_integer_fields": [
                "geometry_id",
                "batch",
                "samples",
                "sum_J",
                "sum_J2",
                "for_each_link_pattern: count,sum_J,sum_J2",
            ],
            "J": "2*number_of_FK_clusters + number_of_open_bonds",
            "covariance_rule": "recompute probabilities and all three measure-score covariances inside each synchronized delete-one replicate",
            "cross_geometry_rule": "share counter domains where geometry coupling is declared; otherwise retain block-diagonal covariance rather than inventing pairing",
        },
        "amplitude_gauge_and_score": {
            "raw_lattice_log_tangent": "z_i=dP_i/P_i + 4*h_prime*log(L_i) + explicit_i/P_i - 2*h_prime*log(K_i)",
            "amplitude_gauge": "z_i -> z_i+c from the Q derivative of the common boundary-field normalization C1(Q)",
            "projection": "subtract z at lambda=1/3; equivalently GLS-project the raw tangent off the ordinary U0 vector",
            "primary_channel": "14|23",
            "frozen_target": high_branch_targets(),
            "residual_vector": "three non-anchor components of anchored z_lattice-z_CLE",
            "primary_score": "r^T Sigma_r^+ r with df=rank(Sigma_r)",
            "secondary_crossing": "anchored 12|34 tangent at lambda equals anchored 14|23 tangent at 1-lambda on the symmetric frozen grid",
        },
        "tiny_square_bond_regression": {
            "graph": "open four-cycle with every vertex a marked boundary point",
            "configurations": len(rows),
            "all_terminal_partition_counts": dict(sorted(Counter(row["partition"] for row in rows).items())),
            "sufficient_statistics": sufficient,
            "derivatives": derivatives,
            "three_link_total_probability_Q1": "7/16",
            "three_link_total_tangent": "-39/256",
            "direct_polynomial_equals_measure_score": True,
        },
        "integer_Q_fallback": {
            "generic_Q_sampler_required_for_plain_connectivities": False,
            "primary_replacement": "sample only Q=1 FK bond configurations and use the exact J/2 measure score",
            "projector_validation_if_later_needed": "at integer Q=2,3,4 color each FK cluster uniformly, regress the declared color-projector polynomial on the terminal partition, then differentiate that symbolic polynomial at Q=1; never finite-difference noisy Q runs",
        },
        "claim_layers": {
            "exact_finite_lattice": [
                "Q=1 critical-manifold measure score J/2",
                "three link-pattern definitions and tiny four-cycle regression",
                "separation of measure, bare/projector, field-renormalization and conformal-prefactor terms",
            ],
            "exact_continuum_from_source_and_parent_oracle": [
                "14|23 is the high V_(3h+1) branch",
                "12|34 is its lambda->1-lambda reflection",
                "the inhomogeneous tangent ODE and universal sqrt(3)/pi exponent derivative",
            ],
            "scaling_hypothesis": "finite square-FK boundary probabilities converge after L^(4h) renormalization to Cai's CLE connectivity vector",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = render()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
