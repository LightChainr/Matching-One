#!/usr/bin/env python3
"""Exact priority-weighted theta/figure-eight direct-birth census.

This reuses the P334 subset homology oracle and the P337 universal-cover
carrier classifier.  It records directed insertion counts by predecessor size
and converts them to uniform-random-permutation probabilities with the exact
Beta weight 1 / (N * binom(N-1,k)).

The output is finite-volume topology, not an asymptotic arm-exponent fit.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from math import comb
from pathlib import Path
from typing import Any

try:
    from p334_birth_age_collision_review_20260830 import enumerate_states
    from p337_direct_birth_arm_topology import carrier_descriptor
except ModuleNotFoundError:
    from scripts.p334_birth_age_collision_review_20260830 import enumerate_states
    from scripts.p337_direct_birth_arm_topology import carrier_descriptor


DEFAULT_OUTPUT = Path("results/exact-direct-birth-carrier-census/latest.json")
GEOMETRIES = ((2, 1), (2, 2), (3, 0), (3, 1), (3, 2), (4, 0), (4, 1))
TYPES = ("one_carrier_theta", "two_carrier_figure_eight")


def priority_mass(counts: list[int], n: int) -> Fraction:
    return sum(
        (Fraction(count, n * comb(n - 1, k)) for k, count in enumerate(counts)),
        Fraction(),
    )


def carrier_census(a: int, b: int) -> dict[str, Any]:
    n, states = enumerate_states(a, b)
    counts = {name: [0] * n for name in TYPES}

    for mask, (old_rank, _) in enumerate(states):
        if old_rank != 0:
            continue
        k = mask.bit_count()
        for vertex in range(n):
            if (mask >> vertex) & 1:
                continue
            if states[mask | (1 << vertex)][0] != 2:
                continue
            descriptor = carrier_descriptor(a, b, mask, vertex)
            topology = descriptor["topology"]
            if topology not in counts:
                raise AssertionError(f"unclassified direct birth: {topology}")
            counts[topology][k] += 1

    masses = {name: priority_mass(row, n) for name, row in counts.items()}
    total_by_k = [sum(counts[name][k] for name in TYPES) for k in range(n)]
    total_mass = sum(masses.values(), Fraction())
    rows = []
    for k in range(n):
        if total_by_k[k] == 0:
            continue
        weight = Fraction(1, n * comb(n - 1, k))
        rows.append(
            {
                "k": k,
                "priority_weight": str(weight),
                "theta_edges": counts["one_carrier_theta"][k],
                "figure_eight_edges": counts["two_carrier_figure_eight"][k],
                "direct_edges": total_by_k[k],
                "theta_mass_contribution": str(counts["one_carrier_theta"][k] * weight),
                "figure_eight_mass_contribution": str(
                    counts["two_carrier_figure_eight"][k] * weight
                ),
            }
        )

    theta_mass = masses["one_carrier_theta"]
    figure_mass = masses["two_carrier_figure_eight"]
    return {
        "a": a,
        "b": b,
        "N": n,
        "direct_edges": sum(total_by_k),
        "theta_edges": sum(counts["one_carrier_theta"]),
        "figure_eight_edges": sum(counts["two_carrier_figure_eight"]),
        "theta_probability": str(theta_mass),
        "theta_probability_float": float(theta_mass),
        "figure_eight_probability": str(figure_mass),
        "figure_eight_probability_float": float(figure_mass),
        "direct_birth_probability": str(total_mass),
        "direct_birth_probability_float": float(total_mass),
        "figure_eight_share_of_direct_probability": (
            str(figure_mass / total_mass) if total_mass else None
        ),
        "rows_by_predecessor_size": rows,
    }


def build_result() -> dict[str, Any]:
    rows = [carrier_census(a, b) for a, b in GEOMETRIES]
    expected_edges = {
        (2, 1): (0, 0, 0),
        (2, 2): (40, 40, 0),
        (3, 0): (45, 36, 9),
        (3, 1): (80, 80, 0),
        (3, 2): (793, 793, 0),
        (4, 0): (4624, 4288, 336),
        (4, 1): (8823, 8704, 119),
    }
    expected_total_mass = {
        (2, 1): "0",
        (3, 0): "3/35",
        (3, 1): "5/63",
        (3, 2): "304/3465",
        (4, 0): "2809/45045",
        (4, 1): "3511/60060",
    }
    for row in rows:
        key = (row["a"], row["b"])
        if (
            row["direct_edges"],
            row["theta_edges"],
            row["figure_eight_edges"],
        ) != expected_edges[key]:
            raise AssertionError(f"carrier edge reference changed for {key}")
        if key in expected_total_mass and row["direct_birth_probability"] != expected_total_mass[key]:
            raise AssertionError(f"priority-weighted total changed for {key}")

    return {
        "schema": "matching-one/p405-direct-birth-carrier-census/v1",
        "status": "exact_finite_volume_priority_weighted",
        "identity": (
            "D_c(N)=sum_k E_c(k)/(N*binom(N-1,k))="
            "sum_v integral_0^1 P_p(E_c(v,A)) dp"
        ),
        "geometries": rows,
        "decision": "theta_and_figure_eight_priority_masses_separated_exactly",
        "claim_boundary": [
            "The census is finite-volume and does not fit or prove either arm exponent.",
            "Raw directed-edge shares differ from uniform-permutation probability shares because predecessor sizes carry unequal Beta weights.",
            "Zero figure-eight mass on one finite quotient is not a universal absence theorem.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--stdout", action="store_true")
    args = parser.parse_args()
    result = build_result()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    if args.stdout:
        print(text, end="")


if __name__ == "__main__":
    main()
