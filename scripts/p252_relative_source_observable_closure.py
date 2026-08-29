#!/usr/bin/env python3
"""Exact observer closure for the three-sector relative topology source.

This extends the Issue #54 source PDE from the partition function to an
arbitrary single observer W.  It shows exactly which higher source rows are
new and which are algebraic repetitions.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Callable

from p54_relative_source_pde import (
    cumulant_engine,
    fraction_record,
    tiny_joint_distribution,
)


Distribution = list[tuple[int, int, Fraction]]
Observer = Callable[[int, int], int]


def exact_rank(matrix: list[list[Fraction]]) -> int:
    if not matrix:
        return 0
    rows = [row[:] for row in matrix]
    rank = 0
    columns = len(rows[0])
    for column in range(columns):
        pivot = next(
            (index for index in range(rank, len(rows)) if rows[index][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = rows[rank][column]
        rows[rank] = [value / scale for value in rows[rank]]
        for index in range(len(rows)):
            if index == rank or not rows[index][column]:
                continue
            factor = rows[index][column]
            rows[index] = [
                value - factor * pivot_value
                for value, pivot_value in zip(rows[index], rows[rank])
            ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def observer_record(
    distribution: Distribution,
    name: str,
    observer: Observer,
) -> dict:
    moments = [
        sum(weight * observer(q, score) * q**power for q, score, weight in distribution)
        for power in range(7)
    ]
    direct_sectors = {
        q_value: sum(
            weight * observer(q, score)
            for q, score, weight in distribution
            if q == q_value
        )
        for q_value in (-1, 0, 1)
    }
    reconstructed = {
        1: (moments[2] + moments[1]) / 2,
        -1: (moments[2] - moments[1]) / 2,
        0: moments[0] - moments[2],
    }

    modified_distribution = [
        (q, observer(q, score), weight) for q, score, weight in distribution
    ]
    kappa = cumulant_engine(modified_distribution)
    mu = sum(weight * q for q, _, weight in distribution)
    raw_second = sum(weight * q * q for q, _, weight in distribution)
    kappa_one = kappa(1, 1)
    kappa_two = kappa(2, 1)
    kappa_three = kappa(3, 1)
    closure = (1 - 3 * raw_second) * kappa_one - 3 * mu * kappa_two

    return {
        "name": name,
        "raw_moments_E_W_q_power_0_to_6": [fraction_record(value) for value in moments],
        "odd_raw_rows_equal": moments[1] == moments[3] == moments[5],
        "positive_even_raw_rows_equal": moments[2] == moments[4] == moments[6],
        "sector_weighted_observer_direct": {
            str(q): fraction_record(value) for q, value in direct_sectors.items()
        },
        "sector_weighted_observer_reconstructed": {
            str(q): fraction_record(value) for q, value in reconstructed.items()
        },
        "sector_reconstruction_exact": direct_sectors == reconstructed,
        "mixed_connected_cumulants": {
            "kappa_W_q": fraction_record(kappa_one),
            "kappa_W_q_q": fraction_record(kappa_two),
            "kappa_W_q_q_q": fraction_record(kappa_three),
            "third_from_closure": fraction_record(closure),
            "closure_exact": kappa_three == closure,
        },
    }


def build_oracle() -> dict:
    distribution = tiny_joint_distribution()
    observers: list[tuple[str, Observer]] = [
        ("one", lambda q, score: 1),
        ("thermal_score", lambda q, score: score),
        ("thermal_score_squared", lambda q, score: score * score),
        ("neutral_sector_mark", lambda q, score: int(q == 0)),
        ("signed_thermal_score", lambda q, score: q * score),
    ]
    records = [
        observer_record(distribution, name, observer)
        for name, observer in observers
    ]
    source_matrix = [
        [
            Fraction(row["numerator"], row["denominator"])
            for row in record["raw_moments_E_W_q_power_0_to_6"]
        ]
        for record in records
    ]

    return {
        "schema": "matching-one.p252-relative-source-observer-closure.v1",
        "issues": [54, 114, 252],
        "exact_algebra": {
            "algebra": "R[q]/(q^3-q)",
            "dimension": 3,
            "basis": ["1", "q", "q^2"],
            "primitive_sector_idempotents": {
                "plus": "(q^2+q)/2",
                "minus": "(q^2-q)/2",
                "zero": "1-q^2",
            },
            "observer_sector_inversion": {
                "W_plus": "(E[W q^2]+E[W q])/2",
                "W_minus": "(E[W q^2]-E[W q])/2",
                "W_zero": "E[W]-E[W q^2]",
            },
            "mixed_connected_closure": "kappa(W,q,q,q)=(1-3 E[q^2]) kappa(W,q)-3 E[q] kappa(W,q,q)",
        },
        "tiny_exact_oracle": {
            "geometry": "gaussian(2,1)",
            "N": 5,
            "p": "1/2",
            "observers": records,
            "observer_by_source_power_matrix_shape": [len(source_matrix), 7],
            "observer_by_source_power_matrix_rank": exact_rank(source_matrix),
            "rank_bound": 3,
        },
        "rank3_consequence": {
            "exact_no_go": "Repeated insertions of the single global relative-topology source cannot create a fourth independent linear response direction; every observer row is fixed by its 0-, 1-, and 2-source moments.",
            "issue_252": "A rank-3 second-energy/Jordan channel cannot be identified merely by taking a third relative-Q-source derivative. It requires an independent thermal, spatial, charged, or multi-point insertion W; after W is supplied, the relative source only sector-resolves that W.",
            "allowed": "The observer W itself may couple to a rank-3 Virasoro module. The theorem removes repeated global-source order as evidence for that coupling; it does not forbid the coupling.",
            "next_score": "Measure the independent W channel first, then store only E[W], E[Wq], and E[Wq^2] (or the equivalent three sector-weighted W means). Generate all higher relative-source rows algebraically.",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(build_oracle(), indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
