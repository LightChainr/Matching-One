#!/usr/bin/env python3
"""Issue #180: exact rank-2 clock algebra and N145/N290 extrapolation."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Mapping, Sequence

from analyze_p48_retrospective import covariance_of_mean
from score_p50_fullcurve_n290 import (
    grouped,
    load_metadata,
    read_one_size,
    rng_group,
    sha256,
    size_statistics,
)


STATE_ORDER = ("I_S", "I_Du", "T_D", "T_Su")


def matmul2(a: Sequence[Sequence[Fraction]], b: Sequence[Sequence[Fraction]]):
    return [
        [sum((a[i][k] * b[k][j] for k in range(2)), Fraction(0)) for j in range(2)]
        for i in range(2)
    ]


def affine_matrix2(
    identity_coefficient: Fraction,
    transfer_coefficient: Fraction,
    transfer: Sequence[Sequence[Fraction]],
):
    return [
        [
            identity_coefficient * Fraction(i == j) + transfer_coefficient * transfer[i][j]
            for j in range(2)
        ]
        for i in range(2)
    ]


def square_clock_coefficients(trace: Fraction, determinant: Fraction) -> tuple[Fraction, Fraction]:
    """A^2=alpha I+beta A from Cayley--Hamilton."""

    return -determinant, trace


def state_from_statistics(stat: Mapping[str, float], n: int) -> list[float]:
    slope = float(stat["mean_slope"])
    if not math.isfinite(slope) or slope == 0:
        raise ValueError(f"N={n}: invalid mean slope")
    n13 = float(n) ** (13.0 / 8.0)
    return [
        float(n) * float(stat["P4_S"]),
        float(n) * float(stat["P4_D_prime"]) / slope,
        n13 * float(stat["P4_D"]),
        n13 * float(stat["P4_S_prime"]) / slope,
    ]


def estimate_state(data, n: int) -> tuple[list[float], list[list[float]]]:
    by_orientation = grouped(data, n)
    point = state_from_statistics(size_statistics(by_orientation, lineage_sign=1.0), n)
    batch_ids = [row.batch for row in by_orientation["first"]]
    deleted = [
        state_from_statistics(
            size_statistics(by_orientation, lineage_sign=1.0, omitted=batch), n
        )
        for batch in batch_ids
    ]
    batches = len(deleted)
    pseudo = [
        [batches * point[j] - (batches - 1) * row[j] for j in range(len(point))]
        for row in deleted
    ]
    return point, covariance_of_mean(pseudo)


def affine_predict(
    parent: Sequence[float], child: Sequence[float], alpha: float, beta: float
) -> list[float]:
    return [alpha * float(x) + beta * float(y) for x, y in zip(parent, child)]


def affine_covariance(
    parent_cov: Sequence[Sequence[float]],
    child_cov: Sequence[Sequence[float]],
    alpha: float,
    beta: float,
) -> list[list[float]]:
    """Prediction covariance for independent parent and child streams."""

    return [
        [
            alpha * alpha * float(parent_cov[i][j])
            + beta * beta * float(child_cov[i][j])
            for j in range(len(parent_cov))
        ]
        for i in range(len(parent_cov))
    ]


def exact_oracle() -> dict[str, object]:
    ordinary = [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(1, 2)]]
    jordan = [[Fraction(1), Fraction(1)], [Fraction(0), Fraction(1)]]
    models = {}
    for name, transfer, trace, determinant in (
        ("ordinary_q2", ordinary, Fraction(3, 2), Fraction(1, 2)),
        ("rank2_Jordan", jordan, Fraction(2), Fraction(1)),
    ):
        alpha, beta = square_clock_coefficients(trace, determinant)
        direct = matmul2(transfer, transfer)
        affine = affine_matrix2(alpha, beta, transfer)
        assert direct == affine
        models[name] = {
            "alpha": str(alpha),
            "beta": str(beta),
            "A": [[str(v) for v in row] for row in transfer],
            "A_squared": [[str(v) for v in row] for row in direct],
        }

    # The exact #249 witness has endpoint A eigenvalues 1,2 and U=A^2 there,
    # hence U=-2 I+3 A.  Its charged entries are A=-1,U=0, which reject the
    # endpoint affine pencil without changing the endpoint rank.
    charged_predicted = Fraction(-2) + Fraction(3) * Fraction(-1)
    charged_observed = Fraction(0)
    return {
        "models": models,
        "enriched_sector_counterexample": {
            "A_charged": "-1",
            "U_charged": "0",
            "endpoint_alpha": "-2",
            "endpoint_beta": "3",
            "endpoint_affine_prediction": str(charged_predicted),
            "defect": str(charged_observed - charged_predicted),
        },
    }


def render(
    parent_hist: Path,
    child_hist: Path,
    parent_meta_path: Path,
    child_meta_path: Path,
) -> dict[str, object]:
    parent_meta = load_metadata(parent_meta_path)
    child_meta = load_metadata(child_meta_path)
    if rng_group(parent_meta) == rng_group(child_meta):
        raise ValueError("N145 and N290 must be independent streams")
    parent, parent_cov = estimate_state(read_one_size(parent_hist, 145), 145)
    child, child_cov = estimate_state(read_one_size(child_hist, 290), 290)

    specifications = {
        "ordinary_q2": (-0.5, 1.5),
        "rank2_Jordan": (-1.0, 2.0),
    }
    predictions = {}
    for name, (alpha, beta) in specifications.items():
        covariance = affine_covariance(parent_cov, child_cov, alpha, beta)
        predictions[name] = {
            "identity_coefficient_alpha": alpha,
            "one_step_coefficient_beta": beta,
            "N580_state_prediction": dict(
                zip(STATE_ORDER, affine_predict(parent, child, alpha, beta))
            ),
            "N580_state_prediction_covariance": covariance,
            "N580_state_prediction_standard_error": dict(
                zip(STATE_ORDER, [math.sqrt(max(0.0, covariance[i][i])) for i in range(4)])
            ),
        }

    return {
        "schema": "matching-one/p180-affine-clock-hankel/v1",
        "issue": 180,
        "status": "exact_rank2_clock_theorem_and_post_reveal_N580_prediction",
        "state_order": list(STATE_ORDER),
        "observed_state": {
            "N145": dict(zip(STATE_ORDER, parent)),
            "N290": dict(zip(STATE_ORDER, child)),
        },
        "predictions": predictions,
        "exact_oracle": exact_oracle(),
        "hankel_gate": {
            "identity": "H_U=alpha*H_I+beta*H_A for every source/readout context",
            "rank": "rank([vec(H_I),vec(H_A),vec(H_U)])<=2",
            "commutator": "U=alpha*I+beta*A implies [A,U]=0",
            "minimal_score": "solve alpha,beta on two noncollinear context entries and predict every remaining entry",
        },
        "claim_boundary": {
            "exact": "Cayley-Hamilton coefficients and the shared-context affine/commutator implications",
            "post_reveal": "N580 state predictions use the already opened N145/N290 means",
            "conditional": "all four lattice readouts must factor through one common two-state transfer realization",
        },
        "provenance": {
            "parent_histogram": str(parent_hist),
            "parent_sha256": sha256(parent_hist),
            "child_histogram": str(child_hist),
            "child_sha256": sha256(child_hist),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--parent-hist",
        type=Path,
        default=Path("results/server-20260829/P50-n145-n290-fullcurve/raw/n145_100m.hist.csv"),
    )
    parser.add_argument(
        "--child-hist",
        type=Path,
        default=Path("results/server-20260829/P50-n145-n290-fullcurve/raw/n290_100m.hist.csv"),
    )
    parser.add_argument(
        "--parent-meta",
        type=Path,
        default=Path("results/server-20260829/P50-n145-n290-fullcurve/raw/n145_100m.metadata.json"),
    )
    parser.add_argument(
        "--child-meta",
        type=Path,
        default=Path("results/server-20260829/P50-n145-n290-fullcurve/raw/n290_100m.metadata.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/full-curve-transfer/p180_n145_n290_affine_clock.json"),
    )
    args = parser.parse_args()
    payload = render(args.parent_hist, args.child_hist, args.parent_meta, args.child_meta)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
