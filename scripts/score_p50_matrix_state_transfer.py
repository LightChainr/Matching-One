#!/usr/bin/env python3
"""Score the frozen low-rank state transfer on the N=145 -> 290 full curves.

This is a prospective *sensitivity* score for the same raw P50 full-curve block
used by the already-frozen center/slope/root scorer.  It must not be counted as
an additional independent evidence block.

The state is

    I_S  = N P4[S]
    I_Du = N P4[D'] / Mbar'
    T_D  = N^(13/8) P4[D]
    T_Su = N^(13/8) P4[S'] / Mbar'

and the target statistic is state(290)-state(145).  The first three components
are frozen to zero under the compact three-state closure.  The fourth component
uses either the source-fitted analytic A+C/N increment or the source-fitted
rank-2 Jordan A+B log N increment frozen in the prediction artifact.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import yaml

from analyze_p48_retrospective import covariance_of_mean, pseudovalues, quadratic
from score_p50_fullcurve_n290 import (
    grouped,
    load_metadata,
    read_one_size,
    rng_group,
    sha256,
    size_statistics,
)

PARENT_N = 145
CHILD_N = 290
STATE_ORDER = ("I_S", "I_Du", "T_D", "T_Su")


def state_from_statistics(stat: Mapping[str, float], n: int) -> list[float]:
    slope = float(stat["mean_slope"])
    if not math.isfinite(slope) or slope == 0.0:
        raise ValueError(f"N={n}: invalid intrinsic-center mean slope {slope}")
    n13 = math.pow(float(n), 13.0 / 8.0)
    return [
        float(n) * float(stat["P4_S"]),
        float(n) * float(stat["P4_D_prime"]) / slope,
        n13 * float(stat["P4_D"]),
        n13 * float(stat["P4_S_prime"]) / slope,
    ]


def estimate_state(data, n: int) -> tuple[list[float], list[list[float]]]:
    by_orientation = grouped(data, n)
    point_stat = size_statistics(by_orientation, lineage_sign=+1.0)
    point = state_from_statistics(point_stat, n)
    batch_ids = [row.batch for row in by_orientation["first"]]
    deleted = [
        state_from_statistics(
            size_statistics(by_orientation, lineage_sign=+1.0, omitted=batch), n
        )
        for batch in batch_ids
    ]
    pseudo = pseudovalues(point, deleted)
    covariance = covariance_of_mean(pseudo)
    return point, covariance


def add_covariances(
    left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]
) -> list[list[float]]:
    if len(left) != len(right):
        raise ValueError("covariance dimensions differ")
    return [
        [float(left[i][j]) + float(right[i][j]) for j in range(len(left))]
        for i in range(len(left))
    ]


def score_model(
    delta_state: Sequence[float],
    target_covariance: Sequence[Sequence[float]],
    model: Mapping[str, object],
) -> dict[str, object]:
    expected_map = model["expected_child_minus_parent"]
    expected = [float(expected_map[name]) for name in STATE_ORDER]
    residual = [float(delta_state[i]) - expected[i] for i in range(len(STATE_ORDER))]
    covariance = [list(map(float, row)) for row in target_covariance]
    source_se = float(model["source_se_T_Su_increment"])
    covariance[3][3] += source_se * source_se
    chi_square = quadratic(residual, covariance)
    t_su_z = residual[3] / math.sqrt(covariance[3][3])
    # Exact survival function for chi-square with four degrees of freedom.
    survival_df4 = math.exp(-chi_square / 2.0) * (1.0 + chi_square / 2.0)
    return {
        "expected_child_minus_parent": dict(zip(STATE_ORDER, expected)),
        "source_se_T_Su_increment": source_se,
        "residual": dict(zip(STATE_ORDER, residual)),
        "residual_covariance": covariance,
        "joint_chi_square": chi_square,
        "degrees_of_freedom": 4,
        "chi_square_survival_df4": survival_df4,
        "marginal_T_Su_signed_z": t_su_z,
    }


def render(
    parent_hist: Path,
    child_hist: Path,
    parent_meta_path: Path,
    child_meta_path: Path,
    prediction_path: Path,
) -> dict[str, object]:
    prediction = yaml.safe_load(prediction_path.read_text(encoding="utf-8"))
    if int(prediction["target_block"]["parent_N"]) != PARENT_N:
        raise ValueError("prediction parent N does not match scorer contract")
    if int(prediction["target_block"]["child_N"]) != CHILD_N:
        raise ValueError("prediction child N does not match scorer contract")
    if tuple(prediction["state_definition"]["state_order"]) != STATE_ORDER:
        raise ValueError("prediction state order does not match scorer contract")

    parent_meta = load_metadata(parent_meta_path)
    child_meta = load_metadata(child_meta_path)
    if rng_group(parent_meta) == rng_group(child_meta):
        raise ValueError(
            "matrix-state scorer assumes independent N145/N290 streams, but metadata share one RNG group"
        )

    parent_data = read_one_size(parent_hist, PARENT_N)
    child_data = read_one_size(child_hist, CHILD_N)
    parent_state, parent_cov = estimate_state(parent_data, PARENT_N)
    child_state, child_cov = estimate_state(child_data, CHILD_N)
    delta_state = [child_state[i] - parent_state[i] for i in range(len(STATE_ORDER))]
    target_covariance = add_covariances(parent_cov, child_cov)

    scores = {
        name: score_model(delta_state, target_covariance, model)
        for name, model in prediction["models"].items()
    }
    best = min(scores, key=lambda name: float(scores[name]["joint_chi_square"]))
    best_chi = float(scores[best]["joint_chi_square"])
    for name in scores:
        scores[name]["delta_chi_square_from_best"] = (
            float(scores[name]["joint_chi_square"]) - best_chi
        )

    return {
        "schema": "matching-one/P50-heldout-matrix-state-transfer/v1",
        "status": "prospective sensitivity score; nonadditive with P50 primary",
        "state_order": list(STATE_ORDER),
        "state": {
            "N145": dict(zip(STATE_ORDER, parent_state)),
            "N290": dict(zip(STATE_ORDER, child_state)),
            "child_minus_parent": dict(zip(STATE_ORDER, delta_state)),
        },
        "target_covariance_child_minus_parent": target_covariance,
        "scores": scores,
        "best_by_joint_chi_square": best,
        "evidence_rule": (
            "This score reuses the P50 N145/N290 raw block.  It is a mechanism sensitivity view, "
            "not an additional primary evidence vote."
        ),
        "provenance": {
            "prediction_artifact": str(prediction_path),
            "prediction_sha256": sha256(prediction_path),
            "parent_histogram": str(parent_hist),
            "parent_histogram_sha256": sha256(parent_hist),
            "child_histogram": str(child_hist),
            "child_histogram_sha256": sha256(child_hist),
            "parent_metadata": str(parent_meta_path),
            "parent_metadata_sha256": sha256(parent_meta_path),
            "child_metadata": str(child_meta_path),
            "child_metadata_sha256": sha256(child_meta_path),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
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
        "--prediction",
        type=Path,
        default=Path("predictions/p50_matrix_state_transfer_20260829.yaml"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/server-20260829/P50-n145-n290-fullcurve/matrix_state_transfer_score.json"),
    )
    args = parser.parse_args()

    payload = render(
        args.parent_hist,
        args.child_hist,
        args.parent_meta,
        args.child_meta,
        args.prediction,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
