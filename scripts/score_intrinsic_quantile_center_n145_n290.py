#!/usr/bin/env python3
"""Score the frozen Issue #101 intrinsic quantile-center N145->N290 law.

This is a secondary post-primary reuse of the independent N145 and N290
full-curve blocks.  The frozen primary score in ``score_p50_fullcurve_n290``
is not changed.  Every size-local delete-one replicate resolves all four
level crossings at u={0.025,0.05} before Q and the scaled widths are formed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import yaml

from analyze_p48_retrospective import covariance_of_mean
from intrinsic_quantile_center import (
    DOUBLING_RATIO,
    FROZEN_U,
    midpoint_difference,
    quantile_levels,
)
from score_p49_fullcurve_doubling import aggregate, orientation_values
from score_p50_fullcurve_n290 import (
    CHILD_N,
    PARENT_N,
    grouped,
    load_metadata,
    read_one_size,
    rng_group,
)


FEATURE_ORDER = (
    "Q",
    "w_0.025_scaled",
    "w_0.05_scaled",
    "c_0.025",
    "c_0.05",
)
RESIDUAL_ORDER = (
    "Q290_minus_frozen_ratio_Q145",
    "scaled_width_drift_0.025",
    "scaled_width_drift_0.05",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def size_coordinates(by_orientation, omitted: int = -1) -> dict[str, object]:
    """Resolve the frozen levels for one full or delete-one size sample."""

    n = by_orientation["first"][0].n
    rows = {
        name: aggregate(by_orientation[name], omitted)
        for name in ("first", "second")
    }

    def mean_matching(p: float) -> float:
        return math.fsum(
            orientation_values(n, rows[name], p)["M"]
            for name in ("first", "second")
        ) / 2.0

    levels = quantile_levels(mean_matching)
    q_value = midpoint_difference(levels)
    payload: dict[str, object] = {
        "N": n,
        "Q": q_value,
        "Q_scaled_N_to_3_over_4": q_value * n**0.75,
        "levels": {},
    }
    level_payload = {}
    for u in FROZEN_U:
        level = levels[u]
        level_payload[str(u)] = {
            "p_minus": level.p_minus,
            "p_plus": level.p_plus,
            "c": level.c,
            "w": level.w,
            "w_scaled_N_to_3_over_8": level.w * n**0.375,
            "M_minus": mean_matching(level.p_minus),
            "M_plus": mean_matching(level.p_plus),
        }
    payload["levels"] = level_payload
    return payload


def feature_vector(coordinates: Mapping[str, object]) -> list[float]:
    levels = coordinates["levels"]
    return [
        float(coordinates["Q"]),
        float(levels["0.025"]["w_scaled_N_to_3_over_8"]),
        float(levels["0.05"]["w_scaled_N_to_3_over_8"]),
        float(levels["0.025"]["c"]),
        float(levels["0.05"]["c"]),
    ]


def pseudovalue_vectors(
    full: Sequence[float], deleted: Sequence[Sequence[float]]
) -> list[list[float]]:
    batches = len(deleted)
    if batches < 2 or any(len(row) != len(full) for row in deleted):
        raise ValueError("delete-one feature matrix is ragged or too short")
    return [
        [batches * full[j] - (batches - 1) * row[j] for j in range(len(full))]
        for row in deleted
    ]


def estimate_size(by_orientation):
    point = size_coordinates(by_orientation)
    batch_ids = [row.batch for row in by_orientation["first"]]
    deleted = [
        feature_vector(size_coordinates(by_orientation, omitted=batch))
        for batch in batch_ids
    ]
    pseudo = pseudovalue_vectors(feature_vector(point), deleted)
    return point, covariance_of_mean(pseudo)


def matmul(left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]):
    return [
        [
            math.fsum(left[i][k] * right[k][j] for k in range(len(right)))
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def transpose(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    return [list(column) for column in zip(*matrix)]


def transformed_covariance(transform, covariance):
    return matmul(matmul(transform, covariance), transpose(transform))


def add_matrices(first, second):
    return [
        [first[i][j] + second[i][j] for j in range(len(first[i]))]
        for i in range(len(first))
    ]


def render(
    parent_data,
    child_data,
    parent_metadata: Mapping[str, object],
    child_metadata: Mapping[str, object],
    freeze: Mapping[str, object],
    provenance: Mapping[str, object],
) -> dict[str, object]:
    if tuple(float(value) for value in freeze["frozen_u"]) != FROZEN_U:
        raise ValueError("prediction does not freeze exactly u={0.025,0.05}")
    target = float(freeze["predictions"]["doubling_ratio"]["float"])
    if not math.isclose(target, DOUBLING_RATIO, rel_tol=0.0, abs_tol=1e-15):
        raise ValueError("prediction doubling ratio differs from 2^{-3/4}")
    if rng_group(parent_metadata) == rng_group(child_metadata):
        raise ValueError("N145 and N290 must be independent random streams")

    parent, parent_cov = estimate_size(grouped(parent_data, PARENT_N))
    child, child_cov = estimate_size(grouped(child_data, CHILD_N))
    parent_vector = feature_vector(parent)
    child_vector = feature_vector(child)

    # Residuals are [Q_child-r*Q_parent, W025_child-W025_parent,
    # W05_child-W05_parent].  The two size blocks are independent.
    parent_transform = (
        (-target, 0.0, 0.0, 0.0, 0.0),
        (0.0, -1.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, -1.0, 0.0, 0.0),
    )
    child_transform = (
        (1.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 1.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0, 0.0, 0.0),
    )
    residual = [
        child_vector[0] - target * parent_vector[0],
        child_vector[1] - parent_vector[1],
        child_vector[2] - parent_vector[2],
    ]
    residual_cov = add_matrices(
        transformed_covariance(parent_transform, parent_cov),
        transformed_covariance(child_transform, child_cov),
    )
    marginal_se = [math.sqrt(max(0.0, residual_cov[i][i])) for i in range(3)]
    marginal_z = [
        residual[i] / marginal_se[i] if marginal_se[i] else None
        for i in range(3)
    ]

    q_parent = parent_vector[0]
    q_child = child_vector[0]
    observed_ratio = q_child / q_parent
    ratio_gradient_parent = -q_child / (q_parent * q_parent)
    ratio_gradient_child = 1.0 / q_parent
    ratio_variance = (
        ratio_gradient_parent**2 * parent_cov[0][0]
        + ratio_gradient_child**2 * child_cov[0][0]
    )

    return {
        "schema": "matching-one/intrinsic-quantile-center-N145-N290-score/v1",
        "status": "secondary post-primary score from pre-target frozen Issue #101 coordinates",
        "primary_score_unchanged": "scripts/score_p50_fullcurve_n290.py and its frozen order are untouched",
        "frozen": {
            "u": list(FROZEN_U),
            "Q": "c_0.05-c_0.025",
            "target_ratio": target,
            "target_formula": "Q_290/Q_145=2^{-3/4}",
        },
        "observations": {
            "N145": parent,
            "N290": child,
        },
        "size_local_feature_order": list(FEATURE_ORDER),
        "size_local_covariance": {
            "N145": parent_cov,
            "N290": child_cov,
            "cross_size": "zero by independent RNG domains",
        },
        "primary_quantile_center_score": {
            "observed_ratio_Q290_over_Q145": observed_ratio,
            "observed_ratio_standard_error_delta_method": math.sqrt(
                max(0.0, ratio_variance)
            ),
            "target_ratio": target,
            "residual_Q290_minus_target_Q145": residual[0],
            "residual_standard_error": marginal_se[0],
            "signed_z": marginal_z[0],
            "chi_square": marginal_z[0] ** 2 if marginal_z[0] is not None else None,
            "degrees_of_freedom": 1,
        },
        "width_drift": {
            "definition": "w_u N^{3/8}(N290)-w_u N^{3/8}(N145)",
            "u": list(FROZEN_U),
            "residual": residual[1:],
            "standard_error": marginal_se[1:],
            "signed_z": marginal_z[1:],
        },
        "joint_residual_order": list(RESIDUAL_ORDER),
        "joint_residual": residual,
        "joint_residual_covariance": residual_cov,
        "covariance_rule": (
            "all crossings are recomputed inside each size-local delete-one replicate; "
            "N145/N290 pseudovalue covariances are combined as independent blocks"
        ),
        "provenance": dict(provenance),
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-hist", type=Path, required=True)
    parser.add_argument("--parent-metadata", type=Path, required=True)
    parser.add_argument("--child-hist", type=Path, required=True)
    parser.add_argument("--child-metadata", type=Path, required=True)
    parser.add_argument(
        "--freeze",
        type=Path,
        default=root / "predictions/intrinsic_quantile_center_20260829.yaml",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    freeze = yaml.safe_load(args.freeze.read_text(encoding="utf-8"))
    if freeze.get("status") != "definition_freeze":
        raise ValueError("Issue #101 artifact is not the frozen definition contract")
    parent_metadata = load_metadata(args.parent_metadata)
    child_metadata = load_metadata(args.child_metadata)
    payload = render(
        read_one_size(args.parent_hist, PARENT_N),
        read_one_size(args.child_hist, CHILD_N),
        parent_metadata,
        child_metadata,
        freeze,
        {
            "freeze": str(args.freeze),
            "freeze_sha256": sha256(args.freeze),
            "freeze_commit": "3762b342b8b376e587df0073044b2c7f6452aa8e",
            "freeze_committed_at_utc": "2026-08-29T06:41:21Z",
            "target_results_first_commit": (
                "9675bce5b406247e15c03bca20abef954f26a3a2"
            ),
            "target_results_first_committed_at_utc": "2026-08-29T06:47:18Z",
            "chronology": "freeze precedes first target-result commit",
            "parent_hist": str(args.parent_hist),
            "parent_hist_sha256": sha256(args.parent_hist),
            "parent_metadata": str(args.parent_metadata),
            "parent_metadata_sha256": sha256(args.parent_metadata),
            "child_hist": str(args.child_hist),
            "child_hist_sha256": sha256(args.child_hist),
            "child_metadata": str(args.child_metadata),
            "child_metadata_sha256": sha256(args.child_metadata),
        },
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
