#!/usr/bin/env python3
"""Correct the Issue #43 matching-even score through the exact channel map."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Sequence

from wrapping_channels import (
    Combination,
    Normalization,
    ObservableDescriptor,
    OrientationOrder,
    ProbabilityCoordinate,
    Quantity,
    TopologyChannel,
    map_observable,
)


PREDICTION_SHA256 = "a370e79a10854341fac3ee75e8c518dbf3533e8c077cba2c2ec1018178144f44"
EXPECTED_PROTOCOL = "Issue #43 prospective N=185/265 two-spin4 full-curve score"

SOURCE_DESCRIPTOR = ObservableDescriptor(
    channel=TopologyChannel.EITHER,
    combination=Combination.EVEN,
    coordinate=ProbabilityCoordinate.P,
    orientation_order=OrientationOrder.FIRST_MINUS_SECOND,
    normalization=Normalization.RAW,
    quantity=Quantity.ORIENTATION_CONTRAST,
)
TARGET_DESCRIPTOR = ObservableDescriptor(
    channel=TopologyChannel.CROSS,
    combination=Combination.EVEN,
    coordinate=ProbabilityCoordinate.P,
    orientation_order=OrientationOrder.FIRST_MINUS_SECOND,
    normalization=Normalization.RAW,
    quantity=Quantity.ORIENTATION_CONTRAST,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def quadratic_2(vector: Sequence[float], covariance: Sequence[Sequence[float]]) -> float:
    a = float(covariance[0][0])
    b = float(covariance[0][1])
    d = float(covariance[1][1])
    determinant = a * d - b * b
    if determinant <= 0.0:
        raise ValueError("score covariance is not positive definite")
    x, y = map(float, vector)
    return (d * x * x - 2.0 * b * x * y + a * y * y) / determinant


def issue43_channel_map(prediction_text: str):
    """Validate source semantics before returning the registered exact map."""

    if "channel: either" not in prediction_text or "sector: even" not in prediction_text:
        raise ValueError("frozen artifact no longer declares the expected either/even source")
    return map_observable(SOURCE_DESCRIPTOR, TARGET_DESCRIPTOR)


def corrected_score(primary: dict, prediction_text: str) -> dict:
    transformation = issue43_channel_map(prediction_text)
    if primary.get("protocol") != EXPECTED_PROTOCOL:
        raise ValueError("input is not the Issue #43 primary score")
    if primary.get("prediction_artifact_sha256") != PREDICTION_SHA256:
        raise ValueError("primary score does not reference the frozen Issue #43 artifact")

    score = primary.get("scores", {}).get("DeltaS")
    if not isinstance(score, dict):
        raise ValueError("primary score lacks DeltaS")
    observed = [float(x) for x in score.get("observed", ())]
    sampling_se = [float(x) for x in score.get("sampling_se", ())]
    frozen_either = [float(x) for x in score.get("frozen_mean", ())]
    source_se = [float(x) for x in score.get("source_coefficient_se", ())]
    if not all(len(x) == 2 for x in (observed, sampling_se, frozen_either, source_se)):
        raise ValueError("Issue #43 correction requires the two N=185/265 endpoints")
    if any(x <= 0.0 or not math.isfinite(x) for x in sampling_se + source_se):
        raise ValueError("standard errors must be finite and positive")

    frozen_cross = [transformation.apply(x) for x in frozen_either]
    mapped_source_se = [transformation.apply_standard_error(x) for x in source_se]
    residual = [observed[i] - frozen_cross[i] for i in range(2)]
    covariance = [
        [
            (sampling_se[i] ** 2 if i == j else 0.0)
            + mapped_source_se[i] * mapped_source_se[j]
            for j in range(2)
        ]
        for i in range(2)
    ]
    chi_square = quadratic_2(residual, covariance)
    z = [residual[i] / math.sqrt(covariance[i][i]) for i in range(2)]

    return {
        "protocol": "Issue #43 exact cross/either channel-map correction",
        "status": "protocol correction; no target refit",
        "source_channel": "either/even",
        "target_channel": "cross/even",
        "source_descriptor": SOURCE_DESCRIPTOR.to_dict(),
        "target_descriptor": TARGET_DESCRIPTOR.to_dict(),
        "applied_transform": transformation.to_dict(),
        "exact_map": "DeltaS_cross = -DeltaS_either",
        "sizes": [185, 265],
        "observed_cross_DeltaS": observed,
        "sampling_se": sampling_se,
        "original_frozen_either_mean": frozen_either,
        "corrected_frozen_cross_mean": frozen_cross,
        "source_coefficient_se": mapped_source_se,
        "source_error_correlation": [[1.0, 1.0], [1.0, 1.0]],
        "residual": residual,
        "residual_covariance": covariance,
        "marginal_signed_z": z,
        "chi_square": chi_square,
        "df": 2,
        "chi_square_survival_df2": math.exp(-0.5 * chi_square),
        "target_refit_parameters": 0,
        "prediction_artifact_sha256": PREDICTION_SHA256,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--prediction", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()

    primary = json.loads(args.primary.read_text(encoding="utf-8"))
    prediction_text = args.prediction.read_text(encoding="utf-8")
    result = corrected_score(primary, prediction_text)
    if sha256(args.prediction) != PREDICTION_SHA256:
        raise SystemExit("frozen prediction artifact SHA-256 mismatch")
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
