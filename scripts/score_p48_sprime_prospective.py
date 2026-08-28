#!/usr/bin/env python3
"""Score fresh N=185/265 P48 S-prime targets against the training-only freeze.

The frozen source model lives in predictions/p48_sprime_correction_20260828.yaml
and uses only N=65,85,130 for coefficient estimation. This scorer never fits a
parameter to N=145,170,185,265.

Input target JSON is in the raw projector convention::

    {
      "sizes": [185, 265],
      "independent_of_training_source": true,
      "P4_S_prime": [value185, value265],
      "covariance_P4_S_prime": [[var185, cov], [cov, var265]]
    }

The score is performed in Y_N=N^(5/4) P4[S_prime] space. Target sampling
covariance is added to the frozen source-prediction covariance. All reported
models have two target degrees of freedom and zero target-fit parameters.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import yaml


TARGET_SIZES = (185, 265)
POWER = 1.25
MODEL_ORDER = (
    "pure_power_baseline",
    "rank2_log_primary_correction",
    "analytic_inverse_N_competitor",
    "zero_effect",
)


def read_json(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("target must be a JSON object")
    return payload


def read_yaml(path: Path) -> dict:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("prediction artifact must be a YAML mapping")
    return payload


def validate_covariance(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    if len(matrix) != 2 or any(len(row) != 2 for row in matrix):
        raise ValueError("covariance must be 2x2")
    result = [[float(value) for value in row] for row in matrix]
    if any(not math.isfinite(value) for row in result for value in row):
        raise ValueError("covariance entries must be finite")
    if not math.isclose(result[0][1], result[1][0], rel_tol=1e-12, abs_tol=1e-18):
        raise ValueError("covariance must be symmetric")
    if result[0][0] <= 0.0 or result[1][1] <= 0.0:
        raise ValueError("covariance diagonal must be positive")
    if result[0][0] * result[1][1] - result[0][1] ** 2 <= 0.0:
        raise ValueError("covariance must be positive definite")
    return result


def validate_artifact(artifact: Mapping[str, object]) -> None:
    if artifact.get("status") != "frozen_before_target_reveal":
        raise ValueError("artifact is not frozen_before_target_reveal")
    training = artifact.get("training_only")
    if not isinstance(training, dict):
        raise ValueError("artifact lacks training_only block")
    if tuple(int(value) for value in training.get("sizes", ())) != (65, 85, 130):
        raise ValueError("training sizes changed")
    if tuple(int(value) for value in training.get("excluded_from_fit", ())) != (145, 170):
        raise ValueError("N145/N170 exclusion contract changed")
    models = artifact.get("models")
    if not isinstance(models, dict):
        raise ValueError("artifact lacks models block")
    for name in MODEL_ORDER[:-1]:
        if name not in models:
            raise ValueError(f"artifact lacks frozen model {name}")


def validate_target(payload: Mapping[str, object]) -> tuple[list[float], list[list[float]]]:
    sizes = tuple(int(value) for value in payload.get("sizes", ()))
    if sizes != TARGET_SIZES:
        raise ValueError(f"target sizes must be exactly {TARGET_SIZES}")
    if payload.get("independent_of_training_source") is not True:
        raise ValueError("target must declare independent_of_training_source=true")
    observed = [float(value) for value in payload.get("P4_S_prime", ())]
    if len(observed) != 2 or any(not math.isfinite(value) for value in observed):
        raise ValueError("P4_S_prime must contain exactly two finite values")
    covariance = validate_covariance(payload.get("covariance_P4_S_prime", ()))
    return observed, covariance


def scale_target(
    observed_raw: Sequence[float], covariance_raw: Sequence[Sequence[float]]
) -> tuple[list[float], list[list[float]], list[float]]:
    scales = [float(n**POWER) for n in TARGET_SIZES]
    observed_y = [scales[i] * float(observed_raw[i]) for i in range(2)]
    covariance_y = [
        [scales[i] * float(covariance_raw[i][j]) * scales[j] for j in range(2)]
        for i in range(2)
    ]
    return observed_y, covariance_y, scales


def source_prediction(model_name: str, artifact: Mapping[str, object]) -> tuple[list[float], list[list[float]]]:
    if model_name == "zero_effect":
        return [0.0, 0.0], [[0.0, 0.0], [0.0, 0.0]]

    models = artifact["models"]
    model = models[model_name]
    targets = model["targets"]
    means = [float(targets[n]["mean"]) for n in TARGET_SIZES]

    if model_name == "pure_power_baseline":
        variance = float(model["A_se"]) ** 2
        covariance = [[variance, variance], [variance, variance]]
    else:
        covariance = [
            [float(value) for value in row]
            for row in model["target_source_covariance"]
        ]
        if len(covariance) != 2 or any(len(row) != 2 for row in covariance):
            raise ValueError(f"invalid source covariance for {model_name}")
    return means, covariance


def add_covariance(left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]) -> list[list[float]]:
    return [
        [float(left[i][j]) + float(right[i][j]) for j in range(2)]
        for i in range(2)
    ]


def quadratic_2(vector: Sequence[float], covariance: Sequence[Sequence[float]]) -> float:
    a = float(covariance[0][0])
    b = float(covariance[0][1])
    d = float(covariance[1][1])
    determinant = a * d - b * b
    if determinant <= 0.0:
        raise ValueError("residual covariance is not positive definite")
    x, y = map(float, vector)
    return (d * x * x - 2.0 * b * x * y + a * y * y) / determinant


def score(target: Mapping[str, object], artifact: Mapping[str, object]) -> dict:
    validate_artifact(artifact)
    observed_raw, covariance_raw = validate_target(target)
    observed_y, covariance_y, scales = scale_target(observed_raw, covariance_raw)

    rows = []
    for name in MODEL_ORDER:
        predicted_y, source_covariance = source_prediction(name, artifact)
        residual = [observed_y[i] - predicted_y[i] for i in range(2)]
        residual_covariance = add_covariance(covariance_y, source_covariance)
        chi_square = quadratic_2(residual, residual_covariance)
        rows.append(
            {
                "name": name,
                "predicted_Y": predicted_y,
                "predicted_P4_S_prime": [predicted_y[i] / scales[i] for i in range(2)],
                "source_prediction_covariance_Y": source_covariance,
                "residual_Y": residual,
                "residual_covariance_Y": residual_covariance,
                "chi_square": chi_square,
                "df": 2,
                "chi_square_survival_df2": math.exp(-0.5 * chi_square),
                "marginal_signed_z": [
                    residual[i] / math.sqrt(residual_covariance[i][i])
                    for i in range(2)
                ],
            }
        )

    return {
        "status": "frozen training-only prospective score; no target refit",
        "sizes": list(TARGET_SIZES),
        "observed_P4_S_prime": observed_raw,
        "observed_Y": observed_y,
        "target_covariance_Y": covariance_y,
        "model_order": list(MODEL_ORDER),
        "results": rows,
        "evidence_guard": (
            "The correction forms were motivated after the retrospective N145/N170 "
            "drift was known even though their coefficients exclude those sizes. "
            "Treat them as prospective for N185/N265, not as pre-drift confirmatory models."
        ),
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument(
        "--artifact",
        type=Path,
        default=root / "predictions/p48_sprime_correction_20260828.yaml",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = score(read_json(args.target), read_yaml(args.artifact))
    rendered = json.dumps(result, indent=2) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
