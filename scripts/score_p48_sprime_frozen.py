#!/usr/bin/env python3
"""Score fresh P48 P4[S'] targets against the frozen chronological manifest.

The target JSON must contain exactly N=185 and N=265 observations with their
2x2 sampling covariance. No model parameter is fit to target data.

This scorer assumes fresh target statistics are independent of the retrospective
P33/P48 source used to freeze amplitudes. If production deliberately reuses
source random counters, a covariance-aware replacement must be preregistered
before the target is inspected.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import List, Mapping, Sequence, Tuple

import yaml


def read_json(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return payload


def read_yaml(path: Path) -> dict:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a YAML mapping")
    return payload


def validate_covariance(
    matrix: Sequence[Sequence[float]], width: int
) -> List[List[float]]:
    if len(matrix) != width or any(len(row) != width for row in matrix):
        raise ValueError(f"covariance must be {width}x{width}")
    result = [[float(value) for value in row] for row in matrix]
    for i in range(width):
        if result[i][i] <= 0.0 or not math.isfinite(result[i][i]):
            raise ValueError("covariance diagonal must be finite and positive")
        for j in range(width):
            if not math.isfinite(result[i][j]):
                raise ValueError("covariance entries must be finite")
            tolerance = 1e-12 * max(
                1.0, abs(result[i][j]), abs(result[j][i])
            )
            if abs(result[i][j] - result[j][i]) > tolerance:
                raise ValueError("covariance must be symmetric")
    if width == 2:
        determinant = result[0][0] * result[1][1] - result[0][1] ** 2
        if determinant <= 0.0:
            raise ValueError("covariance must be positive definite")
    return result


def quadratic_2(
    vector: Sequence[float], covariance: Sequence[Sequence[float]]
) -> float:
    if len(vector) != 2:
        raise ValueError("this frozen scorer requires exactly two targets")
    a, b = map(float, covariance[0])
    _, d = map(float, covariance[1])
    determinant = a * d - b * b
    if determinant <= 0.0:
        raise ValueError("residual covariance is not positive definite")
    x, y = map(float, vector)
    return (d * x * x - 2.0 * b * x * y + a * y * y) / determinant


def matadd(
    left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]
) -> List[List[float]]:
    return [
        [float(left[i][j]) + float(right[i][j]) for j in range(2)]
        for i in range(2)
    ]


def prediction_covariance(
    design: Sequence[Sequence[float]],
    parameter_covariance: Sequence[Sequence[float]],
) -> List[List[float]]:
    width = len(parameter_covariance)
    if any(len(row) != width for row in parameter_covariance):
        raise ValueError("parameter covariance must be square")
    if any(len(row) != width for row in design):
        raise ValueError("design width does not match parameter covariance")
    result = [[0.0, 0.0], [0.0, 0.0]]
    for i in range(2):
        for j in range(2):
            result[i][j] = math.fsum(
                design[i][a]
                * float(parameter_covariance[a][b])
                * design[j][b]
                for a in range(width)
                for b in range(width)
            )
    return result


def basis(n: int, model: Mapping[str, object], power: float) -> List[float]:
    scale = n ** (-power)
    kind = str(model["basis"])
    if kind == "constant":
        return [scale]
    if kind == "constant_plus_inverse_N":
        return [scale, scale / n]
    if kind == "constant_plus_log_N":
        return [scale, scale * math.log(n)]
    if kind == "zero":
        return []
    raise ValueError(f"unknown frozen basis {kind!r}")


def validate_manifest(
    payload: Mapping[str, object]
) -> Tuple[Tuple[int, int], float, list]:
    if payload.get("status") != "frozen_prospective_scoring_manifest":
        raise ValueError("manifest is not frozen_prospective_scoring_manifest")
    sizes = tuple(int(value) for value in payload.get("target_sizes_N", ()))
    if len(sizes) != 2:
        raise ValueError("manifest must freeze exactly two target sizes")
    power = float(payload.get("leading_power_in_N"))
    if not math.isclose(power, 1.25, rel_tol=0.0, abs_tol=1e-15):
        raise ValueError("unexpected P4[S'] leading power")
    models = payload.get("models_in_scoring_order")
    if not isinstance(models, list) or len(models) < 2:
        raise ValueError("manifest must contain an ordered model list")
    names = [str(model.get("name")) for model in models]
    expected = [
        "pure_N^-5/4",
        "zero_effect",
        "q2_even_scalar",
        "rank2_jordan_log",
    ]
    if names != expected:
        raise ValueError(f"frozen scoring order changed: {names!r}")
    return (sizes[0], sizes[1]), power, models


def validate_target(
    payload: Mapping[str, object], target_sizes: Tuple[int, int]
) -> Tuple[List[float], List[List[float]]]:
    sizes = tuple(int(value) for value in payload.get("sizes", ()))
    if sizes != target_sizes:
        raise ValueError(
            f"fresh target sizes must be exactly {target_sizes}, got {sizes}"
        )
    if payload.get("independent_of_retrospective_source") is not True:
        raise ValueError(
            "target must explicitly set independent_of_retrospective_source=true"
        )
    observed = [float(value) for value in payload.get("P4_S_prime", ())]
    if len(observed) != 2 or any(not math.isfinite(value) for value in observed):
        raise ValueError("P4_S_prime must contain exactly two finite observations")
    covariance = validate_covariance(
        payload.get("covariance_P4_S_prime", ()), 2
    )
    return observed, covariance


def score_model(
    model: Mapping[str, object],
    target_sizes: Tuple[int, int],
    power: float,
    observed: Sequence[float],
    target_covariance: Sequence[Sequence[float]],
) -> dict:
    design = [basis(n, model, power) for n in target_sizes]
    parameters = [float(value) for value in model.get("parameters", ())]
    parameter_covariance = model.get("parameter_covariance", ())
    if model["basis"] == "zero":
        if parameters or parameter_covariance:
            raise ValueError("zero model must not carry source parameters")
        predicted = [0.0, 0.0]
        source_covariance = [[0.0, 0.0], [0.0, 0.0]]
    else:
        if any(len(row) != len(parameters) for row in parameter_covariance):
            raise ValueError("parameter covariance width mismatch")
        if len(parameter_covariance) != len(parameters):
            raise ValueError("parameter covariance height mismatch")
        predicted = [
            math.fsum(
                coefficient * parameter
                for coefficient, parameter in zip(row, parameters)
            )
            for row in design
        ]
        source_covariance = prediction_covariance(
            design, parameter_covariance
        )
    residual = [float(observed[i]) - predicted[i] for i in range(2)]
    residual_covariance = matadd(target_covariance, source_covariance)
    chi_square = quadratic_2(residual, residual_covariance)
    return {
        "name": model["name"],
        "basis": model["basis"],
        "source": model.get("source"),
        "parameters": parameters,
        "predicted_P4_S_prime": predicted,
        "source_prediction_covariance": source_covariance,
        "residual": residual,
        "residual_covariance": residual_covariance,
        "chi_square": chi_square,
        "df": 2,
        "signed_z_marginal": [
            residual[i] / math.sqrt(residual_covariance[i][i])
            for i in range(2)
        ],
    }


def score(target: Mapping[str, object], manifest: Mapping[str, object]) -> dict:
    target_sizes, power, models = validate_manifest(manifest)
    observed, target_covariance = validate_target(target, target_sizes)
    results = [
        score_model(
            model,
            target_sizes,
            power,
            observed,
            target_covariance,
        )
        for model in models
    ]
    return {
        "status": "prospective frozen-model score; no target refit",
        "sizes": list(target_sizes),
        "observed_P4_S_prime": observed,
        "target_covariance_P4_S_prime": target_covariance,
        "scoring_order": [result["name"] for result in results],
        "results": results,
        "interpretation_guard": (
            "Report models in frozen chronological order. Do not change source "
            "amplitudes, correction coefficients, target sample count, or model "
            "order after viewing target outcomes."
        ),
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root
        / "predictions/p48_sprime_scoring_manifest_20260828.yaml",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = score(read_json(args.target), read_yaml(args.manifest))
    rendered = json.dumps(result, indent=2, sort_keys=False) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
