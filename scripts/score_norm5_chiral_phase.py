#!/usr/bin/env python3
"""Covariance-aware scorer for the frozen N325 chiral Hecke phase test."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import atan2, exp, pi, sqrt
from pathlib import Path
from typing import Sequence


PRIMARY_ORDER = ("plus_re", "plus_im", "minus_re", "minus_im")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inverse(matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("matrix must be nonempty and square")
    augmented = [
        [float(value) for value in row]
        + [1.0 if column == row_index else 0.0 for column in range(n)]
        for row_index, row in enumerate(matrix)
    ]
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1e-15:
            raise ValueError("matrix is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(augmented[row], augmented[column])
            ]
    return [row[n:] for row in augmented]


def matvec(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    return [sum(value * vector[j] for j, value in enumerate(row)) for row in matrix]


def quadratic(vector: Sequence[float], matrix: Sequence[Sequence[float]]) -> float:
    product = matvec(matrix, vector)
    return sum(value * product[index] for index, value in enumerate(vector))


def gls_model(
    mean: Sequence[float], covariance: Sequence[Sequence[float]], q: complex
) -> dict:
    precision = inverse(covariance)
    design = [
        [q.real, -q.imag],
        [q.imag, q.real],
        [1.0, 0.0],
        [0.0, 1.0],
    ]
    # Normal equations X^T C^-1 X beta = X^T C^-1 y.
    precision_design = [
        [sum(precision[i][k] * design[k][j] for k in range(4)) for j in range(2)]
        for i in range(4)
    ]
    normal = [
        [sum(design[k][i] * precision_design[k][j] for k in range(4)) for j in range(2)]
        for i in range(2)
    ]
    precision_mean = matvec(precision, mean)
    rhs = [sum(design[k][i] * precision_mean[k] for k in range(4)) for i in range(2)]
    normal_inverse = inverse(normal)
    beta = matvec(normal_inverse, rhs)
    fitted = matvec(design, beta)
    residual = [mean[index] - fitted[index] for index in range(4)]
    chi_square = quadratic(residual, precision)
    return {
        "fitted_common_complex_normalization_re_im": beta,
        "normalization_covariance": normal_inverse,
        "fitted_primary_vector": fitted,
        "residual_primary_vector": residual,
        "chi_square": chi_square,
        "degrees_of_freedom": 2,
        "chi_square_per_df": chi_square / 2.0,
    }


def wrapped_degrees(value: float) -> float:
    return (value + 180.0) % 360.0 - 180.0


def score_reflection_null(analysis: dict) -> dict:
    null = analysis["true_reflection_conjugacy_null"]
    point = [float(value) for value in null["point_re_im"]]
    covariance = [[float(value) for value in row] for row in null["covariance_of_mean"]]
    if point == [0.0, 0.0] and covariance == [[0.0, 0.0], [0.0, 0.0]]:
        return {
            "status": "exact_configurationwise_null",
            "point_re_im": point,
            "chi_square": 0.0,
            "degrees_of_freedom": 0,
            "relation": null["relation"],
        }
    try:
        chi_square = quadratic(point, inverse(covariance))
    except ValueError:
        return {
            "status": "singular_covariance_nonexact_null",
            "point_re_im": point,
            "covariance": covariance,
            "relation": null["relation"],
        }
    return {
        "status": "statistical_null",
        "point_re_im": point,
        "covariance": covariance,
        "chi_square": chi_square,
        "degrees_of_freedom": 2,
        "chi_square_per_df": chi_square / 2.0,
        "relation": null["relation"],
    }


def validate_inputs(payload: dict, manifest: dict, contract: dict) -> None:
    if payload.get("schema") != "matching-one.norm5-chiral-fixedp-response.v1":
        raise ValueError("unexpected response schema")
    if tuple(payload["analysis"]["primary_order"]) != PRIMARY_ORDER:
        raise ValueError("response primary order changed")
    if tuple(contract["primary_order"]) != PRIMARY_ORDER:
        raise ValueError("contract primary order changed")
    frozen = manifest["run"]
    observed = payload["run"]
    for key in ("samples", "batches", "workers", "p", "seed", "radius"):
        if observed[key] != frozen[key]:
            raise ValueError(f"response run field {key} differs from manifest")
    if payload["mapping_gate"].get("passed") is not True:
        raise ValueError("cover mapping gate did not pass")
    if payload["analysis"].get("same_parent_plus_minus_conjugacy_is_not_a_null") is not True:
        raise ValueError("same-parent conjugacy semantics changed")


def score_payload(payload: dict, manifest: dict, contract: dict) -> dict:
    validate_inputs(payload, manifest, contract)
    analysis = payload["analysis"]
    mean = [float(value) for value in analysis["primary_point"]]
    covariance = [
        [float(value) for value in row]
        for row in analysis["primary_covariance_of_mean"]
    ]
    models = {}
    observed_minus = complex(mean[2], mean[3])
    observed_ratio = complex(mean[0], mean[1]) / observed_minus if observed_minus else None
    for name, target in contract["targets"].items():
        q = complex(float(Fraction(target["q_re"])), float(Fraction(target["q_im"])))
        result = gls_model(mean, covariance, q)
        target_phase = atan2(q.imag, q.real) * 180.0 / pi
        observed_phase = (
            atan2(observed_ratio.imag, observed_ratio.real) * 180.0 / pi
            if observed_ratio is not None
            else None
        )
        result.update(
            {
                "spin": target["spin"],
                "exact_ratio": {"real": target["q_re"], "imag": target["q_im"]},
                "target_phase_degrees": target_phase,
                "observed_ratio_phase_degrees": observed_phase,
                "wrapped_phase_residual_degrees": (
                    wrapped_degrees(observed_phase - target_phase)
                    if observed_phase is not None
                    else None
                ),
            }
        )
        models[name] = result
    minimum = min(model["chi_square"] for model in models.values())
    relative = {
        name: exp(-0.5 * (model["chi_square"] - minimum))
        for name, model in models.items()
    }
    normalization = sum(relative.values())
    for name, model in models.items():
        model["delta_chi_square"] = model["chi_square"] - minimum
        model["relative_likelihood_to_best"] = relative[name]
        model["normalized_weight_over_frozen_three"] = relative[name] / normalization
    ranking = sorted(models, key=lambda name: models[name]["chi_square"])
    return {
        "schema": "matching-one.norm5-chiral-phase-score.v1",
        "status": "scored_under_pre_reveal_contract",
        "primary_vector": mean,
        "primary_covariance": covariance,
        "joint_model": "R_minus=c, R_plus=q_s*c; one free complex normalization",
        "models": models,
        "ranking": ranking,
        "best_model": ranking[0],
        "minimum_phase_residual_model": min(
            models,
            key=lambda name: abs(models[name]["wrapped_phase_residual_degrees"])
            if models[name]["wrapped_phase_residual_degrees"] is not None
            else float("inf"),
        ),
        "reflection_conjugacy_null": score_reflection_null(analysis),
        "evidence_accounting": "one joint 4D comparison; Re/Im and the two hands are not independent votes",
    }


def synthetic_payload(manifest: dict, contract: dict, target_name: str = "H8") -> dict:
    """Correlated exact-model oracle; off-diagonals make diagonal scoring wrong."""

    target = contract["targets"][target_name]
    q = complex(float(Fraction(target["q_re"])), float(Fraction(target["q_im"])))
    common = complex(2.0, -1.0)
    plus = q * common
    covariance = [
        [0.040, 0.012, 0.009, -0.006],
        [0.012, 0.030, 0.004, 0.008],
        [0.009, 0.004, 0.025, -0.007],
        [-0.006, 0.008, -0.007, 0.035],
    ]
    inverse(covariance)  # positive nonsingular gate for the fixture
    return {
        "schema": "matching-one.norm5-chiral-fixedp-response.v1",
        "status": "synthetic_oracle",
        "mapping_gate": {"passed": True},
        "run": {
            key: manifest["run"][key]
            for key in ("samples", "batches", "workers", "p", "seed", "radius")
        },
        "analysis": {
            "primary_order": list(PRIMARY_ORDER),
            "primary_point": [plus.real, plus.imag, common.real, common.imag],
            "primary_covariance_of_mean": covariance,
            "same_parent_plus_minus_conjugacy_is_not_a_null": True,
            "true_reflection_conjugacy_null": {
                "point_re_im": [0.0, 0.0],
                "covariance_of_mean": [[0.0, 0.0], [0.0, 0.0]],
                "relation": "synthetic exact reflection transport",
            },
        },
        "synthetic_truth": {
            "target": target_name,
            "common_complex_normalization_re_im": [common.real, common.imag],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input", type=Path)
    source.add_argument("--synthetic-target", choices=("H4", "H8", "H12"))
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text())
    manifest_hash = sha256(args.manifest)
    if manifest_hash != contract["production_manifest_sha256"]:
        raise ValueError("production manifest hash differs from frozen scorer contract")
    manifest = json.loads(args.manifest.read_text())
    payload = (
        json.loads(args.input.read_text())
        if args.input is not None
        else synthetic_payload(manifest, contract, args.synthetic_target)
    )
    result = score_payload(payload, manifest, contract)
    result["provenance"] = {
        "input": str(args.input) if args.input is not None else f"synthetic:{args.synthetic_target}",
        "input_sha256": sha256(args.input) if args.input is not None else None,
        "manifest": str(args.manifest),
        "manifest_sha256": manifest_hash,
        "contract": str(args.contract),
        "contract_sha256": sha256(args.contract),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
