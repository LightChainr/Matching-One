#!/usr/bin/env python3
"""Frozen covariance-aware scorer for the Issue #200 N650 mixed join."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp


PRIMARY = ("ES", "ED", "OS", "OD")
SECONDARY = ("ambient_ES", "ambient_ED", "ambient_OS", "ambient_OD")
PREDICTION_SHA256 = "8c043677605f9c8d0a00ee11fcc93539da9ea87bbacb93553e76ea26ec27ba5f"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def load_inputs(batch_path: Path, metadata_path: Path, prediction_path: Path):
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("schema") != "matching-one.p200-n650-mixed-join-run.v1":
        raise ValueError("wrong N650 mixed-join metadata schema")
    if metadata.get("p_ref") != "0.592746050790" or metadata.get("stored_sum_divisor") != 2:
        raise ValueError("p_ref or stored-sum convention changed")
    if metadata.get("state_order") != list(PRIMARY):
        raise ValueError("primary state order changed")
    prediction_hash = file_sha256(prediction_path)
    if prediction_hash != PREDICTION_SHA256:
        raise ValueError("Phase B prediction artifact hash changed")
    prediction = json.loads(prediction_path.read_text(encoding="utf-8"))
    if prediction.get("schema") != "matching-one.p200-n650-mixed-join-phaseB.v1":
        raise ValueError("wrong Phase B prediction schema")

    rows = []
    with batch_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            samples = int(row["samples"])
            primary = [int(row[f"{name}_num_sum"]) for name in PRIMARY]
            secondary = [int(row[f"{name}_num_sum"]) for name in SECONDARY]
            rows.append({"batch": int(row["batch"]), "samples": samples, "primary": primary, "secondary": secondary})
    if len(rows) != metadata.get("batches") or sum(row["samples"] for row in rows) != metadata.get("samples"):
        raise ValueError("batch count/sample total disagrees with metadata")
    if len(rows) < 5:
        raise ValueError("at least five batches are required for a four-state covariance")
    return metadata, prediction_hash, rows


def delete_one_summary(rows: list[dict], key: str) -> tuple[list[mp.mpf], list[list[mp.mpf]], list[list[mp.mpf]]]:
    dimension = 4
    total_samples = sum(row["samples"] for row in rows)
    total = [sum(row[key][index] for row in rows) for index in range(dimension)]
    mean = [mp.mpf(total[index]) / (2 * total_samples) for index in range(dimension)]
    leave_one = []
    for row in rows:
        denominator = 2 * (total_samples - row["samples"])
        leave_one.append([mp.mpf(total[index] - row[key][index]) / denominator for index in range(dimension)])
    center = [mp.fsum(values[index] for values in leave_one) / len(leave_one) for index in range(dimension)]
    factor = mp.mpf(len(rows) - 1) / len(rows)
    covariance = [
        [
            factor * mp.fsum((value[i] - center[i]) * (value[j] - center[j]) for value in leave_one)
            for j in range(dimension)
        ]
        for i in range(dimension)
    ]
    return mean, covariance, leave_one


def gls_zero(mean: list[mp.mpf], covariance: list[list[mp.mpf]]) -> dict:
    matrix = mp.matrix(covariance)
    eigenvalues, eigenvectors = mp.eigsy(matrix)
    values = [mp.mpf(eigenvalues[index]) for index in range(4)]
    maximum = max(values)
    tolerance = max(mp.mpf("1e-40"), maximum * mp.mpf("1e-10"))
    vector = mp.matrix(mean)
    projections = eigenvectors.T * vector
    live = [index for index, value in enumerate(values) if value > tolerance]
    chi_square = mp.fsum(projections[index] ** 2 / values[index] for index in live)
    degrees = len(live)
    p_value = mp.gammainc(mp.mpf(degrees) / 2, chi_square / 2, mp.inf, regularized=True) if degrees else mp.nan
    marginal = []
    for index, name in enumerate(PRIMARY):
        variance = covariance[index][index]
        marginal.append({"name": name, "mean": float(mean[index]), "se": float(mp.sqrt(variance)), "z": float(mean[index] / mp.sqrt(variance)) if variance > 0 else None})
    return {
        "chi_square": float(chi_square),
        "degrees_of_freedom": degrees,
        "p_value": float(p_value),
        "covariance_eigenvalues": [float(value) for value in values],
        "pseudoinverse_tolerance": float(tolerance),
        "marginal_diagnostics": marginal,
    }


def matrix_float(matrix):
    return [[float(value) for value in row] for row in matrix]


def render(batch_path: Path, metadata_path: Path, prediction_path: Path) -> dict:
    mp.mp.dps = 50
    metadata, prediction_hash, rows = load_inputs(batch_path, metadata_path, prediction_path)
    primary_mean, primary_covariance, primary_leave_one = delete_one_summary(rows, "primary")
    secondary_mean, secondary_covariance, _ = delete_one_summary(rows, "secondary")
    score = gls_zero(primary_mean, primary_covariance)
    leading = max(
        range(4),
        key=lambda index: abs(score["marginal_diagnostics"][index]["z"] or 0.0),
    )
    leading_sign = mp.sign(primary_mean[leading])
    sign_stability = sum(mp.sign(value[leading]) == leading_sign for value in primary_leave_one) / len(primary_leave_one)
    stop_recommended = score["p_value"] < 1e-6 and sign_stability >= 0.8
    return {
        "schema": "matching-one.p200-n650-mixed-join-score.v1",
        "issue": 200,
        "status": "frozen_interface_score",
        "input": {
            "batches": str(batch_path),
            "metadata": str(metadata_path),
            "prediction": str(prediction_path),
            "prediction_sha256": prediction_hash,
            "run_git_commit": metadata["git_commit"],
            "samples": metadata["samples"],
            "batch_count": metadata["batches"],
        },
        "primary": {
            "state_order": list(PRIMARY),
            "mean": [float(value) for value in primary_mean],
            "delete_one_covariance": matrix_float(primary_covariance),
            "null": [0.0] * 4,
            "joint_GLS": score,
        },
        "secondary_ambient_H1": {
            "state_order": list(SECONDARY),
            "mean": [float(value) for value in secondary_mean],
            "delete_one_covariance": matrix_float(secondary_covariance),
            "role": "correlated mechanism diagnostic; not an additional primary vote",
        },
        "decision": {
            "stop_if_p_below": 1e-6,
            "stop_if": "the factor-additive zero bridge is already overwhelmingly rejected and the leading channel sign is stable across at least 80% of leave-one batches",
            "continue_if": "joint p is not decisive or the leading S/D/color channel is unstable",
            "interpretation_boundary": "rejection identifies connected C2-by-C5 interaction, not path noncommutativity or Jordan by itself",
            "leading_channel": PRIMARY[leading],
            "leading_leave_one_sign_stability": sign_stability,
            "stop_recommended": stop_recommended,
        },
    }


def report(payload: dict) -> str:
    primary = payload["primary"]
    score = primary["joint_GLS"]
    lines = [
        "# P200 N650 mixed-join score",
        "",
        f"Samples: `{payload['input']['samples']}` in `{payload['input']['batch_count']}` batches.",
        "",
        f"Joint zero GLS: chi2=`{score['chi_square']:.8g}`, df=`{score['degrees_of_freedom']}`, p=`{score['p_value']:.8g}`.",
        "",
        "| state | mean | SE | z |",
        "|---|---:|---:|---:|",
    ]
    for row in score["marginal_diagnostics"]:
        lines.append(f"| {row['name']} | {row['mean']:.9g} | {row['se']:.6g} | {row['z']:.6g} |")
    lines += ["", "A rejection is a connected mixed-factor interaction, not evidence of chronological path memory by itself."]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--prediction", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    payload = render(args.batches, args.metadata, args.prediction)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.report.write_text(report(payload), encoding="utf-8")


if __name__ == "__main__":
    main()
