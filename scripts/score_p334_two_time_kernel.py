#!/usr/bin/env python3
"""Score the frozen P334 two-time homology-rank kernel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.stats import chi2


SCHEMA = "matching-one/p334-two-time-kernel-score/v1"
ORIENTATIONS = ("first", "second")


def load(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "matching-one/p334-two-time-kernel-summary/v1":
        raise ValueError("summary schema changed")
    return payload


def moments(rows: Sequence[dict[str, object]]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    samples = sum(int(row["samples"]) for row in rows)
    mean = sum(np.asarray(row["sum_r"], dtype=float) for row in rows) / samples
    second = sum(np.asarray(row["sum_rr"], dtype=float) for row in rows) / samples
    f1 = sum(np.asarray(row["sum_f1"], dtype=float) for row in rows) / samples
    f2 = sum(np.asarray(row["sum_f2"], dtype=float) for row in rows) / samples
    joint = sum(np.asarray(row["sum_joint"], dtype=float) for row in rows) / samples
    return mean, second, f1, f2, joint


def connected(rows: Sequence[dict[str, object]]) -> np.ndarray:
    mean, second, _, _, _ = moments(rows)
    return second - np.outer(mean, mean)


def upper(matrix: np.ndarray) -> np.ndarray:
    return np.asarray([matrix[i, j] for i in range(len(matrix)) for j in range(i, len(matrix))])


def adjacent_defects(matrix: np.ndarray) -> np.ndarray:
    output = []
    for i in range(len(matrix) - 1):
        denominator = matrix[i, i] * matrix[i + 1, i + 1]
        if denominator <= 0:
            raise ValueError("nonpositive temporal variance")
        output.append(1.0 - matrix[i, i + 1] ** 2 / denominator)
    return np.asarray(output)


def jackknife(rows_by_orientation: dict[str, list[dict[str, object]]], feature) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    batches = len(rows_by_orientation[ORIENTATIONS[0]])
    full = np.concatenate([feature(rows_by_orientation[name]) for name in ORIENTATIONS])
    replicates = []
    for omitted in range(batches):
        pieces = []
        for name in ORIENTATIONS:
            pieces.append(feature([row for row in rows_by_orientation[name] if int(row["batch"]) != omitted]))
        replicates.append(np.concatenate(pieces))
    values = np.asarray(replicates)
    centered = values - values.mean(axis=0)
    covariance = (batches - 1) / batches * centered.T @ centered
    return full, covariance, replicates


def wald(vector: np.ndarray, covariance: np.ndarray, *, rtol: float = 1e-10) -> dict[str, object]:
    values, vectors = np.linalg.eigh((covariance + covariance.T) / 2)
    cutoff = max(values[-1] * rtol, 0.0)
    keep = values > cutoff
    if not np.any(keep):
        raise ValueError("contrast covariance has no resolved modes")
    projected = vectors[:, keep].T @ vector
    statistic = float(np.sum(projected * projected / values[keep]))
    df = int(np.sum(keep))
    return {"chi_square": statistic, "df": df, "p": float(chi2.sf(statistic, df)), "eigen_cutoff": cutoff}


def exact_identity(rows: Sequence[dict[str, object]]) -> float:
    _, second, f1, f2, joint = moments(rows)
    residual = 0.0
    for i in range(len(f1)):
        for j in range(i, len(f1)):
            residual = max(residual, abs(joint[i, j] - (second[i, j] - f1[i] - 2.0 * f2[i])))
    return residual


def score_summary(payload: dict[str, object]) -> dict[str, object]:
    grouped = {
        name: [row for row in payload["batches"] if row["orientation"] == name]
        for name in ORIENTATIONS
    }
    defects, defect_covariance, _ = jackknife(grouped, lambda rows: adjacent_defects(connected(rows)))
    kernels, kernel_covariance, kernel_replicates = jackknife(grouped, lambda rows: upper(connected(rows)))
    primary = wald(defects, defect_covariance)
    by_orientation = {}
    for index, name in enumerate(ORIENTATIONS):
        kernel = connected(grouped[name])
        eigenvalues = np.linalg.eigvalsh(kernel)[::-1]
        by_orientation[name] = {
            "adjacent_rank1_defects": defects[index * 6:(index + 1) * 6].tolist(),
            "kernel_eigenvalues": eigenvalues.tolist(),
            "leading_eigen_fraction": float(eigenvalues[0] / np.sum(eigenvalues)),
            "identity_max_abs_residual": exact_identity(grouped[name]),
        }
    return {
        "N": int(payload["N"]),
        "k0": int(payload["k0"]),
        "z_grid": payload["z_grid"],
        "layers": payload["layers"],
        "rank1_separability": primary,
        "joint_defect_vector": defects.tolist(),
        "joint_defect_covariance": defect_covariance.tolist(),
        "kernel_vector": kernels.tolist(),
        "kernel_covariance": kernel_covariance.tolist(),
        "kernel_jackknife": np.asarray(kernel_replicates).tolist(),
        "orientations": by_orientation,
    }


def orientation_mean_kernel(result: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    vector = np.asarray(result["kernel_vector"], dtype=float)
    dimension = len(vector) // 2
    transform = np.concatenate([0.5 * np.eye(dimension), 0.5 * np.eye(dimension)], axis=1)
    covariance = np.asarray(result["kernel_covariance"], dtype=float)
    return transform @ vector, transform @ covariance @ transform.T


def score(first: dict[str, object], second: dict[str, object]) -> dict[str, object]:
    rows = [score_summary(first), score_summary(second)]
    mean0, cov0 = orientation_mean_kernel(rows[0])
    mean1, cov1 = orientation_mean_kernel(rows[1])
    transfer = wald(mean1 - mean0, cov0 + cov1)
    return {
        "schema": SCHEMA,
        "decision": "two_time_rank1_separability_rejected" if all(row["rank1_separability"]["p"] < 0.01 for row in rows) else "mixed_rank1_result",
        "sizes": rows,
        "scale_transfer_N425_minus_N325": {
            **transfer,
            "contrast": (mean1 - mean0).tolist(),
            "boundary": "diagnostic equality test across different finite quotients, not a universal-collapse theorem"
        },
        "scientific_card": {
            "mechanism_changed": "one separable scalar temporal amplitude is removed for the production homology-rank process",
            "not_proved": "CFT state count, Jordan structure, intrinsic memory after conditioning on full geometry, or universal kernel collapse",
            "observer_sector_source_geometry": "same-permutation two-time ambient-H1 rank process from N325/N425 paired norm-five quotient orientations",
            "dependency_groups": "one paired-orientation block per size; the two sizes are independent",
            "next_lift": "condition the kernel on the newly recorded current-k0 geometry and test which nonleading temporal mode remains"
        }
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("n325", type=Path)
    parser.add_argument("n425", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = score(load(args.n325), load(args.n425))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
