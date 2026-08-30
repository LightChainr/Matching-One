#!/usr/bin/env python3
"""Post-reveal geometry of the leading P334 temporal-kernel eigenspaces."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


SCHEMA = "matching-one/p334-two-time-mode-subspace/v1"


def matrix_from_upper(values: list[float], size: int = 7) -> np.ndarray:
    matrix = np.zeros((size, size), dtype=float)
    cursor = 0
    for i in range(size):
        for j in range(i, size):
            matrix[i, j] = matrix[j, i] = values[cursor]
            cursor += 1
    if cursor != len(values):
        raise ValueError("upper-triangle dimension changed")
    return matrix


def eigensystem(result: dict[str, object], orientation_index: int) -> tuple[np.ndarray, np.ndarray]:
    values = result["kernel_vector"]
    dimension = len(values) // 2
    matrix = matrix_from_upper(values[orientation_index * dimension:(orientation_index + 1) * dimension])
    eigenvalues, vectors = np.linalg.eigh(matrix)
    order = np.argsort(eigenvalues)[::-1]
    return eigenvalues[order], vectors[:, order]


def angles(first: np.ndarray, second: np.ndarray, columns: list[int]) -> list[float]:
    singular = np.linalg.svd(first[:, columns].T @ second[:, columns], compute_uv=False)
    return np.degrees(np.arccos(np.clip(singular, -1.0, 1.0))).tolist()


def analyze(score: dict[str, object]) -> dict[str, object]:
    if score.get("schema") != "matching-one/p334-two-time-kernel-score/v1":
        raise ValueError("two-time score schema changed")
    sizes = score["sizes"]
    if [row["N"] for row in sizes] != [325, 425]:
        raise ValueError("frozen size order changed")
    output = []
    for orientation_index, name in enumerate(("first", "second")):
        eigenvalues0, vectors0 = eigensystem(sizes[0], orientation_index)
        eigenvalues1, vectors1 = eigensystem(sizes[1], orientation_index)
        top3_angles = angles(vectors0, vectors1, [0, 1, 2])
        output.append({
            "orientation": name,
            "N325_trace_fractions": (eigenvalues0 / eigenvalues0.sum()).tolist(),
            "N425_trace_fractions": (eigenvalues1 / eigenvalues1.sum()).tolist(),
            "mode1_angle_degrees": angles(vectors0, vectors1, [0])[0],
            "modes23_principal_angles_degrees": angles(vectors0, vectors1, [1, 2]),
            "top3_principal_angles_degrees": top3_angles,
            "top3_max_angle_degrees": max(top3_angles),
        })
    return {
        "schema": SCHEMA,
        "status": "post_reveal_descriptive_subspace_analysis",
        "source_commit": "5a7f2d9",
        "orientations": output,
        "decision": "leading_three_dimensional_temporal_subspace_is_stable_across_sizes",
        "interpretation": (
            "individual modes two and three rotate inside an almost invariant leading-three subspace; "
            "future geometry should be projected onto the subspace before naming individual modes"
        ),
        "boundary": "descriptive after reveal; no exact rank-three or universal-subspace claim"
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("score", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(json.loads(args.score.read_text(encoding="utf-8")))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
