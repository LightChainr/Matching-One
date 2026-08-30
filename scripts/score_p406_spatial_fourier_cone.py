#!/usr/bin/env python3
"""Known-support positive Fourier-cone score for the P250 spatial autocorrelation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy.optimize import nnls


SCHEMA = "matching-one/p406-spatial-fourier-cone-score/v1"
GROUP_ORDER = 101
HANDS = ("plus", "minus")
CHARGES = (1, 2)
FIELD = re.compile(r"^a([pm])(\d+)_b([pm])(\d+)_r([12])_(plus|minus)_(re|im)$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def decode(name: str) -> tuple[int, int, int, str, str] | None:
    match = FIELD.match(name)
    if not match:
        return None
    a = int(match.group(2)) * (1 if match.group(1) == "p" else -1)
    b = int(match.group(4)) * (1 if match.group(3) == "p" else -1)
    return a, b, int(match.group(5)), match.group(6), match.group(7)


def read_block(path: Path, expected_hash: str) -> dict[tuple[str, int], dict[str, object]]:
    actual = sha256(path)
    if actual != expected_hash:
        raise ValueError(f"input hash changed: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields: dict[tuple[str, int], list[tuple[tuple[int, int], str, str]]] = {
            (hand, charge): [] for hand in HANDS for charge in CHARGES
        }
        for name in reader.fieldnames or []:
            item = decode(name)
            if item:
                a, b, charge, hand, part = item
                fields[(hand, charge)].append(((a, b), part, name))
        rows = list(reader)
    output = {}
    for key, columns in fields.items():
        if not columns:
            continue
        coordinate_order = []
        for coordinate, _, _ in columns:
            if coordinate not in coordinate_order:
                coordinate_order.append(coordinate)
        expected_columns = [
            next(name for coord, part, name in columns if coord == coordinate and part == component)
            for coordinate in coordinate_order for component in ("re", "im")
        ]
        values = np.asarray([[float(row[name]) for name in expected_columns] for row in rows])
        output[key] = {
            "coordinates": coordinate_order,
            "values": values,
            "batches": len(rows),
            "source_sha256": actual,
        }
    return output


def design(coordinates: list[tuple[int, int]]) -> np.ndarray:
    rows = []
    frequency = np.arange(GROUP_ORDER, dtype=float)
    for a, b in coordinates:
        residue = (a - 10 * b) % GROUP_ORDER
        phase = 2.0 * math.pi * frequency * residue / GROUP_ORDER
        rows.append(np.cos(phase))
        rows.append(-np.sin(phase))
    return np.asarray(rows)


def whiten(values: np.ndarray, model: np.ndarray, rtol: float = 1e-10) -> dict[str, object]:
    batches = len(values)
    mean = values.mean(axis=0)
    centered = values - mean
    covariance = centered.T @ centered / (batches * (batches - 1))
    eigenvalues, eigenvectors = np.linalg.eigh((covariance + covariance.T) / 2.0)
    cutoff = max(float(eigenvalues[-1]) * rtol, 0.0)
    keep = eigenvalues > cutoff
    if not np.any(keep):
        raise ValueError("block covariance has no resolved mode")
    transform = eigenvectors[:, keep].T / np.sqrt(eigenvalues[keep])[:, None]
    return {
        "X": transform @ model,
        "y": transform @ mean,
        "resolved_modes": int(np.sum(keep)),
        "covariance_eigen_cutoff": cutoff,
        "mean": mean,
        "covariance": covariance,
    }


def fit_nonnegative(matrix: np.ndarray, vector: np.ndarray) -> tuple[np.ndarray, float]:
    weights, norm = nnls(matrix, vector, maxiter=100 * matrix.shape[1])
    return weights, float(norm * norm)


def score_channel(blocks: list[dict[str, object]], bootstrap: int, seed: int) -> dict[str, object]:
    whitened = []
    aliases: dict[int, list[list[int]]] = {}
    for block in blocks:
        coordinates = block["coordinates"]
        for coordinate in coordinates:
            aliases.setdefault((coordinate[0] - 10 * coordinate[1]) % GROUP_ORDER, []).append(list(coordinate))
        whitened.append(whiten(block["values"], design(coordinates)))
    matrix = np.vstack([row["X"] for row in whitened])
    vector = np.concatenate([row["y"] for row in whitened])
    weights, statistic = fit_nonnegative(matrix, vector)
    unconstrained, _, unconstrained_rank, _ = np.linalg.lstsq(matrix, vector, rcond=1e-10)
    unconstrained_statistic = float(np.sum((matrix @ unconstrained - vector) ** 2))

    fitted = matrix @ weights
    rng = np.random.default_rng(seed)
    reference = []
    for _ in range(bootstrap):
        _, value = fit_nonnegative(matrix, fitted + rng.standard_normal(len(fitted)))
        reference.append(value)
    p = (1 + sum(value >= statistic for value in reference)) / (bootstrap + 1)
    positive = weights[weights > max(weights.max() * 1e-10, 1e-14)]
    mass = float(weights.sum())
    return {
        "minimum_cone_distance_squared": statistic,
        "bootstrap_p": p,
        "bootstrap_replicates": bootstrap,
        "bootstrap_quantiles": {
            "q50": float(np.quantile(reference, 0.50)),
            "q90": float(np.quantile(reference, 0.90)),
            "q99": float(np.quantile(reference, 0.99)),
        },
        "resolved_whitened_coordinates": len(vector),
        "unconstrained_design_rank": int(unconstrained_rank),
        "unconstrained_distance_squared": unconstrained_statistic,
        "cone_increment": statistic - unconstrained_statistic,
        "positive_weight_count_descriptive": len(positive),
        "total_spectral_mass": mass,
        "inverse_participation_effective_modes": float(mass * mass / np.sum(weights * weights)) if mass else 0.0,
        "weights": weights.tolist(),
        "blocks": [
            {
                "coordinates": [list(value) for value in block["coordinates"]],
                "resolved_modes": white["resolved_modes"],
                "covariance_eigen_cutoff": white["covariance_eigen_cutoff"],
            }
            for block, white in zip(blocks, whitened)
        ],
        "alias_classes": [
            {"residue": residue, "coordinates": coordinates}
            for residue, coordinates in sorted(aliases.items()) if len(coordinates) > 1
        ],
        "decision": "cone_not_rejected" if p >= 0.01 else "cone_rejected",
        "nonidentification": "fitted weights are one NNLS completion and are not a unique physical spectrum",
    }


def score(paths: list[tuple[Path, str]], bootstrap: int, seed: int) -> dict[str, object]:
    archives = [read_block(path, digest) for path, digest in paths]
    channels = {}
    for hand_index, hand in enumerate(HANDS):
        for charge in CHARGES:
            key = hand, charge
            blocks = [archive[key] for archive in archives if key in archive]
            channels[f"{hand}_r{charge}"] = score_channel(
                blocks, bootstrap, seed + 1000 * hand_index + charge
            )
    return {
        "schema": SCHEMA,
        "group_order": GROUP_ORDER,
        "coordinate_map": "j(a,b)=a-10b mod 101",
        "channels": channels,
        "decision": (
            "positive_fourier_cone_not_rejected_in_all_channels"
            if all(row["decision"] == "cone_not_rejected" for row in channels.values())
            else "at_least_one_channel_rejects_positive_fourier_cone"
        ),
        "boundary": (
            "This tests the finite spatial autocorrelation cone. It does not identify physical transfer states, "
            "CFT fields, ordered memory, or a unique spectrum."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius4", type=Path, required=True)
    parser.add_argument("--radius4-sha256", required=True)
    parser.add_argument("--radius5", type=Path, required=True)
    parser.add_argument("--radius5-sha256", required=True)
    parser.add_argument("--radius6", type=Path, required=True)
    parser.add_argument("--radius6-sha256", required=True)
    parser.add_argument("--bootstrap", type=int, default=250)
    parser.add_argument("--seed", type=int, default=40610120260830)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = score([
        (args.radius4, args.radius4_sha256),
        (args.radius5, args.radius5_sha256),
        (args.radius6, args.radius6_sha256),
    ], args.bootstrap, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
