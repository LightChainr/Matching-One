#!/usr/bin/env python3
"""Score the frozen N505 adaptive x spatial-spectrum pilot.

The positive Fourier cone is an ordinary spatial-state adversary, not an exact
constraint on the adaptive intervention observable.  The sharper frozen test is
whether its spatial response is exhausted by the independently measured
defined-rate profile.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import nnls
from scipy.stats import chi2


ORDER = 101
HANDS = ("plus", "minus")
RESIDUES = tuple(range(1, ORDER))
SYMMETRIC_RESIDUES = tuple(range(1, (ORDER + 1) // 2))
ALPHA = 1.7873118647707438
ALPHA_SE = 0.010073378702833143
SCHEMA = "matching-one/p250-adaptive-spatial-joint-score/v1"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_profiles(path: Path, expected_sha256: str | None):
    actual = sha256(path)
    if expected_sha256 and actual != expected_sha256:
        raise ValueError(f"batch archive hash changed: {actual} != {expected_sha256}")
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) < 10:
        raise ValueError("at least ten batches are required")
    output = {}
    for hand in HANDS:
        samples = np.asarray([
            [float(row[f"{hand}_j{j}_samples"]) for j in RESIDUES]
            for row in rows
        ])
        if np.any(samples <= 0):
            where = np.argwhere(samples <= 0)[0]
            raise ValueError(
                f"empty spatial cell at batch={where[0]}, residue={where[1] + 1}; "
                "increase samples per batch"
            )
        output[hand] = {
            "response": np.asarray([
                [float(row[f"{hand}_j{j}_sum_Rminus"]) for j in RESIDUES]
                for row in rows
            ]) / samples,
            "defined": np.asarray([
                [float(row[f"{hand}_j{j}_defined"]) for j in RESIDUES]
                for row in rows
            ]) / samples,
            "tie": np.asarray([
                [float(row[f"{hand}_j{j}_ties"]) for j in RESIDUES]
                for row in rows
            ]) / samples,
        }
    return rows, output, actual


def symmetrize(values: np.ndarray) -> np.ndarray:
    return np.column_stack([
        0.5 * (values[:, j - 1] + values[:, ORDER - j - 1])
        for j in SYMMETRIC_RESIDUES
    ])


def covariance_of_mean(values: np.ndarray) -> np.ndarray:
    centered = values - values.mean(axis=0)
    return centered.T @ centered / (len(values) * (len(values) - 1))


def whitening(covariance: np.ndarray, rtol: float = 1e-10):
    covariance = 0.5 * (covariance + covariance.T)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    cutoff = max(float(eigenvalues[-1]) * rtol, 0.0)
    keep = eigenvalues > cutoff
    if not np.any(keep):
        raise ValueError("covariance has no resolved mode")
    transform = eigenvectors[:, keep].T / np.sqrt(eigenvalues[keep])[:, None]
    return transform, eigenvalues, keep, cutoff


def cosine_design() -> np.ndarray:
    frequencies = np.arange((ORDER + 1) // 2, dtype=float)
    return np.asarray([
        np.cos(2.0 * math.pi * frequencies * j / ORDER)
        for j in SYMMETRIC_RESIDUES
    ])


def fit_cone(values_by_hand: dict[str, np.ndarray], bootstrap: int, seed: int):
    values = np.column_stack([values_by_hand[hand] for hand in HANDS])
    mean = values.mean(axis=0)
    covariance = covariance_of_mean(values)
    transform, eigenvalues, keep, cutoff = whitening(covariance)
    base = cosine_design()
    zeros = np.zeros_like(base)
    model = np.block([[base, zeros], [zeros, base]])
    x, y = transform @ model, transform @ mean
    weights, norm = nnls(x, y, maxiter=200 * x.shape[1])
    statistic = float(norm * norm)
    free, _, rank, _ = np.linalg.lstsq(x, y, rcond=1e-10)
    free_statistic = float(np.sum((x @ free - y) ** 2))
    fitted = x @ weights
    rng = np.random.default_rng(seed)
    reference = []
    for _ in range(bootstrap):
        _, boot_norm = nnls(
            x, fitted + rng.standard_normal(len(fitted)), maxiter=200 * x.shape[1]
        )
        reference.append(float(boot_norm * boot_norm))
    p = (1 + sum(value >= statistic for value in reference)) / (bootstrap + 1)
    positive = weights[weights > max(float(weights.max()) * 1e-10, 1e-14)]
    mass = float(weights.sum())
    return {
        "minimum_cone_distance_squared": statistic,
        "unconstrained_distance_squared": free_statistic,
        "cone_increment": statistic - free_statistic,
        "bootstrap_p": p,
        "bootstrap_replicates": bootstrap,
        "bootstrap_quantiles": {
            "q50": float(np.quantile(reference, 0.50)),
            "q90": float(np.quantile(reference, 0.90)),
            "q99": float(np.quantile(reference, 0.99)),
        },
        "resolved_modes": int(np.sum(keep)),
        "covariance_eigen_cutoff": cutoff,
        "unconstrained_design_rank": int(rank),
        "positive_weight_count_descriptive": len(positive),
        "inverse_participation_effective_modes": (
            float(mass * mass / np.sum(weights * weights)) if mass else 0.0
        ),
        "weights": {
            hand: weights[index * 51:(index + 1) * 51].tolist()
            for index, hand in enumerate(HANDS)
        },
        "decision": "ordinary_cone_not_rejected" if p >= 0.01 else "ordinary_cone_rejected",
    }


def definition_only(profiles: dict[str, dict[str, np.ndarray]]):
    residual_batches = np.column_stack([
        symmetrize(profiles[hand]["response"] - ALPHA * profiles[hand]["defined"])
        for hand in HANDS
    ])
    defined_mean = np.concatenate([
        symmetrize(profiles[hand]["defined"]).mean(axis=0) for hand in HANDS
    ])
    mean = residual_batches.mean(axis=0)
    covariance = covariance_of_mean(residual_batches)
    covariance += ALPHA_SE * ALPHA_SE * np.outer(defined_mean, defined_mean)
    transform, eigenvalues, keep, cutoff = whitening(covariance)
    whitened = transform @ mean
    statistic = float(whitened @ whitened)
    degrees = int(np.sum(keep))
    p = float(chi2.sf(statistic, degrees))
    return {
        "source_frozen_alpha": ALPHA,
        "source_alpha_se": ALPHA_SE,
        "chi_square": statistic,
        "degrees_of_freedom_resolved": degrees,
        "asymptotic_p": p,
        "covariance_eigen_cutoff": cutoff,
        "maximum_absolute_spatial_residual": float(np.max(np.abs(mean))),
        "decision": "definition_only_not_rejected" if p >= 0.01 else "definition_only_rejected",
    }


def global_summary(rows):
    result = {}
    for hand in HANDS:
        samples = np.asarray([float(row["samples"]) for row in rows])
        defined = np.asarray([float(row[f"{hand}_defined"]) for row in rows])
        response = np.asarray([float(row[f"{hand}_sum_Rminus"]) for row in rows])
        tie = np.asarray([float(row[f"{hand}_ties"]) for row in rows])
        batch_defined = defined / samples
        batch_unconditional = response / samples
        batch_conditional = response / defined
        result[hand] = {
            "defined_rate": mean_se(batch_defined),
            "tie_rate": mean_se(tie / samples),
            "unconditional_Rminus": mean_se(batch_unconditional),
            "conditional_Rminus": mean_se(batch_conditional),
            "conditional_minus_source_alpha": mean_se(batch_conditional - ALPHA),
        }
    return result


def mean_se(values):
    return {"mean": float(np.mean(values)),
            "se": float(np.std(values, ddof=1) / math.sqrt(len(values)))}


def joint_cosine_archive(profiles):
    base = cosine_design()
    ordering, columns = [], []
    for hand in HANDS:
        for observable in ("response", "defined"):
            transformed = symmetrize(profiles[hand][observable]) @ base
            transformed /= len(SYMMETRIC_RESIDUES)
            for frequency in range(51):
                ordering.append(f"{hand}:{observable}:k{frequency}")
            columns.append(transformed)
    values = np.column_stack(columns)
    return {
        "transform": "cosine coefficients of j<->-j symmetrized nonzero-coordinate profile",
        "zero_coordinate": "not sampled; excluded before transform",
        "ordering": ordering,
        "mean": values.mean(axis=0).tolist(),
        "covariance_of_mean": covariance_of_mean(values).tolist(),
        "resolved_rank_upper_bound": len(values) - 1,
    }


def score(batch_path: Path, expected_hash: str | None, bootstrap: int, seed: int):
    rows, profiles, source_hash = read_profiles(batch_path, expected_hash)
    cone = fit_cone(
        {hand: symmetrize(profiles[hand]["response"]) for hand in HANDS},
        bootstrap, seed,
    )
    definition = definition_only(profiles)
    return {
        "schema": SCHEMA,
        "source": {"batch_path": str(batch_path), "sha256": source_hash,
                   "batches": len(rows)},
        "frozen_source_calibration": {
            "alpha_interpretation": "combined N325/N425 conditional Rminus per defined event",
            "alpha": ALPHA, "se": ALPHA_SE,
        },
        "global": global_summary(rows),
        "ordinary_positive_fourier_cone_adversary": cone,
        "definition_only_model": definition,
        "joint_cosine_periodogram": joint_cosine_archive(profiles),
        "primary_mechanism_decision": (
            "adaptive_spatial_increment" if definition["decision"] == "definition_only_rejected"
            else "definition_rate_exhausts_pilot_resolution"
        ),
        "boundary": (
            "The cone is an ordinary spatial-state adversary, not an exact theorem for an adaptive "
            "intervention. Rejection does not identify a field or memory state. The zero residue was "
            "excluded before reveal to keep the three marked locations distinct."
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--batches-sha256")
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=25050510120263001)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = score(args.batches, args.batches_sha256, args.bootstrap, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
