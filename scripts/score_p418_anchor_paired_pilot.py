#!/usr/bin/env python3
"""Score the single frozen P418 paired-anchor pilot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import chi2

import score_p406_spatial_fourier_cone as p406
import score_p418_crt_degauging as p418


ROOT = Path(__file__).resolve().parents[1]
FREEZE = ROOT / "analysis/p418_anchor_paired_pilot_freeze.json"
OLD4 = ROOT / "results/huawei-20260830/P250-projective-leg-bivariate-80k/response_80k.batches.csv"
DEFAULT_BATCHES = ROOT / "results/local-20260830/P418-anchor-paired-5k/response_5k.batches.csv"
DEFAULT_OUTPUT = ROOT / "results/local-20260830/P418-anchor-paired-5k/score.json"
ESTIMATORS = ("current", "independent", "full")
HANDS = ("plus", "minus")
CHARGES = (1, 2)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_paired(path: Path) -> tuple[list[dict[str, str]], dict[tuple[str, str, int], dict[str, Any]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames or []
    blocks = {}
    for estimator in ESTIMATORS:
        for hand in HANDS:
            for charge in CHARGES:
                columns = []
                for name in fieldnames:
                    prefix = estimator + "__"
                    if not name.startswith(prefix):
                        continue
                    decoded = p406.decode(name[len(prefix):])
                    if decoded and decoded[2] == charge and decoded[3] == hand:
                        columns.append((decoded[:2], decoded[4], name))
                coordinate_order = []
                for coordinate, _, _ in columns:
                    if coordinate not in coordinate_order:
                        coordinate_order.append(coordinate)
                expected = [
                    next(name for coord, part, name in columns if coord == coordinate and part == component)
                    for coordinate in coordinate_order for component in ("re", "im")
                ]
                values = np.asarray([
                    [float(row[name]) / int(row["samples"]) for name in expected]
                    for row in rows
                ])
                blocks[(estimator, hand, charge)] = {
                    "coordinates": coordinate_order,
                    "values": values,
                    "columns": expected,
                }
    return rows, blocks


def fit_cones(exact: dict[str, Any], block: dict[str, Any], hand: str, charge: int, bootstrap: int, seed: int) -> dict[str, Any]:
    coordinates = block["coordinates"]
    raw = p406.whiten(block["values"], p406.design(coordinates))
    masked = p406.whiten(block["values"], p418.masked_design(exact, hand, charge, coordinates))
    output = {}
    noises = np.random.default_rng(seed).standard_normal((bootstrap, len(raw["y"])))
    for name, white in (("raw", raw), ("masked", masked)):
        weights, statistic = p406.fit_nonnegative(white["X"], white["y"])
        fitted = white["X"] @ weights
        reference = [
            p406.fit_nonnegative(white["X"], fitted + noise)[1]
            for noise in noises
        ]
        output[name] = {
            "distance_squared": statistic,
            "bootstrap_p": (1 + sum(value >= statistic for value in reference)) / (bootstrap + 1),
            "resolved_modes": len(white["y"]),
            "bootstrap_q99": float(np.quantile(reference, 0.99)),
        }
    output["masked_minus_raw"] = output["masked"]["distance_squared"] - output["raw"]["distance_squared"]
    return output


def paired_zero(values: np.ndarray, rtol: float = 1e-10) -> dict[str, Any]:
    batches = len(values)
    mean = values.mean(axis=0)
    centered = values - mean
    covariance = centered.T @ centered / (batches * (batches - 1))
    eigenvalues, eigenvectors = np.linalg.eigh((covariance + covariance.T) / 2)
    cutoff = max(float(eigenvalues[-1]) * rtol, 0.0)
    keep = eigenvalues > cutoff
    transform = eigenvectors[:, keep].T / np.sqrt(eigenvalues[keep])[:, None]
    whitened = transform @ mean
    statistic = float(whitened @ whitened)
    degrees = int(np.sum(keep))
    return {
        "mahalanobis_squared": statistic,
        "resolved_modes": degrees,
        "chi2_reference_p": float(chi2.sf(statistic, degrees)),
        "rms_coordinate_difference": float(np.sqrt(np.mean(mean * mean))),
    }


def flattened_correlation(first: np.ndarray, second: np.ndarray) -> float:
    left = first - first.mean(axis=0)
    right = second - second.mean(axis=0)
    denominator = math.sqrt(float(np.sum(left * left) * np.sum(right * right)))
    return float(np.sum(left * right) / denominator) if denominator else 0.0


def historical_replay(rows: list[dict[str, str]], blocks: dict[tuple[str, str, int], dict[str, Any]]) -> dict[str, Any]:
    with OLD4.open(newline="") as handle:
        old_rows = list(csv.DictReader(handle))[:25]
    maximum = 0.0
    comparisons = 0
    for hand in HANDS:
        for charge in CHARGES:
            block = blocks[("current", hand, charge)]
            for old_index, old in enumerate(old_rows):
                for coordinate_index, name in enumerate(block["columns"]):
                    original_name = name.split("__", 1)[1]
                    new_sum = sum(
                        float(rows[4 * old_index + offset][name]) for offset in range(4)
                    )
                    maximum = max(maximum, abs(new_sum - float(old[original_name])))
                    comparisons += 1
    return {
        "old_radius4_first_batches": len(old_rows),
        "new_batches_per_old_batch": 4,
        "coordinates_compared": comparisons,
        "maximum_absolute_sum_error": maximum,
        "passed": maximum < 1e-9,
    }


def build_score(path: Path, bootstrap: int, seed: int) -> dict[str, Any]:
    freeze = json.loads(FREEZE.read_text())
    rows, blocks = read_paired(path)
    exact = p418.exact_section_and_masks()
    channels = {}
    decisions = []
    for hand_index, hand in enumerate(HANDS):
        for charge in CHARGES:
            key = f"{hand}_r{charge}"
            estimator_scores = {}
            for estimator_index, estimator in enumerate(ESTIMATORS):
                estimator_scores[estimator] = fit_cones(
                    exact,
                    blocks[(estimator, hand, charge)],
                    hand,
                    charge,
                    bootstrap,
                    seed + 10_000 * estimator_index + 1_000 * hand_index + charge,
                )
            paired = {}
            for left, right in (("current", "independent"), ("current", "full"), ("independent", "full")):
                left_values = blocks[(left, hand, charge)]["values"]
                right_values = blocks[(right, hand, charge)]["values"]
                paired[f"{left}_minus_{right}"] = {
                    **paired_zero(left_values - right_values),
                    "batch_fluctuation_correlation": flattened_correlation(left_values, right_values),
                }
            full_rejects = estimator_scores["full"]["masked"]["bootstrap_p"] < 0.01
            current_rejects = estimator_scores["current"]["masked"]["bootstrap_p"] < 0.01
            independent_rejects = estimator_scores["independent"]["masked"]["bootstrap_p"] < 0.01
            if full_rejects:
                decision = "full_anchor_rejects_anchor_sampling_not_causal"
            elif current_rejects or independent_rejects:
                decision = "one_anchor_only_tension"
            else:
                decision = "no_estimator_rejects_at_5k"
            decisions.append(decision)
            channels[key] = {
                "estimators": estimator_scores,
                "paired_differences": paired,
                "decision": decision,
            }
    if all(value == "full_anchor_rejects_anchor_sampling_not_causal" for value in decisions):
        decision = "scorer_or_model_assembly_persists_after_full_anchor_average"
    elif any(value == "one_anchor_only_tension" for value in decisions):
        decision = "anchor_sampling_or_counter_correlation_candidate"
    else:
        decision = "5k_underpowered_or_no_rejection"
    replay = historical_replay(rows, blocks)
    if not replay["passed"]:
        raise AssertionError("current estimator does not replay the old radius4 stream")
    return {
        "schema": "matching-one/p418-anchor-paired-pilot-score/v1",
        "status": "single_authorized_5k_paired_pilot_scored",
        "issues": [418, 250],
        "freeze": str(FREEZE.relative_to(ROOT)),
        "input": {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)},
        "bootstrap": {"replicates": bootstrap, "seed": seed},
        "historical_current_estimator_replay": replay,
        "batch_covariance": {
            "batches": len(rows),
            "coordinates_per_batch": sum(len(blocks[(estimator, hand, charge)]["coordinates"]) * 2 for estimator in ESTIMATORS for hand in HANDS for charge in CHARGES),
            "full_cross_estimator_channel_covariance_recoverable_from_batch_rows": True,
        },
        "channels": channels,
        "decision": decision,
        "claim_boundary": freeze["claim_boundary"],
    }


def markdown(result: dict[str, Any]) -> str:
    lines = [
        "# P418 paired-anchor 5k pilot",
        "",
        f"Decision: `{result['decision']}`.",
        "",
        "| channel | current masked d2 (p) | independent masked d2 (p) | full masked d2 (p) | channel decision |",
        "|---|---:|---:|---:|---|",
    ]
    for key, row in result["channels"].items():
        cells = []
        for estimator in ESTIMATORS:
            score = row["estimators"][estimator]["masked"]
            cells.append(f"{score['distance_squared']:.6g} ({score['bootstrap_p']:.4g})")
        lines.append(f"| {key} | {' | '.join(cells)} | {row['decision']} |")
    replay = result["historical_current_estimator_replay"]
    lines += [
        "",
        f"The current estimator replays {replay['old_radius4_first_batches']} historical radius-4 batches with maximum summed-coordinate error `{replay['maximum_absolute_sum_error']:.3g}`.",
        "",
        "All 100 batch rows retain the joint current/independent/full × hand × charge coordinates, so the complete paired covariance remains reconstructible. This is the only authorized 5k pilot; no production extension was run.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=Path, default=DEFAULT_BATCHES)
    parser.add_argument("--bootstrap", type=int, default=100)
    parser.add_argument("--seed", type=int, default=41850510120260831)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = build_score(args.batches, args.bootstrap, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.output.with_suffix(".md").write_text(markdown(result))
    print(json.dumps({"decision": result["decision"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
