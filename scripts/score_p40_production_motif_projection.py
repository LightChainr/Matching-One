#!/usr/bin/env python3
"""Frozen out-of-fold score for Issue #40's paired four-motif projection."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable, Sequence

from control_variate_estimator import _solve, invertibility_test


LABELS = ("delta_q", "delta_nn_edge", "delta_diagonal_pair", "delta_face", "delta_right_angle")
FOLDS = 5


def read_jsonl(paths: Iterable[Path]) -> list[dict]:
    rows: list[dict] = []
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            rows.extend(json.loads(line) for line in handle if line.strip())
    if not rows:
        raise ValueError("motif inputs are empty")
    return rows


def _transform(names: Sequence[str]) -> list[list[float]]:
    lookup = {name: index for index, name in enumerate(names)}
    required = ("q", "E_mc", "F0_mc", "nnn_pos_mc", "nnn_neg_mc", "right_angle_mc")
    missing = [name for name in required if name not in lookup]
    if missing:
        raise ValueError("missing frozen observables: " + ",".join(missing))
    transform = [[0.0] * len(names) for _ in LABELS]
    transform[0][lookup["q"]] = 1.0
    transform[1][lookup["E_mc"]] = 1.0
    transform[2][lookup["nnn_pos_mc"]] = 1.0
    transform[2][lookup["nnn_neg_mc"]] = 1.0
    transform[3][lookup["F0_mc"]] = 1.0
    transform[4][lookup["right_angle_mc"]] = 1.0
    return transform


def _matvec(transform: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    return [math.fsum(weight * value for weight, value in zip(row, vector)) for row in transform]


def _congruence(transform: Sequence[Sequence[float]], matrix: Sequence[Sequence[float]]) -> list[list[float]]:
    return [[
        math.fsum(
            transform[left][i] * matrix[i][j] * transform[right][j]
            for i in range(len(matrix)) for j in range(len(matrix))
        )
        for right in range(len(transform))]
        for left in range(len(transform))]


def batch_moments(row: dict) -> dict:
    names = list(row["names"])
    width = len(names)
    cross = row.get("cross_gram")
    if cross is None or row.get("cross_gram_semantics") != "sum first[i]*second[j] over same replicas":
        raise ValueError("joint first-second replica Gram is required")
    if len(cross) != width or any(len(item) != width for item in cross):
        raise ValueError("cross-Gram width mismatch")
    first, second = row["first"], row["second"]
    delta_sum = [float(a) - float(b) for a, b in zip(first["sum"], second["sum"])]
    delta_gram = [[
        float(first["gram"][i][j]) + float(second["gram"][i][j])
        - float(cross[i][j]) - float(cross[j][i])
        for j in range(width)] for i in range(width)]
    transform = _transform(names)
    return {
        "n": int(row["n"]),
        "batch": int(row["batch"]),
        "samples": int(row["samples"]),
        "sum": _matvec(transform, delta_sum),
        "gram": _congruence(transform, delta_gram),
        "first": [int(first["a"]), int(first["b"])],
        "second": [int(second["a"]), int(second["b"])],
        "identity_l1": float(first["identity_l1"]) + float(second["identity_l1"]),
        "wrapping_l1": float(first["wrapping_l1"]) + float(second["wrapping_l1"]),
    }


def pool(moments: Sequence[dict]) -> dict:
    if not moments:
        raise ValueError("cannot pool an empty moment set")
    width = len(LABELS)
    samples = sum(item["samples"] for item in moments)
    total = [math.fsum(item["sum"][i] for item in moments) for i in range(width)]
    gram = [[math.fsum(item["gram"][i][j] for item in moments) for j in range(width)] for i in range(width)]
    means = [value / samples for value in total]
    covariance = [[
        (gram[i][j] - total[i] * total[j] / samples) / (samples - 1)
        for j in range(width)] for i in range(width)]
    return {"samples": samples, "sum": total, "gram": gram, "means": means, "covariance": covariance}


def fit_beta(covariance: Sequence[Sequence[float]]) -> list[float]:
    controls = [[covariance[i][j] for j in range(1, 5)] for i in range(1, 5)]
    invertibility_test(controls)
    return _solve(controls, [covariance[i][0] for i in range(1, 5)])


def adjusted_variance(covariance: Sequence[Sequence[float]], beta: Sequence[float]) -> float:
    weights = [1.0, *(-value for value in beta)]
    return math.fsum(
        weights[i] * covariance[i][j] * weights[j]
        for i in range(5) for j in range(5)
    )


def score_size(moments: Sequence[dict]) -> dict:
    ordered = sorted(moments, key=lambda item: item["batch"])
    batches = [item["batch"] for item in ordered]
    if batches != list(range(len(ordered))):
        raise ValueError("batches must be the contiguous frozen range starting at zero")
    if len(ordered) < 2 * FOLDS:
        raise ValueError("too few batches for frozen five-fold score")
    full = pool(ordered)
    fold_rows = []
    adjusted_sum = 0.0
    weighted_raw_variance = 0.0
    weighted_adjusted_variance = 0.0
    for fold in range(FOLDS):
        evaluation_rows = [item for item in ordered if item["batch"] % FOLDS == fold]
        training_rows = [item for item in ordered if item["batch"] % FOLDS != fold]
        training = pool(training_rows)
        evaluation = pool(evaluation_rows)
        beta = fit_beta(training["covariance"])
        raw_variance = evaluation["covariance"][0][0]
        residual_variance = adjusted_variance(evaluation["covariance"], beta)
        eval_adjusted_sum = evaluation["sum"][0] - math.fsum(
            value * evaluation["sum"][index] for index, value in enumerate(beta, start=1)
        )
        adjusted_sum += eval_adjusted_sum
        weighted_raw_variance += evaluation["samples"] * raw_variance
        weighted_adjusted_variance += evaluation["samples"] * residual_variance
        fold_rows.append({
            "fold": fold,
            "training_batches": len(training_rows),
            "evaluation_batches": len(evaluation_rows),
            "training_samples": training["samples"],
            "evaluation_samples": evaluation["samples"],
            "beta": beta,
            "heldout_raw_variance_per_replica": raw_variance,
            "heldout_residual_variance_per_replica": residual_variance,
            "heldout_projection_variance_per_replica": raw_variance - residual_variance,
            "heldout_variance_reduction": raw_variance / residual_variance if residual_variance > 0 else None,
            "heldout_adjusted_mean": eval_adjusted_sum / evaluation["samples"],
        })
    total_samples = full["samples"]
    raw_variance = weighted_raw_variance / total_samples
    residual_variance = weighted_adjusted_variance / total_samples
    return {
        "N": moments[0]["n"],
        "geometries": {"first": moments[0]["first"], "second": moments[0]["second"]},
        "samples": total_samples,
        "batches": len(ordered),
        "labels": list(LABELS),
        "full_replica_means": full["means"],
        "full_target_control_covariance": full["covariance"],
        "exact_gates": {
            "identity_l1": math.fsum(item["identity_l1"] for item in ordered),
            "wrapping_l1": math.fsum(item["wrapping_l1"] for item in ordered),
            "control_mean_max_abs": max(abs(value) for value in full["means"][1:]),
        },
        "cross_fit": {
            "fold_rule": "batch modulo 5; fit on four folds and score only the omitted fold",
            "folds": fold_rows,
            "oof_adjusted_mean": adjusted_sum / total_samples,
            "oof_raw_variance_per_replica": raw_variance,
            "oof_residual_variance_per_replica": residual_variance,
            "oof_projection_variance_per_replica": raw_variance - residual_variance,
            "oof_unexplained_fraction": residual_variance / raw_variance,
            "oof_variance_reduction": raw_variance / residual_variance,
            "oof_conditional_standard_error": math.sqrt(weighted_adjusted_variance) / total_samples,
            "boundary": "all variance ratios are held-out fold scores; no training covariance is reported as performance",
        },
    }


def build_report(rows: Sequence[dict], metadata_paths: Sequence[Path]) -> dict:
    moments = [batch_moments(row) for row in rows]
    by_n: dict[int, list[dict]] = {}
    for item in moments:
        by_n.setdefault(item["n"], []).append(item)
    metadata = [json.loads(path.read_text(encoding="utf-8")) for path in metadata_paths]
    return {
        "schema": "matching-one/p40-production-motif-projection-score/v1",
        "issue": 40,
        "status": "post-freeze production score",
        "control_contract": {
            "order": list(LABELS[1:]),
            "definitions": {
                "delta_nn_edge": "Delta E_mc",
                "delta_diagonal_pair": "Delta(nnn_pos_mc+nnn_neg_mc)",
                "delta_face": "Delta F0_mc",
                "delta_right_angle": "Delta right_angle_mc for translated {i,i+x,i+y}",
            },
            "fixed_K_conditional_means": "exactly zero for every control contrast",
        },
        "metadata": metadata,
        "by_N": {str(n): score_size(items) for n, items in sorted(by_n.items())},
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--motifs", type=Path, nargs="+", required=True)
    parser.add_argument("--metadata", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    report = build_report(read_jsonl(args.motifs), args.metadata)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("wrote " + str(args.output))
    for n, item in report["by_N"].items():
        score = item["cross_fit"]
        print(f"N={n} OOF VR={score['oof_variance_reduction']:.6g} residual={score['oof_residual_variance_per_replica']:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
