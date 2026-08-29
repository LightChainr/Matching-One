#!/usr/bin/env python3
"""Render and, after reveal, score the frozen rank-gap target predictions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Mapping, Sequence

import mpmath as mp
import yaml

from analyze_rank_gap_thermal_window import (
    jackknife_se,
    pooled_statistics,
    read_run,
    sha256,
)


DEFAULT_MANIFEST = "predictions/rank_gap_boundary_correction_targets_20260829.yaml"


def load_freeze(path: Path) -> dict:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if payload.get("status") != "source_fit_frozen_before_N325_N425_reveal":
        raise ValueError("rank-gap target freeze status changed")
    exponent = payload["model"]["exponent_in_N"]
    if (int(exponent["numerator"]), int(exponent["denominator"]), exponent["fitted"]) != (
        5, 8, False,
    ):
        raise ValueError("target model must keep the unfitted 5/8 exponent")
    if payload["freeze_state"] != {
        "target_moment_files": "unseen_zero_bytes",
        "source_only_fit": True,
        "target_values_used": False,
    }:
        raise ValueError("target freeze-state declaration changed")
    return payload


def boundary_gls(
    sizes: Sequence[int], points: Sequence[mp.mpf], covariance: Sequence[Sequence[mp.mpf]],
    target_sizes: Sequence[int],
) -> dict[str, object]:
    exponent = mp.mpf(5) / 8
    design = mp.matrix([[mp.power(n, exponent), 1] for n in sizes])
    vector = mp.matrix(points)
    inverse = mp.matrix(covariance) ** -1
    parameter_covariance = (design.T * inverse * design) ** -1
    parameters = parameter_covariance * design.T * inverse * vector
    residual = vector - design * parameters
    chi_square = (residual.T * inverse * residual)[0]
    degrees = len(sizes) - 2
    survival = mp.gammainc(mp.mpf(degrees) / 2, chi_square / 2, mp.inf) / mp.gamma(
        mp.mpf(degrees) / 2
    )
    target_design = mp.matrix([[mp.power(n, exponent), 1] for n in target_sizes])
    prediction = target_design * parameters
    prediction_covariance = target_design * parameter_covariance * target_design.T
    return {
        "parameters": parameters,
        "parameter_covariance": parameter_covariance,
        "source_residual": residual,
        "source_chi_square": chi_square,
        "source_degrees_of_freedom": degrees,
        "source_chi_square_survival": survival,
        "target_prediction": prediction,
        "target_prediction_covariance": prediction_covariance,
    }


def _close(actual: mp.mpf, frozen: object, tolerance: mp.mpf = mp.mpf("1e-25")) -> None:
    if abs(actual - mp.mpf(str(frozen))) > tolerance * max(1, abs(actual)):
        raise ValueError(f"recomputed value {actual} differs from frozen value {frozen}")


def validate_freeze(fit: Mapping[str, object], manifest: Mapping[str, object]) -> None:
    frozen = manifest["frozen_source_fit"]
    parameters = fit["parameters"]
    _close(parameters[0], frozen["A"])
    _close(parameters[1], frozen["B"])
    for actual_row, frozen_row in zip(
        fit["parameter_covariance"].tolist(), frozen["parameter_covariance"]
    ):
        for actual, expected in zip(actual_row, frozen_row):
            _close(actual, expected)
    _close(fit["source_chi_square"], frozen["source_chi_square"])
    _close(fit["source_chi_square_survival"], frozen["source_chi_square_survival"])
    targets = manifest["targets"]
    for actual, entry in zip(fit["target_prediction"], targets["entries"]):
        _close(actual, entry["frozen_gap_mean_prediction"])
    for actual_row, frozen_row in zip(
        fit["target_prediction_covariance"].tolist(),
        targets["frozen_prediction_covariance_from_source_fit"],
    ):
        for actual, expected in zip(actual_row, frozen_row):
            _close(actual, expected)


def source_fit(source_score: Mapping[str, object], manifest: Mapping[str, object]) -> dict:
    if source_score["size_order"] != manifest["model"]["source_size_order"]:
        raise ValueError("source size order differs from freeze")
    sizes = source_score["size_order"]
    points = [
        mp.mpf(source_score["by_size"][str(n)]["metrics"]["gap_mean"]["point"])
        for n in sizes
    ]
    covariance = [[mp.mpf(value) for value in row] for row in source_score["gap_mean_covariance"]]
    target_sizes = manifest["targets"]["order"]
    fit = boundary_gls(sizes, points, covariance, target_sizes)
    validate_freeze(fit, manifest)
    return fit


def score_targets(observed: Sequence[mp.mpf], observation_covariance, fit) -> dict:
    prediction = fit["target_prediction"]
    residual = mp.matrix(observed) - prediction
    total_covariance = mp.matrix(observation_covariance) + fit["target_prediction_covariance"]
    chi_square = (residual.T * total_covariance**-1 * residual)[0]
    return {
        "residual": residual,
        "total_covariance": total_covariance,
        "signed_marginal_z": [
            residual[i] / mp.sqrt(total_covariance[i, i]) for i in range(len(observed))
        ],
        "joint_chi_square": chi_square,
        "degrees_of_freedom": len(observed),
        "chi_square_survival": mp.gammainc(
            mp.mpf(len(observed)) / 2, chi_square / 2, mp.inf
        ) / mp.gamma(mp.mpf(len(observed)) / 2),
    }


def matrix_strings(matrix, digits: int = 25):
    return [[mp.nstr(value, digits) for value in row] for row in matrix.tolist()]


def render_source_fit(fit, manifest_path: Path, source_path: Path) -> dict:
    prediction = fit["target_prediction"]
    targets = [int(value) for value in yaml.safe_load(manifest_path.read_text())["targets"]["order"]]
    return {
        "schema": "matching-one/rank-gap-boundary-targets/v1",
        "status": "frozen_source_fit_targets_unseen",
        "model": "E[G]=A*N^(5/8)+B; exponent fixed, not fitted",
        "source_fit": {
            "A": mp.nstr(fit["parameters"][0], 30),
            "B": mp.nstr(fit["parameters"][1], 30),
            "parameter_covariance": matrix_strings(fit["parameter_covariance"], 30),
            "chi_square": mp.nstr(fit["source_chi_square"], 30),
            "degrees_of_freedom": fit["source_degrees_of_freedom"],
            "chi_square_survival": mp.nstr(fit["source_chi_square_survival"], 25),
        },
        "target_order": targets,
        "target_gap_mean_prediction": [mp.nstr(value, 30) for value in prediction],
        "target_scaled_amplitude_prediction": [
            mp.nstr(prediction[i] * mp.power(n, -mp.mpf(5) / 8), 30)
            for i, n in enumerate(targets)
        ],
        "target_prediction_covariance_from_source_fit": matrix_strings(
            fit["target_prediction_covariance"], 30
        ),
        "scoring_rule": (
            "r=y_target-y_frozen; Cov(r)=Cov_target_delete_one+Cov_frozen_prediction; "
            "score r^T Cov(r)^-1 r against chi-square_2 and report signed marginal z"
        ),
        "provenance": {
            "manifest": str(manifest_path),
            "manifest_sha256": sha256(manifest_path),
            "source_score": str(source_path),
            "source_score_sha256": sha256(source_path),
        },
    }


def parse_run(specification: str):
    fields = specification.split(":", 2)
    if len(fields) != 3:
        raise argparse.ArgumentTypeError("target run must be N:MOMENTS:METADATA")
    return int(fields[0]), Path(fields[1]), Path(fields[2])


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=root / DEFAULT_MANIFEST)
    parser.add_argument("--source-score", type=Path, required=True)
    parser.add_argument("--target-run", action="append", type=parse_run)
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    manifest = load_freeze(args.manifest)
    if sha256(args.source_score) != manifest["source_score_sha256"]:
        raise ValueError("source score hash differs from freeze")
    source_score = json.loads(args.source_score.read_text(encoding="utf-8"))
    fit = source_fit(source_score, manifest)
    payload = render_source_fit(fit, args.manifest, args.source_score)
    if args.target_run:
        expected = manifest["targets"]["entries"]
        if [n for n, _, _ in args.target_run] != manifest["targets"]["order"]:
            raise ValueError("target runs must follow frozen order")
        runs = [read_run(n, moments, metadata) for n, moments, metadata in args.target_run]
        observed = []
        variances = []
        for run, frozen in zip(runs, expected):
            metadata = run.metadata
            design = metadata["designs"][0]
            representations_match = all(
                [int(value) for value in design[orientation]]
                == [int(value) for value in frozen["representations"][orientation]]
                for orientation in ("first", "second")
            )
            if any([
                int(metadata["seed"]) != int(manifest["targets"]["common_seed"]),
                int(metadata["replica_counter_first"]) != int(frozen["replica_counter_first"]),
                int(metadata["replica_counter_last_exclusive"]) != int(frozen["replica_counter_last_exclusive"]),
                int(metadata["samples_per_pair"]) != int(frozen["samples_per_orientation"]),
                int(metadata["batches"]) != int(frozen["batches"]),
                not representations_match,
            ]):
                raise ValueError(f"N={run.n}: target production metadata differs from freeze")
            statistics = pooled_statistics(run)
            deleted = [pooled_statistics(run, batch)["gap_mean"] for batch in range(int(metadata["batches"]))]
            observed.append(statistics["gap_mean"])
            variances.append(jackknife_se(deleted) ** 2)
        observation_covariance = mp.diag(variances)
        scored = score_targets(observed, observation_covariance, fit)
        payload["status"] = "targets_revealed_and_scored_against_frozen_source_fit"
        payload["target_observed_gap_mean"] = [mp.nstr(value, 30) for value in observed]
        payload["target_observation_covariance"] = matrix_strings(observation_covariance, 30)
        payload["score"] = {
            "residual": [mp.nstr(value, 30) for value in scored["residual"]],
            "residual_covariance": matrix_strings(scored["total_covariance"], 30),
            "signed_marginal_z": [mp.nstr(value, 20) for value in scored["signed_marginal_z"]],
            "joint_chi_square": mp.nstr(scored["joint_chi_square"], 25),
            "degrees_of_freedom": scored["degrees_of_freedom"],
            "chi_square_survival": mp.nstr(scored["chi_square_survival"], 20),
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
