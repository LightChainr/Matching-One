#!/usr/bin/env python3
"""Score the frozen shared-prefix branching-continuation pilot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

SIZES = ("N325", "N425")
ORIENTATIONS = ("first", "second")
ENVIRONMENTS = tuple((size, orientation) for size in SIZES
                     for orientation in ORIENTATIONS)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def covariance(delete_one):
    values = np.asarray(delete_one, dtype=float)
    mean = values.mean(axis=0)
    centered = values - mean
    return (len(values) - 1) / len(values) * centered.T @ centered


def load_rows(path: Path):
    output = []
    with path.open(newline="") as stream:
        for raw in csv.DictReader(stream):
            output.append({
                "size": f"N{int(raw['n'])}", "orientation": raw["orientation"],
                "batch": int(raw["batch"]), "replica": int(raw["replica"]),
                "b1_safe": int(raw["checkpoint_b1_safe_count"]),
                "common_safe": int(raw["branch_common_safe"]),
                "y1": int(raw["branch_clone1_survives"]),
                "y2": int(raw["branch_clone2_survives"]),
                "both": int(raw["branch_both_survive"]),
            })
    return output


def environment_estimate(rows):
    if not rows:
        raise ValueError("no at-risk rows")
    y1 = np.asarray([row["y1"] for row in rows], dtype=float)
    y2 = np.asarray([row["y2"] for row in rows], dtype=float)
    both = np.asarray([row["both"] for row in rows], dtype=float)
    b1 = np.asarray([row["b1_safe"] for row in rows], dtype=float)
    common = np.asarray([row["common_safe"] for row in rows], dtype=float)
    b2 = float(np.mean((y1 + y2) / 2.0))
    gap = float(np.mean(both) - np.mean(y1) * np.mean(y2))
    safe = common == 1
    conditional_gap = float(
        np.mean(both[safe]) - np.mean(y1[safe]) * np.mean(y2[safe]))
    return {
        "at_risk_rows": len(rows),
        "b1_safe_mean": float(np.mean(b1)),
        "defined_rate_given_risk": float(np.mean(common)),
        "clone1_survival": float(np.mean(y1)),
        "clone2_survival": float(np.mean(y2)),
        "suffix_symmetry_difference": float(np.mean(y1) - np.mean(y2)),
        "b2_survival_estimate": b2,
        "branch_success": float(np.mean(both)),
        "clone_dependence_gap": gap,
        "successor_heterogeneity_gap_given_common_safe": conditional_gap,
    }


def vector(grouped, omit_size=None, omit_batch=None):
    estimates = {}
    for key in ENVIRONMENTS:
        rows = [row for row in grouped[key]
                if not (key[0] == omit_size and row["batch"] == omit_batch)]
        estimates[key] = environment_estimate(rows)
    return ([estimates[key]["b2_survival_estimate"] for key in ENVIRONMENTS] +
            [estimates[key]["clone_dependence_gap"] for key in ENVIRONMENTS],
            estimates)


def common_gap(gaps, cov):
    cov = np.asarray(cov, dtype=float)
    gaps = np.asarray(gaps, dtype=float)
    inverse = np.linalg.pinv(cov, rcond=1e-12)
    ones = np.ones(len(gaps))
    precision = float(ones @ inverse @ ones)
    mean = float(ones @ inverse @ gaps / precision)
    se = math.sqrt(1.0 / precision)
    return {"estimate": mean, "se": se, "z": mean / se,
            "positive_3sigma": mean > 0 and mean / se >= 3.0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", type=Path, required=True)
    parser.add_argument("--n325-csv", type=Path, required=True)
    parser.add_argument("--n325-metadata", type=Path, required=True)
    parser.add_argument("--n425-csv", type=Path, required=True)
    parser.add_argument("--n425-metadata", type=Path, required=True)
    parser.add_argument("--runner-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    freeze = json.loads(args.freeze.read_text())
    runner = Path(__file__).parents[1] / "src/threshold_rank_integer_period_mc.cpp"
    if sha256(runner) != freeze["runner_source_sha256"]:
        raise ValueError("runner source differs from frozen source")

    grouped = {key: [] for key in ENVIRONMENTS}
    inputs = {}
    for size, csv_path, metadata_path in (
            ("N325", args.n325_csv, args.n325_metadata),
            ("N425", args.n425_csv, args.n425_metadata)):
        contract = freeze["runs"][size]
        metadata = json.loads(metadata_path.read_text())
        checks = {
            "commit": metadata["git_commit"] == args.runner_commit,
            "branching_clones": metadata["branching_clones"] is True,
            "samples": metadata["samples_per_pair"] == contract["samples"],
            "batches": metadata["batches"] == contract["batches"],
            "seed": metadata["seed"] == contract["seed"],
            "replica_first": metadata["replica_counter_first"] ==
                             contract["replica_counter_first"],
            "replica_last": metadata["replica_counter_last_exclusive"] ==
                            contract["replica_counter_last_exclusive"],
            "k0": metadata["geometry_pilot_k0"] == contract["k0"],
        }
        if not all(checks.values()):
            raise ValueError(f"{size} metadata contract failed: {checks}")
        rows = load_rows(csv_path)
        for orientation in ORIENTATIONS:
            grouped[(size, orientation)] = [row for row in rows
                                             if row["orientation"] == orientation]
        inputs[size] = {"csv_sha256": sha256(csv_path),
                        "metadata_sha256": sha256(metadata_path),
                        "at_risk_rows": len(rows), "checks": checks}

    full_vector, estimates = vector(grouped)
    cov = np.zeros((8, 8))
    for size in SIZES:
        deleted = [vector(grouped, size, batch)[0]
                   for batch in range(freeze["runs"][size]["batches"])]
        cov += covariance(deleted)
    standard_errors = np.sqrt(np.maximum(np.diag(cov), 0.0))

    by_environment = {}
    for index, key in enumerate(ENVIRONMENTS):
        item = dict(estimates[key])
        item["risk_rate"] = item["at_risk_rows"] / freeze["runs"][key[0]]["samples"]
        item["b2_survival_se"] = float(standard_errors[index])
        item["clone_dependence_gap_se"] = float(standard_errors[4 + index])
        item["clone_dependence_gap_z"] = (
            item["clone_dependence_gap"] / item["clone_dependence_gap_se"])
        by_environment[f"{key[0]}_{key[1]}"] = item

    size_common = {}
    for si, size in enumerate(SIZES):
        gap = full_vector[4 + 2 * si:4 + 2 * si + 2]
        block = cov[4 + 2 * si:4 + 2 * si + 2,
                    4 + 2 * si:4 + 2 * si + 2]
        size_common[size] = common_gap(gap, block)
    all_positive = all(value > 0 for value in full_vector[4:])
    extend = all_positive and all(size_common[size]["positive_3sigma"] for size in SIZES)

    conditional_vector = [
        estimates[key]["successor_heterogeneity_gap_given_common_safe"]
        for key in ENVIRONMENTS]
    conditional_cov = np.zeros((4, 4))
    for size in SIZES:
        deleted = []
        for batch in range(freeze["runs"][size]["batches"]):
            _, values = vector(grouped, size, batch)
            deleted.append([
                values[key]["successor_heterogeneity_gap_given_common_safe"]
                for key in ENVIRONMENTS])
        conditional_cov += covariance(deleted)

    payload = {
        "schema": "matching-one/p429-branching-continuation-score/v1",
        "freeze_sha256": sha256(args.freeze), "runner_commit": args.runner_commit,
        "inputs": inputs,
        "vector_order": ([f"b2_survival:{s}:{o}" for s, o in ENVIRONMENTS] +
                         [f"clone_gap:{s}:{o}" for s, o in ENVIRONMENTS]),
        "vector": full_vector, "batch_jackknife_covariance": cov.tolist(),
        "standard_errors": standard_errors.tolist(),
        "environments": by_environment, "size_common_gap": size_common,
        "secondary_successor_heterogeneity": {
            "definition": "E[y1*y2|common safe]-E[y1|common safe]E[y2|common safe]",
            "vector_order": [f"conditional_gap:{s}:{o}" for s, o in ENVIRONMENTS],
            "vector": conditional_vector,
            "batch_jackknife_covariance": conditional_cov.tolist(),
            "standard_errors": np.sqrt(
                np.maximum(np.diag(conditional_cov), 0.0)).tolist(),
            "claim_boundary": "secondary decomposition; the frozen unconditional gap remains primary"
        },
        "extension_decision": "extend_to_100k" if extend else "stop_at_20k",
        "claim_boundary": freeze["claim_boundary"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"b2": full_vector[:4], "gap": full_vector[4:],
                      "size_common_gap": size_common,
                      "extension_decision": payload["extension_decision"]}, indent=2))


if __name__ == "__main__":
    main()
