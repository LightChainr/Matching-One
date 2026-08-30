#!/usr/bin/env python3
"""Score the frozen N260/N340 multiradius thermal-null pilot for Issue #155."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


ALPHA = 3.0 / 64.0
RADII = (2, 4, 8)
DESIGNS = {
    "N260_16_2": (260, (16, -2, 2, 16)),
    "N340_18_4": (340, (18, -4, 4, 18)),
}
BATCHES = 100
ALLOWED_SAMPLES = (20_000, 100_000)
ABS_Z_MIN = 3.0
CONDITION_MAX = 50.0
EXPANSION_FACTOR = 5.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def determinant(matrix: list[list[float]]) -> float:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def metrics(vector: list[float]) -> dict:
    matrix = [vector[:2], vector[2:]]
    det = determinant(matrix)
    norm2 = math.fsum(value * value for row in matrix for value in row)
    disc = max(0.0, norm2 * norm2 - 4.0 * det * det)
    singular = [
        math.sqrt(max(0.0, (norm2 + math.sqrt(disc)) / 2.0)),
        math.sqrt(max(0.0, (norm2 - math.sqrt(disc)) / 2.0)),
    ]
    condition = math.inf if singular[1] == 0.0 else singular[0] / singular[1]
    return {
        "matrix": matrix,
        "determinant": det,
        "singular_values": singular,
        "condition_number": condition,
    }


def jackknife_se(values: list[float]) -> float:
    mean = math.fsum(values) / len(values)
    return math.sqrt(
        (len(values) - 1) / len(values)
        * math.fsum((value - mean) ** 2 for value in values)
    )


def jackknife_covariance(vectors: list[list[float]]) -> list[list[float]]:
    count = len(vectors)
    means = [math.fsum(row[index] for row in vectors) / count
             for index in range(len(vectors[0]))]
    return [[
        (count - 1) / count * math.fsum(
            (row[first] - means[first]) * (row[second] - means[second])
            for row in vectors
        )
        for second in range(len(means))
    ] for first in range(len(means))]


def finite(value: float) -> float | None:
    return value if math.isfinite(value) else None


def read_inputs(batch_path: Path, metadata_path: Path):
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("schema") != "matching-one/c4-general-period-thermal-null-score-stream/v1":
        raise ValueError("score-stream schema mismatch")
    samples = metadata.get("samples_per_design")
    if samples not in ALLOWED_SAMPLES or metadata.get("batches") != BATCHES:
        raise ValueError("pilot permits 20k or one 100k stream with 100 batches")
    if metadata.get("radii") != list(RADII) or metadata.get("alpha_star") != "3/64":
        raise ValueError("radii or alpha changed after freeze")

    grouped = {(label, radius): [] for label in DESIGNS for radius in RADII}
    with batch_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (row["label"], int(row["radius"]))
            if key not in grouped:
                raise ValueError(f"unexpected design/radius {key}")
            size, matrix = DESIGNS[key[0]]
            observed_matrix = tuple(int(row[name]) for name in ("m00", "m01", "m10", "m11"))
            if int(row["n"]) != size or observed_matrix != matrix:
                raise ValueError(f"geometry changed for {key[0]}")
            grouped[key].append(row)
    for key, rows in grouped.items():
        rows.sort(key=lambda row: int(row["batch"]))
        if [int(row["batch"]) for row in rows] != list(range(BATCHES)):
            raise ValueError(f"{key} lacks 100 contiguous batches")
        if sum(int(row["samples"]) for row in rows) != samples:
            raise ValueError(f"sample total changed for {key}")
    for label in DESIGNS:
        for batch in range(BATCHES):
            aligned = {
                (int(grouped[(label, radius)][batch]["counter_first"]),
                 int(grouped[(label, radius)][batch]["counter_last_exclusive"]))
                for radius in RADII
            }
            if len(aligned) != 1:
                raise ValueError(f"radii not paired for {label}, batch {batch}")
    return metadata, grouped


def response_vector(rows: list[dict[str, str]]) -> list[float]:
    samples = sum(int(row["samples"]) for row in rows)
    global_t = math.fsum(int(row["global_twice_score_t"]) for row in rows) / (2 * samples)
    global_lambda = math.fsum(int(row["global_twice_score_lambda"]) for row in rows) / (2 * samples)
    local_t = math.fsum(int(row["local_twice_score_t"]) for row in rows) / (2 * samples)
    local_lambda = math.fsum(int(row["local_twice_score_lambda"]) for row in rows) / (2 * samples)
    epsilon_t = math.fsum(int(row["epsilon_sign_score_t"]) for row in rows) / (4 * samples)
    epsilon_lambda = math.fsum(int(row["epsilon_sign_score_lambda"]) for row in rows) / (4 * samples)
    return [
        global_t,
        global_lambda,
        local_t + ALPHA * epsilon_t,
        local_lambda + ALPHA * epsilon_lambda,
    ]


def analyze_rows(rows: list[dict[str, str]]) -> dict:
    point_vector = response_vector(rows)
    point = metrics(point_vector)
    deleted_vectors = [response_vector(rows[:batch] + rows[batch + 1:])
                       for batch in range(BATCHES)]
    deleted_metrics = [metrics(vector) for vector in deleted_vectors]
    covariance = jackknife_covariance(deleted_vectors)
    response_se = [math.sqrt(max(0.0, covariance[index][index])) for index in range(4)]
    determinant_se = jackknife_se([row["determinant"] for row in deleted_metrics])
    determinant_z = point["determinant"] / determinant_se if determinant_se else math.inf
    lambda_z = point_vector[3] / response_se[3] if response_se[3] else math.inf
    passes = (
        abs(lambda_z) >= ABS_Z_MIN
        and abs(determinant_z) >= ABS_Z_MIN
        and point["condition_number"] <= CONDITION_MAX
    )
    samples = sum(int(row["samples"]) for row in rows)
    fisher = {
        "tt": math.fsum(int(row["sum_score_t2"]) for row in rows) / samples,
        "lambda_lambda": math.fsum(int(row["sum_score_lambda2"]) for row in rows) / samples,
        "t_lambda": math.fsum(int(row["sum_score_cross"]) for row in rows) / samples,
    }
    return {
        "samples": samples,
        "response_order": [
            "global_cross_d_t", "global_cross_d_lambda",
            "thermal_null_d_t", "thermal_null_d_lambda",
        ],
        "point": {
            "matrix": point["matrix"],
            "determinant": point["determinant"],
            "singular_values": point["singular_values"],
            "condition_number": finite(point["condition_number"]),
        },
        "response_covariance": covariance,
        "response_standard_errors": response_se,
        "thermal_null_lambda_z": finite(lambda_z),
        "determinant_standard_error": determinant_se,
        "determinant_z": finite(determinant_z),
        "gate": {
            "absolute_thermal_null_lambda_z_min": ABS_Z_MIN,
            "absolute_determinant_z_min": ABS_Z_MIN,
            "condition_number_max": CONDITION_MAX,
            "passes": passes,
        },
        "empirical_fisher": fisher,
        "exact_fisher": {"tt": 4 * int(rows[0]["n"]),
                         "lambda_lambda": 4 * int(rows[0]["n"]),
                         "t_lambda": 0},
        "delete_one": [
            {
                "excluded_batch": batch,
                "matrix": row["matrix"],
                "determinant": row["determinant"],
                "singular_values": row["singular_values"],
                "condition_number": finite(row["condition_number"]),
            }
            for batch, row in enumerate(deleted_metrics)
        ],
    }


def expansion_decision(metadata: dict, results: dict) -> dict:
    if metadata["samples_per_design"] != 20_000:
        return {"eligible": False, "reason": "already_expanded_or_not_20k"}
    if float(metadata["elapsed_seconds"]) >= 120.0:
        return {"eligible": False, "reason": "20k_runtime_not_below_120_seconds"}
    primary = [results[label]["8"] for label in DESIGNS]
    current_pass = all(row["gate"]["passes"] for row in primary)
    projected = []
    for row in primary:
        projected.append({
            "label": next(label for label in DESIGNS if results[label]["8"] is row),
            "projected_abs_lambda_z": math.sqrt(EXPANSION_FACTOR)
                * abs(row["thermal_null_lambda_z"] or 0.0),
            "projected_abs_determinant_z": math.sqrt(EXPANSION_FACTOR)
                * abs(row["determinant_z"] or 0.0),
            "central_condition_number": row["point"]["condition_number"],
        })
    capable = all(
        row["projected_abs_lambda_z"] >= ABS_Z_MIN
        and row["projected_abs_determinant_z"] >= ABS_Z_MIN
        and row["central_condition_number"] is not None
        and row["central_condition_number"] <= CONDITION_MAX
        for row in projected
    )
    return {
        "eligible": capable and not current_pass,
        "rule": (
            "expand once to 100k iff runtime<120s, R8 central condition<=50 at both sizes, "
            "and sqrt(5)-projected |lambda z| and |det z| both reach 3 at both sizes"
        ),
        "already_passes_without_expansion": current_pass,
        "projection": projected,
        "reason": "clear_variance_limited_candidate" if capable and not current_pass
                  else "central_condition_or_sqrtN_projection_cannot_cross_both_size_gates",
    }


def render(batch_path: Path, metadata_path: Path) -> dict:
    metadata, grouped = read_inputs(batch_path, metadata_path)
    results = {
        label: {str(radius): analyze_rows(grouped[(label, radius)]) for radius in RADII}
        for label in DESIGNS
    }
    return {
        "schema": "matching-one/p155-general-period-thermal-null-pilot/v1",
        "issue": 155,
        "design": {
            "sizes": [
                {"label": label, "N": size, "period_matrix": [[matrix[0], matrix[1]], [matrix[2], matrix[3]]]}
                for label, (size, matrix) in DESIGNS.items()
            ],
            "radii": list(RADII),
            "cutoff": "euclidean",
            "alpha_star": "3/64",
            "primary_radius": 8,
            "diagnostic_radii": [2, 4],
            "observable_rows": [
                "global_cross_half_difference",
                "local_pivotal_H4_half_difference + (3/64) epsilon_cell",
            ],
            "score_columns": ["S_t", "S_lambda"],
        },
        "results": results,
        "primary_passes_both_sizes": all(
            results[label]["8"]["gate"]["passes"] for label in DESIGNS
        ),
        "expansion_decision": expansion_decision(metadata, results),
        "claim_boundary": {
            "answers": (
                "whether the frozen UV thermal-null counterterm exposes a nonzero staggered "
                "response and a conditioned second response direction at R=8"
            ),
            "does_not_answer": [
                "whether the finite-R readout remains exactly thermal-null away from N=10,R=1",
                "whether the second direction is an RG eigenoperator or has a CFT field identity",
            ],
        },
        "inputs": {
            "batches": str(batch_path),
            "batches_sha256": sha256(batch_path),
            "metadata": str(metadata_path),
            "metadata_sha256": sha256(metadata_path),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render(args.batches, args.metadata)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
