#!/usr/bin/env python3
"""Score the frozen P255 ordered-filtration proxy on an exact N650 replay."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


SOURCE_BATCH_SHA256 = "db5be1d870135053691e34605703f15e99e95df2da88dc279c0a55e26130d0af"
SOURCE_METADATA_SHA256 = "8f3d3605d5a99c043e2d465be0d2e4753d719d6388822048f6f4d365c05abbd7"
SEED = 2026102003
COUNTER = (18000000000, 18000020000)
SAMPLES = 20000
BATCHES = 100
DIVISOR = 2
CHANNELS = ("ES", "ED", "OS", "OD")
LEGACY_FIELDS = (
    "ES_num_sum", "ED_num_sum", "OS_num_sum", "OD_num_sum",
    "ambient_ES_num_sum", "ambient_ED_num_sum",
    "ambient_OS_num_sum", "ambient_OD_num_sum",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> tuple[list[dict[str, int]], list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"missing header in {path}")
        rows = [{key: int(value) for key, value in row.items()} for row in reader]
        return rows, list(reader.fieldnames)


def validate_replay(replay: list[dict[str, int]], source: list[dict[str, int]],
                    metadata: dict[str, Any]) -> list[str]:
    if metadata.get("schema") != "matching-one.p255-n650-ordered-proxy-replay.v1":
        raise ValueError("replay metadata schema mismatch")
    if (metadata.get("samples"), metadata.get("batches"), metadata.get("seed")) != (
            SAMPLES, BATCHES, SEED):
        raise ValueError("replay sample/batch/seed mismatch")
    if (metadata.get("replica_offset"), metadata.get("replica_offset") + SAMPLES) != COUNTER:
        raise ValueError("replay counter mismatch")
    if metadata.get("stored_sum_divisor") != DIVISOR:
        raise ValueError("stored divisor mismatch")
    if metadata.get("replay_source_batch_sha256") != SOURCE_BATCH_SHA256:
        raise ValueError("replay source hash declaration mismatch")
    if len(replay) != BATCHES or len(source) != BATCHES:
        raise ValueError("batch count mismatch")
    source_by_batch = {row["batch"]: row for row in source}
    if len(source_by_batch) != BATCHES:
        raise ValueError("duplicate source batch")
    for row in replay:
        batch = row["batch"]
        if batch not in source_by_batch:
            raise ValueError("unexpected replay batch")
        old = source_by_batch[batch]
        for field in (
            "counter_first", "counter_last_exclusive", "samples", *LEGACY_FIELDS
        ):
            if row[field] != old[field]:
                raise ValueError(f"legacy replay mismatch batch={batch} field={field}")
    ordered_fields = [
        field for field in replay[0]
        if field.endswith("_num_sum") and field not in LEGACY_FIELDS
    ]
    expected = [
        f"{name}_{channel}_num_sum"
        for name in metadata["ordered_fields"]
        for channel in CHANNELS
    ]
    if ordered_fields != expected:
        raise ValueError("ordered field order mismatch")
    return ordered_fields


def covariance_of_mean(batch_vectors: np.ndarray) -> np.ndarray:
    return np.cov(batch_vectors, rowvar=False, ddof=1) / batch_vectors.shape[0]


def psd_inverse(matrix: np.ndarray, rcond: float = 1e-11) -> tuple[np.ndarray, int]:
    symmetric = (matrix + matrix.T) / 2
    values, vectors = np.linalg.eigh(symmetric)
    scale = max(float(np.max(np.abs(values))), 1e-300)
    if float(np.min(values)) < -1e-8 * scale:
        raise ValueError("covariance has a material negative eigenvalue")
    keep = values > rcond * scale
    return (vectors[:, keep] / values[keep]) @ vectors[:, keep].T, int(keep.sum())


def chi2_sf_even(value: float, degrees: int) -> float:
    if degrees <= 0 or degrees % 2:
        raise ValueError("positive even degrees required")
    half = value / 2
    return math.exp(-half) * sum(
        half ** power / math.factorial(power) for power in range(degrees // 2)
    )


def score(replay_path: Path, metadata_path: Path, source_batch_path: Path,
          source_metadata_path: Path) -> dict[str, Any]:
    if sha256(source_batch_path) != SOURCE_BATCH_SHA256:
        raise ValueError("locked source batch hash mismatch")
    if sha256(source_metadata_path) != SOURCE_METADATA_SHA256:
        raise ValueError("locked source metadata hash mismatch")
    replay, _ = read_rows(replay_path)
    source, _ = read_rows(source_batch_path)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    ordered_fields = validate_replay(replay, source, metadata)
    replay.sort(key=lambda row: row["batch"])
    per_batch = replay[0]["samples"]
    if any(row["samples"] != per_batch for row in replay):
        raise ValueError("unequal batch sizes")
    vectors = np.asarray([
        [row[field] / (DIVISOR * per_batch) for field in ordered_fields]
        for row in replay
    ], dtype=float)
    point = vectors.mean(axis=0)
    covariance = covariance_of_mean(vectors)
    primary_names = [f"Cact_{channel}_num_sum" for channel in CHANNELS]
    primary_indices = [ordered_fields.index(name) for name in primary_names]
    primary = point[primary_indices]
    primary_covariance = covariance[np.ix_(primary_indices, primary_indices)]
    precision, rank = psd_inverse(primary_covariance)
    chi2 = float(primary @ precision @ primary)
    return {
        "schema": "matching-one.p255-ordered-context-proxy-score.v1",
        "status": "existing_archive_replay_scored",
        "input": {
            "replay_batch": str(replay_path),
            "replay_metadata": str(metadata_path),
            "replay_batch_sha256": sha256(replay_path),
            "replay_metadata_sha256": sha256(metadata_path),
            "source_batch_sha256": SOURCE_BATCH_SHA256,
            "source_metadata_sha256": SOURCE_METADATA_SHA256,
            "implementation_commit": metadata["git_commit"],
            "samples": SAMPLES,
            "batches": BATCHES,
            "seed": SEED,
            "counter": list(COUNTER),
            "legacy_batch_replay_exact": True,
        },
        "all_ordered_state": {
            "order": [field.removesuffix("_num_sum") for field in ordered_fields],
            "point": point.tolist(),
            "covariance_of_mean": covariance.tolist(),
        },
        "primary_Cact": {
            "order": [f"Cact_{channel}" for channel in CHANNELS],
            "point": primary.tolist(),
            "standard_error": np.sqrt(np.diag(primary_covariance)).tolist(),
            "covariance": primary_covariance.tolist(),
            "chi_square": chi2,
            "degrees_of_freedom": rank,
            "p_value": chi2_sf_even(chi2, rank),
            "frozen_gate_alpha": 0.01,
            "frozen_gate_critical_df4": 13.276704135987622,
            "signal_gate_passed": rank == 4 and chi2 > 13.276704135987622,
        },
        "decision_boundary": {
            "result_is": "same-endpoint Gaussian-join ordered-filtration proxy on an existing dependency block",
            "result_is_not": "physical Gaussian-by-annulus AU versus UA",
            "missing_real_AU_UA_state": "typed annular boundary-connectivity partition plus Gaussian lift/pushforward map",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replay", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--source-batch", type=Path, required=True)
    parser.add_argument("--source-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    options = parser.parse_args()
    result = score(
        options.replay, options.metadata,
        options.source_batch, options.source_metadata,
    )
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["primary_Cact"], indent=2))


if __name__ == "__main__":
    main()
