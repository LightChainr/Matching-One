#!/usr/bin/env python3
"""Analyze batch sufficient statistics from the marked-pivotal H4 pilot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Sequence


COUNT_FIELDS = (
    "samples", "primal_pivotal", "matching_pivotal", "primal_axis",
    "primal_diagonal", "primal_both", "primal_landed", "primal_h4",
    "matching_axis", "matching_diagonal", "matching_both",
    "matching_landed", "matching_h4",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_batches(path: Path) -> dict[str, list[dict[str, int]]]:
    output = {"first": [], "second": []}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"n", "a", "b", "orientation", "batch", *COUNT_FIELDS}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("batch CSV missing " + ", ".join(sorted(missing)))
        seen = set()
        for raw in reader:
            orientation = raw["orientation"]
            batch = int(raw["batch"])
            if orientation not in output or (orientation, batch) in seen:
                raise ValueError("invalid or duplicate orientation/batch")
            seen.add((orientation, batch))
            if int(raw["n"]) != 65:
                raise ValueError("pilot is frozen at N=65")
            expected = (8, 1) if orientation == "first" else (7, 4)
            if (int(raw["a"]), int(raw["b"])) != expected:
                raise ValueError("Gaussian representation differs from freeze")
            row = {"batch": batch, **{field: int(raw[field]) for field in COUNT_FIELDS}}
            for side in ("primal", "matching"):
                if not (
                    0 <= row[f"{side}_both"]
                    <= min(row[f"{side}_axis"], row[f"{side}_diagonal"])
                    <= row[f"{side}_landed"]
                    <= row[f"{side}_pivotal"]
                    <= row["samples"]
                ):
                    raise ValueError(f"{orientation} batch {batch}: invalid {side} counts")
                if row[f"{side}_h4"] != row[f"{side}_axis"] - row[f"{side}_diagonal"]:
                    raise ValueError(f"{orientation} batch {batch}: H4 identity failed")
                expected_landed = (
                    row[f"{side}_axis"] + row[f"{side}_diagonal"]
                    - row[f"{side}_both"]
                )
                if row[f"{side}_landed"] != expected_landed:
                    raise ValueError(f"{orientation} batch {batch}: landing union failed")
            output[orientation].append(row)
    for orientation, rows in output.items():
        rows.sort(key=lambda row: row["batch"])
        if [row["batch"] for row in rows] != list(range(len(rows))):
            raise ValueError(f"{orientation}: batches are incomplete")
    if len(output["first"]) != len(output["second"]):
        raise ValueError("orientations are not batch aligned")
    return output


def state(rows: Sequence[dict[str, int]], omitted: int | None = None) -> dict[str, float]:
    selected = [row for row in rows if row["batch"] != omitted]
    totals = {field: sum(row[field] for row in selected) for field in COUNT_FIELDS}
    signed = totals["primal_h4"] + totals["matching_h4"]
    landed = totals["primal_landed"] + totals["matching_landed"]
    pivotal = totals["primal_pivotal"] + totals["matching_pivotal"]
    if not landed or not pivotal:
        raise ValueError("insufficient landed pivotal events")
    return {
        "mu0_russo_control": 65.0 * pivotal / totals["samples"],
        "mu4": 65.0 * signed / totals["samples"],
        "a4": signed / landed,
        "landing_acceptance_given_pivotal": landed / pivotal,
        "primal_pivotal_probability": totals["primal_pivotal"] / totals["samples"],
        "matching_pivotal_probability": totals["matching_pivotal"] / totals["samples"],
        "landed_pivotal_events": landed,
        "signed_h4_sum": signed,
    }


def jackknife_se(values: Sequence[float]) -> float:
    mean = math.fsum(values) / len(values)
    return math.sqrt(
        (len(values) - 1) / len(values)
        * math.fsum((value - mean) ** 2 for value in values)
    )


def jackknife_covariance(vectors: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(vectors)
    means = [math.fsum(row[j] for row in vectors) / count for j in range(len(vectors[0]))]
    return [
        [
            (count - 1) / count * math.fsum(
                (row[i] - means[i]) * (row[j] - means[j]) for row in vectors
            )
            for j in range(len(means))
        ]
        for i in range(len(means))
    ]


def analyze(batch_path: Path, metadata_path: Path, exact_path: Path) -> dict:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    expected = {
        "schema": "matching-one/marked-pivotal-h4-pilot/v1",
        "N": 65,
        "representations": [[8, 1], [7, 4]],
        "radius": 3,
        "p": "0.59274605079",
        "samples": 200000,
        "batches": 100,
        "seed": 2026106201,
        "replica_counter_first": 12000000000,
        "replica_counter_last_exclusive": 12000200000,
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(f"metadata {key} differs from frozen pilot")
    if not str(metadata.get("git_commit", "")):
        raise ValueError("metadata lacks the generation commit")
    exact = json.loads(exact_path.read_text(encoding="utf-8"))
    if any(exact["symmetry_violations"].values()):
        raise ValueError("exact symmetry oracle has violations")
    for difference in ("primal_difference", "matching_difference"):
        if float(exact["russo_control"][difference]) != 0.0:
            raise ValueError("fixed-root Russo oracle failed")
    batches = read_batches(batch_path)
    count = len(batches["first"])
    full = {orientation: state(rows) for orientation, rows in batches.items()}
    deleted = {
        orientation: [state(rows, batch) for batch in range(count)]
        for orientation, rows in batches.items()
    }
    metrics = ("mu0_russo_control", "mu4", "a4", "landing_acceptance_given_pivotal")
    by_orientation = {}
    for orientation in ("first", "second"):
        by_orientation[orientation] = {
            metric: {
                "point": full[orientation][metric],
                "delete_one_batch_se": jackknife_se(
                    [row[metric] for row in deleted[orientation]]
                ),
            }
            for metric in metrics
        }
        by_orientation[orientation]["event_counts"] = {
            "landed_pivotal_events": full[orientation]["landed_pivotal_events"],
            "signed_h4_sum": full[orientation]["signed_h4_sum"],
        }
    differences = {}
    for metric in ("mu4", "a4"):
        values = [
            deleted["first"][batch][metric] - deleted["second"][batch][metric]
            for batch in range(count)
        ]
        differences[f"first_minus_second_{metric}"] = {
            "point": full["first"][metric] - full["second"][metric],
            "delete_one_batch_se": jackknife_se(values),
        }
    vector_order = ["first_mu4", "first_a4", "second_mu4", "second_a4"]
    deleted_vectors = [
        [
            deleted["first"][batch]["mu4"], deleted["first"][batch]["a4"],
            deleted["second"][batch]["mu4"], deleted["second"][batch]["a4"],
        ]
        for batch in range(count)
    ]
    return {
        "schema": "matching-one/marked-pivotal-h4-analysis/v1",
        "status": "bounded N65 engineering pilot; no exponent fit",
        "primary_observables": ["mu4", "a4"],
        "russo_control_role": "mu0 is a regression control, not independent evidence",
        "by_orientation": by_orientation,
        "paired_differences": differences,
        "covariance_order": vector_order,
        "delete_one_batch_covariance": jackknife_covariance(deleted_vectors),
        "exact_oracle": exact,
        "provenance": {
            "batches": str(batch_path), "batches_sha256": sha256(batch_path),
            "metadata": str(metadata_path), "metadata_sha256": sha256(metadata_path),
            "exact": str(exact_path), "exact_sha256": sha256(exact_path),
        },
        "interpretation_guard": (
            "The landing mark is new; the unmarked pivotal total, conditional ratio, "
            "and paired views are correlated and must not be counted additively."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--exact", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = analyze(args.batches, args.metadata, args.exact)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
