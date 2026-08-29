#!/usr/bin/env python3
"""Analyze the paired threshold-rank gap as a thermal-window observable.

The ordinary matching curve retains the marginal laws of K_minus and K_plus.
The paired gap G=K_plus-K_minus uses the joint moments emitted by the
threshold-rank engine and is therefore genuinely additional information.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import mpmath as mp
import yaml


DEFAULT_MANIFEST = "predictions/rank_gap_thermal_window_20260829.yaml"
FIELDS = (
    "sum_kminus", "sum_kplus", "sum_kminus2", "sum_kplus2",
    "sum_product", "sum_gap", "sum_gap2",
)


@dataclass(frozen=True)
class Run:
    n: int
    moments_path: Path
    metadata_path: Path
    metadata: Mapping[str, object]
    rows: Mapping[str, Sequence[Mapping[str, int]]]

    @property
    def group_key(self) -> tuple[int, int, int]:
        return (
            int(self.metadata["seed"]),
            int(self.metadata["replica_counter_first"]),
            int(self.metadata["replica_counter_last_exclusive"]),
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_manifest(path: Path) -> dict:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("rank-gap manifest must be a mapping")
    if payload.get("status") != "retrospective_observable_analysis_not_prospective":
        raise ValueError("rank-gap manifest status changed")
    exponent = payload.get("fixed_exponent_in_N", {})
    if (int(exponent.get("numerator", 0)), int(exponent.get("denominator", 0))) != (5, 8):
        raise ValueError("rank-gap exponent must remain fixed at 5/8")
    sizes = tuple(int(row["N"]) for row in payload.get("runs", ()))
    if not sizes or len(set(sizes)) != len(sizes):
        raise ValueError("rank-gap manifest requires a unique ordered size list")
    return payload


def read_run(n: int, moments_path: Path, metadata_path: Path) -> Run:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if not isinstance(metadata, dict):
        raise ValueError(f"N={n}: metadata must be an object")
    required_metadata = (
        "seed", "replica_counter_first", "replica_counter_last_exclusive",
        "samples_per_pair", "batches", "designs", "git_commit",
    )
    missing_metadata = [name for name in required_metadata if name not in metadata]
    if missing_metadata:
        raise ValueError(f"N={n}: metadata lacks {', '.join(missing_metadata)}")
    designs = metadata["designs"]
    if not isinstance(designs, list) or len(designs) != 1 or int(designs[0]["N"]) != n:
        raise ValueError(f"N={n}: metadata design mismatch")
    batches = int(metadata["batches"])
    by_orientation: dict[str, list[dict[str, int]]] = {"first": [], "second": []}
    required = {"n", "a", "b", "orientation", "batch", "samples", *FIELDS}
    with moments_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("moments CSV missing: " + ", ".join(sorted(missing)))
        seen: set[tuple[str, int]] = set()
        for raw in reader:
            if int(raw["n"]) != n:
                raise ValueError(f"N={n}: moments file contains a different size")
            orientation = raw["orientation"]
            batch = int(raw["batch"])
            key = (orientation, batch)
            if orientation not in by_orientation or key in seen:
                raise ValueError(f"N={n}: invalid or duplicate orientation/batch")
            seen.add(key)
            row = {
                "batch": batch,
                "samples": int(raw["samples"]),
                "a": int(raw["a"]),
                "b": int(raw["b"]),
                **{name: int(raw[name]) for name in FIELDS},
            }
            if row["samples"] <= 0 or row["sum_gap"] < 0 or row["sum_gap2"] < 0:
                raise ValueError(f"N={n}: invalid moment totals")
            if row["sum_gap"] != row["sum_kplus"] - row["sum_kminus"]:
                raise ValueError(f"N={n}: first gap-moment identity failed")
            expected_gap2 = (
                row["sum_kplus2"] + row["sum_kminus2"] - 2 * row["sum_product"]
            )
            if row["sum_gap2"] != expected_gap2:
                raise ValueError(f"N={n}: squared-gap identity failed")
            by_orientation[orientation].append(row)
    for orientation, rows in by_orientation.items():
        rows.sort(key=lambda row: row["batch"])
        if [row["batch"] for row in rows] != list(range(batches)):
            raise ValueError(f"N={n}: {orientation} batches are incomplete")
        expected_rep = tuple(int(value) for value in designs[0][orientation])
        if {(row["a"], row["b"]) for row in rows} != {expected_rep}:
            raise ValueError(f"N={n}: {orientation} representation mismatch")
        if sum(row["samples"] for row in rows) != int(metadata["samples_per_pair"]):
            raise ValueError(f"N={n}: {orientation} sample total mismatch")
    first_signature = [(row["batch"], row["samples"]) for row in by_orientation["first"]]
    second_signature = [(row["batch"], row["samples"]) for row in by_orientation["second"]]
    if first_signature != second_signature:
        raise ValueError(f"N={n}: orientations are not batch aligned")
    counter_count = (
        int(metadata["replica_counter_last_exclusive"])
        - int(metadata["replica_counter_first"])
    )
    if counter_count != int(metadata["samples_per_pair"]):
        raise ValueError(f"N={n}: counter interval length mismatch")
    return Run(n, moments_path, metadata_path, metadata, by_orientation)


def pooled_statistics(run: Run, omitted_batch: int | None = None) -> dict[str, mp.mpf]:
    rows = [
        row
        for orientation in ("first", "second")
        for row in run.rows[orientation]
        if row["batch"] != omitted_batch
    ]
    samples = sum(row["samples"] for row in rows)
    totals = {name: sum(row[name] for row in rows) for name in FIELDS}
    mean_minus = mp.mpf(totals["sum_kminus"]) / samples
    mean_plus = mp.mpf(totals["sum_kplus"]) / samples
    mean_gap = mp.mpf(totals["sum_gap"]) / samples
    var_minus = mp.mpf(totals["sum_kminus2"]) / samples - mean_minus**2
    var_plus = mp.mpf(totals["sum_kplus2"]) / samples - mean_plus**2
    var_gap = mp.mpf(totals["sum_gap2"]) / samples - mean_gap**2
    covariance = mp.mpf(totals["sum_product"]) / samples - mean_minus * mean_plus
    if min(var_minus, var_plus, var_gap) < 0:
        raise ValueError(f"N={run.n}: negative variance reconstructed from moments")
    correlation = covariance / mp.sqrt(var_minus * var_plus)
    cv = mp.sqrt(var_gap) / mean_gap
    return {
        "gap_mean": mean_gap,
        "gap_variance": var_gap,
        "gap_cv": cv,
        "rank_correlation": correlation,
    }


def jackknife_se(values: Sequence[mp.mpf]) -> mp.mpf:
    count = len(values)
    if count < 2:
        raise ValueError("at least two delete-one values are required")
    mean = mp.fsum(values) / count
    return mp.sqrt(
        mp.mpf(count - 1) / count
        * mp.fsum((value - mean) ** 2 for value in values)
    )


def jackknife_covariance(left: Sequence[mp.mpf], right: Sequence[mp.mpf]) -> mp.mpf:
    if len(left) != len(right) or len(left) < 2:
        raise ValueError("jackknife vectors are not aligned")
    count = len(left)
    mean_left = mp.fsum(left) / count
    mean_right = mp.fsum(right) / count
    return mp.mpf(count - 1) / count * mp.fsum(
        (x - mean_left) * (y - mean_right) for x, y in zip(left, right)
    )


def covariance_groups(runs: Sequence[Run]) -> dict[tuple[int, int, int], list[Run]]:
    groups: dict[tuple[int, int, int], list[Run]] = {}
    for run in runs:
        groups.setdefault(run.group_key, []).append(run)
    for members in groups.values():
        reference = [
            (row["batch"], row["samples"]) for row in members[0].rows["first"]
        ]
        for run in members[1:]:
            current = [(row["batch"], row["samples"]) for row in run.rows["first"]]
            if current != reference:
                raise ValueError(f"counter group {[item.n for item in members]} is not aligned")
    for index, left in enumerate(runs):
        for right in runs[index + 1:]:
            if left.group_key[0] != right.group_key[0]:
                continue
            left_first, left_last = left.group_key[1:]
            right_first, right_last = right.group_key[1:]
            if max(left_first, right_first) < min(left_last, right_last):
                if left.group_key != right.group_key:
                    raise ValueError(f"N={left.n},N={right.n}: partial counter overlap")
    return groups


def fixed_exponent_score(
    points: Sequence[mp.mpf], covariance: Sequence[Sequence[mp.mpf]], exponent: mp.mpf,
    sizes: Sequence[int],
) -> dict[str, object]:
    factors = [mp.power(n, -exponent) for n in sizes]
    scaled = [point * factor for point, factor in zip(points, factors)]
    scaled_covariance = [
        [covariance[i][j] * factors[i] * factors[j] for j in range(len(sizes))]
        for i in range(len(sizes))
    ]
    matrix = mp.matrix(scaled_covariance)
    inverse = matrix**-1
    one = mp.matrix([1 for _ in sizes])
    vector = mp.matrix(scaled)
    amplitude = (one.T * inverse * vector)[0] / (one.T * inverse * one)[0]
    residual = vector - amplitude * one
    chi_square = (residual.T * inverse * residual)[0]
    diagonal = mp.diag([scaled_covariance[i][i] for i in range(len(sizes))])
    diagonal_inverse = diagonal**-1
    diagonal_amplitude = (
        (one.T * diagonal_inverse * vector)[0] / (one.T * diagonal_inverse * one)[0]
    )
    diagonal_residual = vector - diagonal_amplitude * one
    diagonal_chi_square = (diagonal_residual.T * diagonal_inverse * diagonal_residual)[0]
    degrees = len(sizes) - 1
    survival = mp.gammainc(mp.mpf(degrees) / 2, chi_square / 2, mp.inf) / mp.gamma(
        mp.mpf(degrees) / 2
    )
    return {
        "fixed_exponent_in_N": mp.nstr(exponent, 20),
        "scaled_order": list(sizes),
        "scaled_point": [mp.nstr(value, 25) for value in scaled],
        "scaled_covariance": [
            [mp.nstr(value, 18) for value in row] for row in scaled_covariance
        ],
        "common_amplitude_gls": mp.nstr(amplitude, 25),
        "chi_square": mp.nstr(chi_square, 20),
        "degrees_of_freedom": degrees,
        "chi_square_survival": mp.nstr(survival, 15),
        "diagonal_sensitivity_amplitude": mp.nstr(diagonal_amplitude, 25),
        "diagonal_sensitivity_chi_square": mp.nstr(diagonal_chi_square, 20),
    }


def render(runs: Sequence[Run], manifest_path: Path, manifest: Mapping[str, object]) -> dict:
    sizes = [run.n for run in runs]
    expected_sizes = [int(row["N"]) for row in manifest["runs"]]
    if sizes != expected_sizes:
        raise ValueError(f"runs must follow frozen order {expected_sizes}")
    exponent = mp.mpf(5) / 8
    full = {run.n: pooled_statistics(run) for run in runs}
    deleted = {
        run.n: [
            pooled_statistics(run, int(row["batch"]))
            for row in run.rows["first"]
        ]
        for run in runs
    }
    groups = covariance_groups(runs)
    index = {n: position for position, n in enumerate(sizes)}
    covariance = [[mp.mpf(0) for _ in sizes] for _ in sizes]
    for members in groups.values():
        for left in members:
            for right in members:
                covariance[index[left.n]][index[right.n]] = jackknife_covariance(
                    [row["gap_mean"] for row in deleted[left.n]],
                    [row["gap_mean"] for row in deleted[right.n]],
                )
    points = [full[n]["gap_mean"] for n in sizes]
    by_size = {}
    for run in runs:
        metrics = {}
        for name, point in full[run.n].items():
            values = [row[name] for row in deleted[run.n]]
            metrics[name] = {
                "point": mp.nstr(point, 25),
                "delete_one_batch_se": mp.nstr(jackknife_se(values), 18),
            }
        scale = mp.power(run.n, -exponent)
        metrics["gap_mean_scaled_N_minus_5_8"] = {
            "point": mp.nstr(full[run.n]["gap_mean"] * scale, 25),
            "delete_one_batch_se": mp.nstr(
                jackknife_se([row["gap_mean"] for row in deleted[run.n]]) * scale, 18
            ),
        }
        variance_scale = mp.power(run.n, -2 * exponent)
        metrics["gap_variance_scaled_N_minus_5_4"] = {
            "point": mp.nstr(full[run.n]["gap_variance"] * variance_scale, 25),
            "delete_one_batch_se": mp.nstr(
                jackknife_se([row["gap_variance"] for row in deleted[run.n]])
                * variance_scale,
                18,
            ),
        }
        by_size[str(run.n)] = {
            "representations": {
                orientation: [
                    int(run.rows[orientation][0]["a"]),
                    int(run.rows[orientation][0]["b"]),
                ]
                for orientation in ("first", "second")
            },
            "samples_per_orientation": int(run.metadata["samples_per_pair"]),
            "batches": int(run.metadata["batches"]),
            "metrics": metrics,
        }
    transfers = []
    for transfer in manifest.get("fixed_ratio_diagnostics", ()):
        parent = int(transfer["parent"])
        child = int(transfer["child"])
        size_ratio = mp.mpf(child) / parent
        point = full[child]["gap_mean"] / full[parent]["gap_mean"] / mp.power(
            size_ratio, exponent
        )
        if runs[index[parent]].group_key != runs[index[child]].group_key:
            raise ValueError("frozen transfer diagnostic requires aligned counters")
        values = [
            deleted[child][batch]["gap_mean"]
            / deleted[parent][batch]["gap_mean"]
            / mp.power(size_ratio, exponent)
            for batch in range(len(deleted[parent]))
        ]
        transfers.append({
            "parent": parent,
            "child": child,
            "size_ratio": mp.nstr(size_ratio, 15),
            "expected_raw_gap_ratio": mp.nstr(mp.power(size_ratio, exponent), 20),
            "observed_over_expected": mp.nstr(point, 25),
            "delete_one_batch_se": mp.nstr(jackknife_se(values), 18),
        })
    score = fixed_exponent_score(points, covariance, exponent, sizes)
    return {
        "schema": "matching-one/rank-gap-thermal-window/v1",
        "status": "retrospective fixed-exponent observable bridge; exponent not fitted",
        "observable": "orientation-pooled G=K_plus-K_minus in rank units",
        "exact_joint_moment_identities": "passed for every input row",
        "size_order": sizes,
        "covariance_groups": [
            {
                "seed": key[0], "first": key[1], "last_exclusive": key[2],
                "sizes": [run.n for run in members],
                "rule": "aligned_delete_one" if len(members) > 1 else "independent",
            }
            for key, members in groups.items()
        ],
        "by_size": by_size,
        "gap_mean_covariance": [
            [mp.nstr(value, 18) for value in row] for row in covariance
        ],
        "fixed_ratio_diagnostics": transfers,
        "fixed_exponent_score": score,
        "provenance": {
            "manifest": str(manifest_path),
            "manifest_sha256": sha256(manifest_path),
            "inputs": [
                {
                    "N": run.n,
                    "moments": str(run.moments_path),
                    "moments_sha256": sha256(run.moments_path),
                    "metadata": str(run.metadata_path),
                    "metadata_sha256": sha256(run.metadata_path),
                    "source_commit": run.metadata["git_commit"],
                }
                for run in runs
            ],
        },
        "interpretation_guard": (
            "The 5/8 exponent is a no-fit bridge from nu=4/3. The common amplitude is "
            "a nuisance fit, not an exponent fit. The paired gap is not reconstructible "
            "from the two marginal threshold histograms alone."
        ),
    }


def parse_run(specification: str) -> tuple[int, Path, Path]:
    fields = specification.split(":", 2)
    if len(fields) != 3:
        raise argparse.ArgumentTypeError("run must be N:MOMENTS:METADATA")
    try:
        n = int(fields[0])
    except ValueError as exc:
        raise argparse.ArgumentTypeError("run N must be an integer") from exc
    return n, Path(fields[1]), Path(fields[2])


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", type=parse_run, required=True)
    parser.add_argument("--manifest", type=Path, default=root / DEFAULT_MANIFEST)
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    manifest = load_manifest(args.manifest)
    runs = [read_run(n, moments, metadata) for n, moments, metadata in args.run]
    payload = render(runs, args.manifest, manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
