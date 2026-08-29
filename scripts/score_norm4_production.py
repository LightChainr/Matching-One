#!/usr/bin/env python3
"""Score the frozen norm-4 q2/Jordan closure and correlated secondary views."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from analyze_norm4_variance_pilot import (
    covariance_of_mean,
    quadratic_2,
    u_value,
)
from analyze_p48_retrospective import Histogram, project_size, read_histograms
from score_p50_fullcurve_n290 import grouped
from score_p49_fullcurve_doubling import aggregate, orientation_values


SOURCE_ORDER = (65, 130, 85, 170)
SIZE_ORDER = (65, 130, 260, 85, 170, 340)
METRICS = ("U", "P4_D", "root_gap")
PRODUCTION_COMMIT = "bfab0330f5f56ca4d746b45d737f1607e3d229a0"


@dataclass(frozen=True)
class SourceRun:
    n: int
    base_histogram: Path
    base_metadata: Path
    extension_histogram: Path
    extension_metadata: Path


@dataclass(frozen=True)
class TargetRun:
    n: int
    histogram: Path
    metadata: Path


def parse_source(text: str) -> SourceRun:
    fields = text.split(":", 4)
    if len(fields) != 5:
        raise argparse.ArgumentTypeError(
            "source run must be N:BASE_HIST:BASE_META:EXT_HIST:EXT_META"
        )
    return SourceRun(int(fields[0]), *(Path(value) for value in fields[1:]))


def parse_target(text: str) -> TargetRun:
    fields = text.split(":", 2)
    if len(fields) != 3:
        raise argparse.ArgumentTypeError("target run must be N:HISTOGRAM:METADATA")
    return TargetRun(int(fields[0]), Path(fields[1]), Path(fields[2]))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def metadata(path: Path) -> Mapping[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path}: metadata must be an object")
    return payload


def rng_signature(payload: Mapping[str, object]):
    return (
        int(payload["seed"]),
        int(payload["replica_counter_first"]),
        int(payload["replica_counter_last_exclusive"]),
        int(payload["samples_per_pair"]),
        int(payload["batches"]),
    )


def validate_provenance(sources, targets) -> None:
    base = {rng_signature(metadata(run.base_metadata)) for run in sources.values()}
    extension = {
        rng_signature(metadata(run.extension_metadata)) for run in sources.values()
    }
    if base != {(2026104501, 5000000000, 5100000000, 100000000, 100)}:
        raise ValueError("source base block differs from the frozen aligned block")
    if extension != {(2026104501, 5100000000, 7000000000, 1900000000, 100)}:
        raise ValueError("source extension differs from the frozen aligned block")
    expected_targets = {
        260: (2026105401, 8200000000, 9200000000, 1000000000, 100),
        340: (2026105402, 8200000000, 9200000000, 1000000000, 100),
    }
    for n, run in targets.items():
        payload = metadata(run.metadata)
        if rng_signature(payload) != expected_targets[n]:
            raise ValueError(f"N={n}: target block differs from the frozen block")
        if str(payload["git_commit"]) != PRODUCTION_COMMIT:
            raise ValueError(f"N={n}: target generation commit differs from freeze")
    for run in sources.values():
        payload = metadata(run.extension_metadata)
        if str(payload["git_commit"]) != PRODUCTION_COMMIT:
            raise ValueError(f"N={run.n}: source extension commit differs from freeze")


def merge_histogram_blocks(base: Path, extension: Path, n: int):
    first = read_histograms(base)
    second = read_histograms(extension)
    if set(first) != set(second) or {key[0] for key in first} != {n}:
        raise ValueError(f"N={n}: base and extension histogram keys differ")
    merged = {}
    for key in first:
        left = first[key]
        right = second[key]
        if (left.a, left.b, left.orientation, left.batch) != (
            right.a,
            right.b,
            right.orientation,
            right.batch,
        ):
            raise ValueError(f"N={n}: histogram descriptors differ")
        merged[key] = Histogram(
            n=n,
            a=left.a,
            b=left.b,
            orientation=left.orientation,
            batch=left.batch,
            samples=left.samples + right.samples,
            minus=[a + b for a, b in zip(left.minus, right.minus)],
            plus=[a + b for a, b in zip(left.plus, right.plus)],
        )
    return merged


def solve_root(function) -> float:
    lower, upper = 0.4, 0.75
    if not function(lower) <= 0.0 <= function(upper):
        raise ValueError("physical root is not bracketed by [0.4,0.75]")
    for _ in range(56):
        midpoint = (lower + upper) / 2.0
        if function(midpoint) < 0.0:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def state(by_orientation, *, omitted: int = -1) -> dict[str, float]:
    n = by_orientation["first"][0].n
    projected = project_size(by_orientation, omitted)
    rows = {
        name: aggregate(by_orientation[name], omitted)
        for name in ("first", "second")
    }
    p0 = projected["p0"]
    slopes = [orientation_values(n, rows[name], p0)["M_prime"] for name in rows]
    center_state = {
        "mean_slope": math.fsum(slopes) / 2.0,
        "P4_S_prime": float(projected["P4_S_prime"]),
    }
    roots = {
        name: solve_root(lambda p, name=name: orientation_values(n, rows[name], p)["M"])
        for name in rows
    }
    return {
        "U": u_value(n, center_state),
        "P4_D": float(projected["P4_D"]),
        "root_gap": roots["first"] - roots["second"],
        "p0": p0,
    }


def estimate_aligned(groups: Mapping[int, object], sizes: Sequence[int]):
    batch_ids = [row.batch for row in groups[sizes[0]]["first"]]
    for n in sizes[1:]:
        if [row.batch for row in groups[n]["first"]] != batch_ids:
            raise ValueError("aligned source batch ids differ")
    points = {n: state(groups[n]) for n in sizes}
    full = [points[n][metric] for n in sizes for metric in METRICS]
    pseudovalue_rows = []
    batches = len(batch_ids)
    for batch in batch_ids:
        deleted = [
            state(groups[n], omitted=batch)[metric]
            for n in sizes
            for metric in METRICS
        ]
        pseudovalue_rows.append(
            [batches * a - (batches - 1) * b for a, b in zip(full, deleted)]
        )
    return points, covariance_of_mean(pseudovalue_rows)


def estimate_one(group):
    points, covariance = estimate_aligned({group["first"][0].n: group}, [group["first"][0].n])
    return next(iter(points.values())), covariance


def assemble_covariance(source_cov, target_covariances):
    width = len(SIZE_ORDER) * len(METRICS)
    output = [[0.0] * width for _ in range(width)]
    source_positions = [SIZE_ORDER.index(n) for n in SOURCE_ORDER]
    for i, position_i in enumerate(source_positions):
        for j, position_j in enumerate(source_positions):
            for a in range(len(METRICS)):
                for b in range(len(METRICS)):
                    output[position_i * 3 + a][position_j * 3 + b] = source_cov[
                        i * 3 + a
                    ][j * 3 + b]
    for n, covariance in target_covariances.items():
        position = SIZE_ORDER.index(n)
        for a in range(3):
            for b in range(3):
                output[position * 3 + a][position * 3 + b] = covariance[a][b]
    return output


def transform_covariance(matrix, transform):
    return [
        [
            math.fsum(
                transform[i][a] * matrix[a][b] * transform[j][b]
                for a in range(len(matrix))
                for b in range(len(matrix))
            )
            for j in range(len(transform))
        ]
        for i in range(len(transform))
    ]


def score(vector, covariance, transform, label):
    residual = [math.fsum(row[i] * vector[i] for i in range(len(vector))) for row in transform]
    residual_covariance = transform_covariance(covariance, transform)
    chi_square = quadratic_2(residual, residual_covariance)
    return {
        "label": label,
        "residual": residual,
        "covariance": residual_covariance,
        "z_marginal": [
            residual[i] / math.sqrt(residual_covariance[i][i]) for i in range(2)
        ],
        "chi_square": chi_square,
        "df": 2,
        "p_value": math.exp(-chi_square / 2.0),
    }


def row(metric: str, coefficients: Sequence[float]):
    output = [0.0] * (len(SIZE_ORDER) * len(METRICS))
    offset = METRICS.index(metric)
    for position, coefficient in enumerate(coefficients):
        output[position * 3 + offset] = coefficient
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", action="append", type=parse_source, required=True)
    parser.add_argument("--target-run", action="append", type=parse_target, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sources = {run.n: run for run in args.source_run}
    targets = {run.n: run for run in args.target_run}
    if tuple(sources) != SOURCE_ORDER or tuple(targets) != (260, 340):
        raise SystemExit("source order must be 65,130,85,170 and targets 260,340")
    validate_provenance(sources, targets)

    source_groups = {}
    for n, run in sources.items():
        source_groups[n] = grouped(
            merge_histogram_blocks(run.base_histogram, run.extension_histogram, n), n
        )
    source_points, source_covariance = estimate_aligned(source_groups, SOURCE_ORDER)
    target_points = {}
    target_covariances = {}
    for n, run in targets.items():
        target_points[n], target_covariances[n] = estimate_one(
            grouped(read_histograms(run.histogram), n)
        )

    all_points = {**source_points, **target_points}
    vector = [all_points[n][metric] for n in SIZE_ORDER for metric in METRICS]
    covariance = assemble_covariance(source_covariance, target_covariances)
    q2 = [
        row("U", (1.0, -3.0, 2.0, 0.0, 0.0, 0.0)),
        row("U", (0.0, 0.0, 0.0, 1.0, -3.0, 2.0)),
    ]
    jordan = [
        row("U", (1.0, -2.0, 1.0, 0.0, 0.0, 0.0)),
        row("U", (0.0, 0.0, 0.0, 1.0, -2.0, 1.0)),
    ]
    d_factor = math.pow(4.0, -13.0 / 8.0)
    central_d = [
        row("P4_D", (-d_factor, 0.0, 1.0, 0.0, 0.0, 0.0)),
        row("P4_D", (0.0, 0.0, 0.0, -d_factor, 0.0, 1.0)),
    ]
    root_character = [
        row("root_gap", (-1.0 / 16.0, 0.0, 1.0, 0.0, 0.0, 0.0)),
        row("root_gap", (0.0, 0.0, 0.0, -1.0 / 16.0, 0.0, 1.0)),
    ]
    payload = {
        "schema": "matching-one/norm4-production-score/v1",
        "status": "prospective score under predictions/norm4_two_generator_transfer_20260829.yaml",
        "size_order": list(SIZE_ORDER),
        "metric_order_within_size": list(METRICS),
        "point": {str(n): all_points[n] for n in SIZE_ORDER},
        "primary_models_in_frozen_order": [
            score(vector, covariance, q2, "analytic_q2"),
            score(vector, covariance, jordan, "rank2_Jordan"),
        ],
        "correlated_secondary_views": [
            score(vector, covariance, central_d, "central_D_4_to_1"),
            score(vector, covariance, root_character, "root_character_plus_1_over_16"),
        ],
        "evidence_guard": (
            "The scalar primary and secondary views share histograms and are not "
            "additive evidence rows. The full-jet even-generator vector is scored separately."
        ),
        "inputs": {
            "source": [
                {
                    "N": run.n,
                    "base_histogram": str(run.base_histogram),
                    "base_sha256": sha256(run.base_histogram),
                    "base_metadata": metadata(run.base_metadata),
                    "extension_histogram": str(run.extension_histogram),
                    "extension_sha256": sha256(run.extension_histogram),
                    "extension_metadata": metadata(run.extension_metadata),
                }
                for run in sources.values()
            ],
            "target": [
                {
                    "N": run.n,
                    "histogram": str(run.histogram),
                    "sha256": sha256(run.histogram),
                    "metadata": metadata(run.metadata),
                }
                for run in targets.values()
            ],
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
