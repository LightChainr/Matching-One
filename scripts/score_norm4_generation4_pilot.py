#!/usr/bin/env python3
"""Score the frozen N520/N680 fourth-generation even-mode recurrence."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

from analyze_p48_retrospective import read_histograms
from score_p50_fullcurve_n290 import generalized_covariance_score, grouped
from score_norm4_production import (
    METRICS,
    estimate_aligned as scalar_estimate_aligned,
    estimate_one as scalar_estimate_one,
    merge_histogram_blocks as scalar_merge_histograms,
)
from score_norm4_thermal_jet import (
    ORDERS,
    estimate_aligned as jet_estimate_aligned,
    grouped_histograms,
    merge_histogram_blocks as jet_merge_histograms,
    read_gap_batches,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE_ORDER = (65, 130, 85, 170)
SIZE_ORDER = (65, 130, 260, 520, 85, 170, 340, 680)
LINEAGE_POSITIONS = ((0, 1, 2, 3), (4, 5, 6, 7))
PRODUCTION_COMMIT = "bfab0330f5f56ca4d746b45d737f1607e3d229a0"
EXPECTED_NEW = {
    520: {
        "seed": 2026154520,
        "first": [22, 6],
        "second": [18, 14],
        "smith": [2, 260],
    },
    680: {
        "seed": 2026154680,
        "first": [26, 2],
        "second": [22, 14],
        "smith": [2, 340],
    },
}
SOURCE_FILES = {
    65: (
        "results/server-20260828/P45-root-amplitude/n65",
        "results/server-20260829/P154-norm4-production/raw/n65_1900m",
    ),
    130: (
        "results/server-20260828/P49-fullcurve-doubling-100m/raw/n130",
        "results/server-20260829/P154-norm4-production/raw/n130_1900m",
    ),
    85: (
        "results/server-20260828/P45-root-amplitude/n85",
        "results/server-20260829/P154-norm4-production/raw/n85_1900m",
    ),
    170: (
        "results/server-20260828/P49-fullcurve-doubling-100m/raw/n170",
        "results/server-20260829/P154-norm4-production/raw/n170_1900m",
    ),
}
OLD_TARGET_FILES = {
    260: "results/server-20260829/P154-norm4-production/raw/n260_1b",
    340: "results/server-20260829/P154-norm4-production/raw/n340_1b",
}


@dataclass(frozen=True)
class NewTarget:
    n: int
    histogram: Path
    moments: Path
    metadata: Path


def parse_target(text: str) -> NewTarget:
    fields = text.split(":", 3)
    if len(fields) != 4:
        raise argparse.ArgumentTypeError("target must be N:HIST:MOMENTS:METADATA")
    return NewTarget(int(fields[0]), *(Path(value) for value in fields[1:]))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_new_target(run: NewTarget) -> Mapping[str, object]:
    expected = EXPECTED_NEW[run.n]
    payload = json.loads(run.metadata.read_text(encoding="utf-8"))
    signature = (
        int(payload["samples_per_pair"]),
        int(payload["batches"]),
        int(payload["seed"]),
        int(payload["replica_counter_first"]),
        int(payload["replica_counter_last_exclusive"]),
        str(payload["git_commit"]),
    )
    wanted = (
        100000000,
        100,
        expected["seed"],
        9300000000,
        9400000000,
        PRODUCTION_COMMIT,
    )
    if signature != wanted:
        raise ValueError(f"N={run.n}: production signature differs from freeze")
    design = payload["designs"][0]
    if (
        design["N"] != run.n
        or design["first"] != expected["first"]
        or design["second"] != expected["second"]
        or design["first_smith_invariants"] != expected["smith"]
        or design["second_smith_invariants"] != expected["smith"]
    ):
        raise ValueError(f"N={run.n}: geometry/order differs from freeze")
    return payload


def place_covariance(source_cov, independent, width: int):
    dimension = len(SIZE_ORDER) * width
    output = [[0.0] * dimension for _ in range(dimension)]
    for i, n_i in enumerate(SOURCE_ORDER):
        for j, n_j in enumerate(SOURCE_ORDER):
            p_i, p_j = SIZE_ORDER.index(n_i), SIZE_ORDER.index(n_j)
            for a in range(width):
                for b in range(width):
                    output[p_i * width + a][p_j * width + b] = source_cov[
                        i * width + a
                    ][j * width + b]
    for n, covariance in independent.items():
        position = SIZE_ORDER.index(n)
        for a in range(width):
            for b in range(width):
                output[position * width + a][position * width + b] = covariance[a][b]
    return output


def recurrence_matrix(width: int, eigenvalue: float):
    coefficients = (-eigenvalue, 1 + 2 * eigenvalue, -2 - eigenvalue, 1.0)
    matrix = []
    for positions in LINEAGE_POSITIONS:
        for component in range(width):
            row = [0.0] * (len(SIZE_ORDER) * width)
            for position, coefficient in zip(positions, coefficients):
                row[position * width + component] = coefficient
            matrix.append(row)
    return matrix


def transformed_score(vector, covariance, width: int, eigenvalue: float, labels):
    matrix = recurrence_matrix(width, eigenvalue)
    residual = [
        sum(row[index] * vector[index] for index in range(len(vector)))
        for row in matrix
    ]
    residual_covariance = [
        [
            sum(
                matrix[i][a] * covariance[a][b] * matrix[j][b]
                for a in range(len(vector))
                for b in range(len(vector))
            )
            for j in range(len(matrix))
        ]
        for i in range(len(matrix))
    ]
    return {
        "secondary_eigenvalue": eigenvalue,
        "residual_order": labels,
        "residual": residual,
        "covariance": residual_covariance,
        "score": generalized_covariance_score(residual, residual_covariance),
    }


def render(new_targets: Sequence[NewTarget], prediction_path: Path) -> dict:
    targets = {run.n: run for run in new_targets}
    if tuple(targets) != (520, 680):
        raise ValueError("new target order must be N520 then N680")
    metadata = {n: validate_new_target(run) for n, run in targets.items()}
    prediction = json.loads(prediction_path.read_text(encoding="utf-8"))
    if prediction.get("schema") != "matching-one/norm4-even-rank2-transfer/v1":
        raise ValueError("unexpected prediction schema")
    eigenvalue = float(prediction["future_branch"]["frozen_secondary_eigenvalue"])
    if eigenvalue != 0.5:
        raise ValueError("generation-4 scorer is frozen to lambda=1/2")

    scalar_source_groups = {}
    jet_source_groups = {}
    jet_source_gaps = {}
    for n, (base_text, extension_text) in SOURCE_FILES.items():
        base, extension = ROOT / base_text, ROOT / extension_text
        scalar_source_groups[n] = grouped(
            scalar_merge_histograms(
                base.with_suffix(".hist.csv"), extension.with_suffix(".hist.csv"), n
            ),
            n,
        )
        jet_source_groups[n] = grouped_histograms(
            jet_merge_histograms(
                (base.with_suffix(".hist.csv"), extension.with_suffix(".hist.csv")), n
            ),
            n,
        )
        jet_source_gaps[n] = read_gap_batches(
            (base.with_suffix(".moments.csv"), extension.with_suffix(".moments.csv")), n
        )
    scalar_source_points, scalar_source_cov = scalar_estimate_aligned(
        scalar_source_groups, SOURCE_ORDER
    )
    jet_source_points, jet_source_cov = jet_estimate_aligned(
        jet_source_groups, jet_source_gaps, SOURCE_ORDER
    )

    scalar_independent_points = {}
    scalar_independent_cov = {}
    jet_independent_points = {}
    jet_independent_cov = {}
    paths = {
        **{n: (ROOT / text).with_suffix(".hist.csv") for n, text in OLD_TARGET_FILES.items()},
        **{n: run.histogram for n, run in targets.items()},
    }
    moments = {
        **{n: (ROOT / text).with_suffix(".moments.csv") for n, text in OLD_TARGET_FILES.items()},
        **{n: run.moments for n, run in targets.items()},
    }
    for n in (260, 340, 520, 680):
        scalar_point, scalar_cov = scalar_estimate_one(
            grouped(read_histograms(paths[n]), n)
        )
        scalar_independent_points[n] = scalar_point
        scalar_independent_cov[n] = scalar_cov
        jet_group = grouped_histograms(jet_merge_histograms((paths[n],), n), n)
        jet_gap = read_gap_batches((moments[n],), n)
        jet_point, jet_cov = jet_estimate_aligned({n: jet_group}, {n: jet_gap}, (n,))
        jet_independent_points[n] = jet_point[n]
        jet_independent_cov[n] = jet_cov

    scalar_points = {**scalar_source_points, **scalar_independent_points}
    jet_points = {**jet_source_points, **jet_independent_points}
    scalar_vector = [
        scalar_points[n][metric] for n in SIZE_ORDER for metric in METRICS
    ]
    jet_vector = [value for n in SIZE_ORDER for value in jet_points[n]]
    scalar_covariance = place_covariance(scalar_source_cov, scalar_independent_cov, 3)
    jet_covariance = place_covariance(jet_source_cov, jet_independent_cov, 5)

    scalar_labels = [f"N{start}_U" for start in (65, 85)]
    jet_labels = [f"N{start}_r{order}" for start in (65, 85) for order in ORDERS]
    scalar_u_indices = [index * 3 for index in range(len(SIZE_ORDER))]
    scalar_u_vector = [scalar_vector[index] for index in scalar_u_indices]
    scalar_u_covariance = [
        [scalar_covariance[i][j] for j in scalar_u_indices] for i in scalar_u_indices
    ]

    return {
        "schema": "matching-one/norm4-generation4-pilot-score/v1",
        "status": "prospective_pilot_score_of_post_reveal_frozen_recurrence",
        "prediction": str(prediction_path),
        "prediction_sha256": sha256(prediction_path),
        "size_order": list(SIZE_ORDER),
        "scalar_U_point": {str(n): scalar_points[n]["U"] for n in SIZE_ORDER},
        "thermal_jet_point": {str(n): jet_points[n] for n in SIZE_ORDER},
        "primary_lambda_half": {
            "scalar_U": transformed_score(
                scalar_u_vector, scalar_u_covariance, 1, 0.5, scalar_labels
            ),
            "thermal_jet": transformed_score(
                jet_vector, jet_covariance, 5, 0.5, jet_labels
            ),
        },
        "fixed_diagnostics_after_primary": {
            "lambda_0_pure_Jordan": {
                "scalar_U": transformed_score(
                    scalar_u_vector, scalar_u_covariance, 1, 0.0, scalar_labels
                ),
                "thermal_jet": transformed_score(
                    jet_vector, jet_covariance, 5, 0.0, jet_labels
                ),
            },
            "lambda_1_persistent_curvature": {
                "scalar_U": transformed_score(
                    scalar_u_vector, scalar_u_covariance, 1, 1.0, scalar_labels
                ),
                "thermal_jet": transformed_score(
                    jet_vector, jet_covariance, 5, 1.0, jet_labels
                ),
            },
        },
        "new_target_inputs": [
            {
                "N": n,
                "histogram": str(targets[n].histogram),
                "histogram_sha256": sha256(targets[n].histogram),
                "moments": str(targets[n].moments),
                "moments_sha256": sha256(targets[n].moments),
                "metadata": metadata[n],
            }
            for n in (520, 680)
        ],
        "evidence_guard": (
            "Scalar U and thermal jet reuse curves and are not additive. The "
            "lambda=0 and lambda=1 rows are fixed diagnostics after the frozen "
            "lambda=1/2 primary; no eigenvalue is fitted to this pilot."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--new-target", action="append", type=parse_target, required=True)
    parser.add_argument("--prediction", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render(args.new_target, args.prediction)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
