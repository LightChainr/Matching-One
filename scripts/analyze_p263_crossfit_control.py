#!/usr/bin/env python3
"""Audit revealed #263 controls and score the next-run exact bond control.

The Phase-E raw data contain event counts and J moments but no auxiliary with
a known finite-lattice expectation.  The revealed-data path therefore proves
that event-category conditioning is an algebraic no-op.  When a future stream
contains ``sum_b``, the same program applies the frozen even/odd two-fold bond
count control without changing the target expectation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

import mpmath as mp

from score_p263_boundary_qscore_pilot import GEOMETRY_ORDER, _gls


ANCHOR_INDEX = 1
ACTIVE_INDICES = (0, 2, 3)
CATEGORIES = ("1234", "12_34", "14_23", "other")


def read_rows(paths: Iterable[Path]) -> list[dict[str, int | str]]:
    rows: list[dict[str, int | str]] = []
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            for raw in csv.DictReader(handle):
                rows.append(
                    {
                        key: value if key == "geometry_id" else int(value)
                        for key, value in raw.items()
                    }
                )
    return rows


def _category_totals(rows: list[dict[str, int | str]]) -> dict[str, dict[str, int]]:
    output = {
        category: {"count": 0, "sum_J": 0} for category in CATEGORIES
    }
    for row in rows:
        samples = int(row["samples"])
        event_count = 0
        event_sum_j = 0
        for category in CATEGORIES[:-1]:
            count = int(row[f"count_{category}"])
            sum_j = int(row[f"sum_J_{category}"])
            output[category]["count"] += count
            output[category]["sum_J"] += sum_j
            event_count += count
            event_sum_j += sum_j
        output["other"]["count"] += samples - event_count
        output["other"]["sum_J"] += int(row["sum_J"]) - event_sum_j
    return output


def _direct_dlog(rows: list[dict[str, int | str]]) -> float:
    samples = sum(int(row["samples"]) for row in rows)
    count = sum(int(row["count_14_23"]) for row in rows)
    sum_j = sum(int(row["sum_J"]) for row in rows)
    sum_event_j = sum(int(row["sum_J_14_23"]) for row in rows)
    return sum_event_j / (2 * count) - sum_j / (2 * samples)


def conditional_category_noop(rows: list[dict[str, int | str]]) -> dict:
    """Cross-fit category means and verify exact reconstruction of dlog P."""

    by_geometry = {
        geometry: [row for row in rows if row["geometry_id"] == geometry]
        for geometry in GEOMETRY_ORDER
    }
    results = []
    for geometry, selected in by_geometry.items():
        folds = {
            parity: [row for row in selected if int(row["batch"]) % 2 == parity]
            for parity in (0, 1)
        }
        total_categories = _category_totals(selected)
        total_event_count = total_categories["14_23"]["count"]
        total_samples = sum(int(row["samples"]) for row in selected)
        conditional_reconstructed = 0.0
        overall_reconstructed = 0.0
        fold_details = []
        for heldout in (0, 1):
            training = 1 - heldout
            train_totals = _category_totals(folds[training])
            held_totals = _category_totals(folds[heldout])
            means = {}
            for category in CATEGORIES:
                if train_totals[category]["count"] == 0:
                    raise ValueError(
                        f"zero training count for {geometry} category {category}"
                    )
                means[category] = (
                    train_totals[category]["sum_J"]
                    / (2 * train_totals[category]["count"])
                )
            held_samples = sum(item["count"] for item in held_totals.values())
            held_event_count = held_totals["14_23"]["count"]
            conditional_residual = (
                held_totals["14_23"]["sum_J"] / (2 * held_event_count)
                - means["14_23"]
            )
            conditional_explained = means["14_23"]
            overall_residual = sum(
                item["sum_J"] / 2 - means[category] * item["count"]
                for category, item in held_totals.items()
            ) / held_samples
            overall_explained = sum(
                means[category] * item["count"]
                for category, item in held_totals.items()
            ) / held_samples
            conditional_reconstructed += (
                held_event_count / total_event_count
            ) * (conditional_residual + conditional_explained)
            overall_reconstructed += (held_samples / total_samples) * (
                overall_residual + overall_explained
            )
            fold_details.append(
                {
                    "heldout_batch_parity": heldout,
                    "training_batch_parity": training,
                    "heldout_14_23_count": held_event_count,
                    "conditional_residual": conditional_residual,
                    "conditional_explained": conditional_explained,
                    "overall_residual": overall_residual,
                    "overall_explained": overall_explained,
                }
            )
        reconstructed = conditional_reconstructed - overall_reconstructed
        direct = _direct_dlog(selected)
        results.append(
            {
                "geometry_id": geometry,
                "direct_d_log_probability": direct,
                "crossfit_conditional_reconstruction": reconstructed,
                "absolute_difference": abs(reconstructed - direct),
                "folds": fold_details,
            }
        )
    return {
        "identity": (
            "E[J/2|H]-E[J/2] = "
            "{E[J/2-m(X)|H]-E[J/2-m(X)]} + "
            "{m(H)-E[m(X)]} for any category function m"
        ),
        "geometry_results": results,
        "maximum_absolute_difference": max(
            result["absolute_difference"] for result in results
        ),
        "variance_ratio_to_primary": 1.0,
        "reason_no_reduction": (
            "The residual and explained terms use the same held-out event "
            "counts and sum exactly to the original estimator in every fold."
        ),
    }


def _mean(values: list[list[float]]) -> list[float]:
    return [
        sum(row[index] for row in values) / len(values)
        for index in range(len(values[0]))
    ]


def _mean_covariance(values: list[list[float]]) -> list[list[float]]:
    center = _mean(values)
    count = len(values)
    return [
        [
            sum(
                (row[first] - center[first]) * (row[second] - center[second])
                for row in values
            )
            / (count * (count - 1))
            for second in range(len(center))
        ]
        for first in range(len(center))
    ]


def _pseudoinverse(matrix: list[list[float]]) -> list[list[float]]:
    work = mp.matrix(matrix)
    values, vectors = mp.eigsy(work)
    maximum = max(abs(values[index]) for index in range(len(matrix)))
    cutoff = maximum * mp.mpf("1e-10")
    inverse = mp.zeros(len(matrix))
    for index in range(len(matrix)):
        if values[index] > cutoff:
            column = vectors[:, index]
            inverse += (column * column.T) / values[index]
    return [
        [float(inverse[row, column]) for column in range(len(matrix))]
        for row in range(len(matrix))
    ]


def _fit_beta(outcomes: list[list[float]], controls: list[list[float]]) -> list[list[float]]:
    mean_y = _mean(outcomes)
    mean_w = _mean(controls)
    count = len(outcomes)
    covariance_ww = [
        [
            sum(
                (row[first] - mean_w[first]) * (row[second] - mean_w[second])
                for row in controls
            )
            / (count - 1)
            for second in range(len(mean_w))
        ]
        for first in range(len(mean_w))
    ]
    covariance_yw = [
        [
            sum(
                (outcomes[row][target] - mean_y[target])
                * (controls[row][control] - mean_w[control])
                for row in range(count)
            )
            / (count - 1)
            for control in range(len(mean_w))
        ]
        for target in range(len(mean_y))
    ]
    inverse = _pseudoinverse(covariance_ww)
    return [
        [
            sum(covariance_yw[target][middle] * inverse[middle][control]
                for middle in range(len(mean_w)))
            for control in range(len(mean_w))
        ]
        for target in range(len(mean_y))
    ]


def _matrix_vector(matrix: list[list[float]], vector: list[float]) -> list[float]:
    return [
        sum(coefficient * value for coefficient, value in zip(row, vector))
        for row in matrix
    ]


def crossfit_bond_control(rows: list[dict[str, int | str]], score_payload: dict) -> dict:
    """Apply the exact-zero open-bond control to a next-run batch stream."""

    if not rows or "sum_b" not in rows[0]:
        raise ValueError("sum_b is required for the exact bond control")
    batch_sets = {
        geometry: {
            int(row["batch"]) for row in rows if row["geometry_id"] == geometry
        }
        for geometry in GEOMETRY_ORDER
    }
    if any(batch_sets[geometry] != batch_sets[GEOMETRY_ORDER[0]] for geometry in GEOMETRY_ORDER):
        raise ValueError("geometries do not share synchronized batch ids")
    batch_ids = sorted(batch_sets[GEOMETRY_ORDER[0]])
    row_map = {
        (str(row["geometry_id"]), int(row["batch"])): row for row in rows
    }
    batch_count = len(batch_ids)

    total_samples = {
        geometry: sum(
            int(row_map[(geometry, batch)]["samples"]) for batch in batch_ids
        )
        for geometry in GEOMETRY_ORDER
    }
    total_events = {
        geometry: sum(
            int(row_map[(geometry, batch)]["count_14_23"]) for batch in batch_ids
        )
        for geometry in GEOMETRY_ORDER
    }

    def controls_for(batch: int) -> list[float]:
        values = []
        for geometry in GEOMETRY_ORDER:
            row = row_map[(geometry, batch)]
            samples = int(row["samples"])
            edges = int(row["edges"])
            centered_twice = 2 * int(row["sum_b"]) - samples * edges
            values.append(centered_twice / math.sqrt(samples * edges))
        return values

    def full_geometry_contribution(geometry: str, batch: int) -> float:
        row = row_map[(geometry, batch)]
        return batch_count * (
            int(row["sum_J_14_23"]) / (2 * total_events[geometry])
            - int(row["sum_J"]) / (2 * total_samples[geometry])
        )

    dlog = [
        sum(full_geometry_contribution(geometry, batch) for batch in batch_ids)
        / batch_count
        for geometry in GEOMETRY_ORDER
    ]
    primary_residual = [float(value) for value in score_payload["residual"]]
    deterministic = [
        primary_residual[position] - (dlog[index] - dlog[ANCHOR_INDEX])
        for position, index in enumerate(ACTIVE_INDICES)
    ]
    raw_contributions = []
    controls = []
    for batch in batch_ids:
        geometry_values = [
            full_geometry_contribution(geometry, batch)
            for geometry in GEOMETRY_ORDER
        ]
        raw_contributions.append(
            [
                geometry_values[index]
                - geometry_values[ANCHOR_INDEX]
                + deterministic[position]
                for position, index in enumerate(ACTIVE_INDICES)
            ]
        )
        controls.append(controls_for(batch))

    def training_outcomes(training_batches: list[int]) -> list[list[float]]:
        training_samples = {
            geometry: sum(
                int(row_map[(geometry, batch)]["samples"])
                for batch in training_batches
            )
            for geometry in GEOMETRY_ORDER
        }
        training_events = {
            geometry: sum(
                int(row_map[(geometry, batch)]["count_14_23"])
                for batch in training_batches
            )
            for geometry in GEOMETRY_ORDER
        }
        fold_count = len(training_batches)
        output = []
        for batch in training_batches:
            values = []
            for geometry in GEOMETRY_ORDER:
                row = row_map[(geometry, batch)]
                values.append(
                    fold_count
                    * (
                        int(row["sum_J_14_23"])
                        / (2 * training_events[geometry])
                        - int(row["sum_J"]) / (2 * training_samples[geometry])
                    )
                )
            output.append(
                [values[index] - values[ANCHOR_INDEX] for index in ACTIVE_INDICES]
            )
        return output

    fold_positions = {
        parity: [position for position, batch in enumerate(batch_ids) if batch % 2 == parity]
        for parity in (0, 1)
    }
    beta = {}
    for parity in (0, 1):
        positions = fold_positions[parity]
        training_batches = [batch_ids[position] for position in positions]
        beta[parity] = _fit_beta(
            training_outcomes(training_batches),
            [controls[position] for position in positions],
        )

    adjusted = []
    for position, batch in enumerate(batch_ids):
        training_parity = 1 - (batch % 2)
        correction = _matrix_vector(beta[training_parity], controls[position])
        adjusted.append(
            [
                raw_contributions[position][index] - correction[index]
                for index in range(len(ACTIVE_INDICES))
            ]
        )

    raw_mean = _mean(raw_contributions)
    adjusted_mean = _mean(adjusted)
    raw_covariance = _mean_covariance(raw_contributions)
    adjusted_covariance = _mean_covariance(adjusted)
    raw_gls = _gls(raw_mean, raw_covariance)
    adjusted_gls = _gls(adjusted_mean, adjusted_covariance)
    raw_trace = sum(raw_covariance[index][index] for index in range(3))
    adjusted_trace = sum(adjusted_covariance[index][index] for index in range(3))
    return {
        "fold_rule": "even versus odd synchronized batch ids",
        "control_order": list(GEOMETRY_ORDER),
        "control": "(2*sum_b-samples*edges)/sqrt(samples*edges)",
        "control_expectation": 0,
        "target_preservation": (
            "Each beta is trained on the opposite fold; conditional on beta, "
            "the held-out correction has exact expectation zero at p=1/2."
        ),
        "beta_trained_on_even": beta[0],
        "beta_trained_on_odd": beta[1],
        "raw_batch_linearized": {
            "residual": raw_mean,
            "covariance": raw_covariance,
            "joint_gls": raw_gls,
        },
        "crossfit_bond_control": {
            "residual": adjusted_mean,
            "covariance": adjusted_covariance,
            "joint_gls": adjusted_gls,
        },
        "variance_comparison": {
            "diagonal_ratio_control_over_raw": [
                adjusted_covariance[index][index] / raw_covariance[index][index]
                for index in range(3)
            ],
            "trace_ratio_control_over_raw": adjusted_trace / raw_trace,
        },
    }


def analyze_level(rows: list[dict[str, int | str]], score_payload: dict) -> dict:
    noop = conditional_category_noop(rows)
    has_bond_control = bool(rows) and "sum_b" in rows[0]
    result = {
        "revealed_conditional_score": noop,
        "frozen_primary_unchanged": {
            "residual": score_payload["residual"],
            "residual_covariance": score_payload["residual_covariance"],
            "joint_gls": score_payload["joint_gls"],
        },
        "strict_control_available": has_bond_control,
    }
    if has_bond_control:
        result["strict_crossfit_bond_control"] = crossfit_bond_control(
            rows, score_payload
        )
    else:
        result["strict_crossfit_bond_control"] = {
            "status": "not_identifiable_from_revealed_sufficient_statistics",
            "missing_minimal_field": "sum_b per synchronized geometry/batch",
            "missing_training_cross_matrix": (
                "Sigma_YW cannot be formed because the exact-zero bond-count "
                "control W is absent."
            ),
        }
    return result


def render(
    level1_batches: Path,
    level1_score: Path,
    level2_batches: Path,
    level2_score: Path,
) -> dict:
    levels = {}
    for name, batches, score in (
        ("level1_200k", level1_batches, level1_score),
        ("level2_500k", level2_batches, level2_score),
    ):
        levels[name] = analyze_level(
            read_rows([batches]), json.loads(score.read_text(encoding="utf-8"))
        )
    return {
        "schema": "matching-one.p263-crossfit-control-audit.v1",
        "issue": 263,
        "status": "revealed_data_no_go_and_next_run_schema",
        "source_commit": "90da884",
        "levels": levels,
        "strict_no_go": {
            "available_auxiliaries": [
                "J and J^2 moments with unknown finite-lattice expectations",
                "three mutually exclusive event counts with unknown finite-lattice probabilities",
            ],
            "why_opposite_fold_centering_fails": (
                "For a fixed beta, the two equal-fold corrections "
                "beta*(X_A-X_B) and beta*(X_B-X_A) cancel exactly. Allowing "
                "fold-specific betas leaves a noisy beta-difference times "
                "mean-difference product, not an exact known-mean control."
            ),
            "minimal_missing_field": "sum_b",
            "minimal_exact_control": (
                "W_g=(2*sum_b_g-n_g*edges_g)/sqrt(n_g*edges_g), E[W_g]=0"
            ),
            "batch_crossfit_extra_cross_moments_required": [],
            "sample_level_analytic_beta_extra_fields": [
                "for every geometry pair (g,h): sum J_g*b_h",
                "for every geometry pair (g,h) and pattern p: sum I_p,g*J_g*b_h",
            ],
            "sample_level_control_covariance": (
                "Sigma_WW is already exact from edge counts and the shared-edge "
                "CRN enumeration; b_g*b_h is only an optional audit moment."
            ),
        },
        "next_run_contract": (
            "experiments/p263_boundary_qscore_control_phaseF_20260829.yaml"
        ),
        "claim_boundary": [
            "No new variance-reduced chi-square is claimed from the revealed Phase-E raw data.",
            "The unchanged primary chi-squares are the exact result of the available conditional decomposition.",
            "The bond-control scorer becomes executable only on a new stream containing sum_b.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level1-batches", type=Path, required=True)
    parser.add_argument("--level1-score", type=Path, required=True)
    parser.add_argument("--level2-batches", type=Path, required=True)
    parser.add_argument("--level2-score", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render(
        args.level1_batches,
        args.level1_score,
        args.level2_batches,
        args.level2_score,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
