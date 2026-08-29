#!/usr/bin/env python3
"""Post-reveal context/morphism map from the existing P200 N650 calibration."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp

from score_p200_n650_mixed_join import (
    PRIMARY,
    delete_one_summary,
    load_inputs,
    render as frozen_score,
)


ROOT = Path(__file__).resolve().parents[1]


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repository_path(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def transform(
    mean: list[float], covariance: list[list[float]], matrix: list[list[float]]
) -> tuple[list[float], list[list[float]]]:
    output_mean = [sum(row[index] * mean[index] for index in range(4)) for row in matrix]
    output_covariance = [
        [
            sum(
                matrix[first][i] * covariance[i][j] * matrix[second][j]
                for i in range(4) for j in range(4)
            )
            for second in range(len(matrix))
        ]
        for first in range(len(matrix))
    ]
    return output_mean, output_covariance


def zero_test(names: list[str], mean: list[float], covariance: list[list[float]]) -> dict:
    mp.mp.dps = 80
    values, vectors = mp.eigsy(mp.matrix(covariance))
    maximum = max(abs(values[index]) for index in range(len(mean)))
    cutoff = max(mp.mpf("1e-50"), maximum * mp.mpf("1e-10"))
    inverse = mp.zeros(len(mean))
    rank = 0
    for index in range(len(mean)):
        if values[index] > cutoff:
            column = vectors[:, index]
            inverse += column * column.T / values[index]
            rank += 1
    vector = mp.matrix(mean)
    chi_square = (vector.T * inverse * vector)[0]
    p_value = (
        mp.gammainc(mp.mpf(rank) / 2, chi_square / 2, mp.inf, regularized=True)
        if rank else mp.nan
    )
    marginal = []
    for index, name in enumerate(names):
        se = math.sqrt(covariance[index][index])
        marginal.append({"name": name, "mean": mean[index], "se": se, "z": mean[index] / se})
    return {
        "state_order": names,
        "mean": mean,
        "covariance": covariance,
        "chi_square": float(chi_square),
        "degrees_of_freedom": rank,
        "p_value": float(p_value),
        "log10_p_value": float(mp.log10(p_value)),
        "marginal": marginal,
    }


def _group_mean(rows: list[dict], batches: set[int], indices: tuple[int, ...]) -> list[float]:
    selected = [row for row in rows if row["batch"] in batches]
    samples = sum(row["samples"] for row in selected)
    return [
        sum(row["primary"][index] for row in selected) / (2 * samples)
        for index in indices
    ]


def analyze(
    batch_path: Path,
    metadata_path: Path,
    prediction_path: Path,
    committed_score_path: Path,
) -> dict:
    computed = frozen_score(batch_path, metadata_path, prediction_path)
    committed = json.loads(committed_score_path.read_text(encoding="utf-8"))
    for section in ("primary", "secondary_ambient_H1"):
        if computed[section]["mean"] != committed[section]["mean"]:
            raise ValueError(f"committed {section} mean does not recompute")
        if computed[section]["delete_one_covariance"] != committed[section]["delete_one_covariance"]:
            raise ValueError(f"committed {section} covariance does not recompute")

    metadata, _, rows = load_inputs(batch_path, metadata_path, prediction_path)
    prediction = json.loads(prediction_path.read_text(encoding="utf-8"))
    source_fibers = int(prediction["acquisition"]["source_fibers_per_orientation"])
    if source_fibers != 65 or metadata["samples"] != 20000:
        raise ValueError("this map is frozen to the revealed N650 20k calibration")

    primary_mean = computed["primary"]["mean"]
    primary_covariance = computed["primary"]["delete_one_covariance"]
    ambient_mean = computed["secondary_ambient_H1"]["mean"]
    ambient_covariance = computed["secondary_ambient_H1"]["delete_one_covariance"]

    s_matrix = [[1, 0, 0, 0], [0, 0, 1, 0]]
    d_matrix = [[0, 1, 0, 0], [0, 0, 0, 1]]
    orientation_color_matrix = [
        [0.5, 0.5, 0.5, 0.5],
        [0.5, 0.5, -0.5, -0.5],
        [0.5, -0.5, 0.5, -0.5],
        [0.5, -0.5, -0.5, 0.5],
    ]
    common_color_matrix = [
        [0.5, 0, 0.5, 0],
        [0.5, 0, -0.5, 0],
    ]
    cross_geometry_matrix = [[0, -2, 0, 0], [0, 0, 0, -2]]

    primary_s = transform(primary_mean, primary_covariance, s_matrix)
    primary_d = transform(primary_mean, primary_covariance, d_matrix)
    primary_by_context = transform(
        primary_mean, primary_covariance, orientation_color_matrix
    )
    common_color = transform(primary_mean, primary_covariance, common_color_matrix)
    cross_geometry = transform(primary_mean, primary_covariance, cross_geometry_matrix)
    ambient_s = transform(ambient_mean, ambient_covariance, s_matrix)
    ambient_d = transform(ambient_mean, ambient_covariance, d_matrix)
    ambient_by_context = transform(
        ambient_mean, ambient_covariance, orientation_color_matrix
    )

    local = prediction["toy_exact_normalization"]["p_ref"]
    local_black_density = float(local["mean_decimal"][0])
    local_white_density = float(local["mean_decimal"][1])
    local_covariance = [[float(value) for value in row] for row in local["covariance_decimal"]]
    local_even_variance = local_covariance[3][3]
    local_odd_variance = local_covariance[2][2]
    local_even_scale = math.sqrt(source_fibers * local_even_variance)
    local_odd_scale = math.sqrt(source_fibers * local_odd_variance)

    color_names = ["first_black", "first_white", "second_black", "second_white"]
    color_rows = []
    for index, name in enumerate(color_names):
        residual = primary_by_context[0][index]
        se = math.sqrt(primary_by_context[1][index][index])
        baseline_density = local_black_density if "black" in name else local_white_density
        baseline = source_fibers * baseline_density
        color_rows.append({
            "context": name,
            "exact_local_incidence_baseline_mean": baseline,
            "connected_residual_mean": residual,
            "connected_residual_se": se,
            "implied_full_mixed_join_mean": baseline + residual,
            "per_source_fiber": {
                "local_incidence_baseline": baseline_density,
                "connected_residual": residual / source_fibers,
                "implied_full_mixed_join": baseline_density + residual / source_fibers,
            },
        })

    common_density = [value / source_fibers for value in common_color[0]]
    common_density_se = [
        math.sqrt(common_color[1][index][index]) / source_fibers
        for index in range(2)
    ]
    target_source_fibers = 130
    target_black = target_source_fibers * common_density[0]
    target_white = target_source_fibers * common_density[1]

    _, _, primary_leave_one = delete_one_summary(rows, "primary")
    sign_stability = {
        PRIMARY[index]: sum(
            mp.sign(values[index]) == mp.sign(primary_mean[index])
            for values in primary_leave_one
        ) / len(primary_leave_one)
        for index in (0, 2)
    }
    batch_ids = {row["batch"] for row in rows}
    even_batches = {batch for batch in batch_ids if batch % 2 == 0}
    odd_batches = batch_ids.difference(even_batches)

    card = [
        "Question: Does the revealed N650 signal require ordered path/state memory, or only a static mixed-factor interaction?",
        "Exact boundary: the runner subtracts local C2xC5 incidence cycles configurationwise, and the two final joins commute.",
        "Observed: residual common modes ES and OS survive overwhelmingly, while geometry-difference ED/OD jointly remain compatible with zero.",
        "Not identified: the raw batches contain no intermediate h0,h2,h5,h25 state or chronology, so ES/OS cannot diagnose path memory.",
        "Freeze: transport rho_B=-0.67716808 and rho_W=-0.14910692 per source fiber only as a future geometry/scale challenge, with scale discrepancy unestimated.",
    ]

    return {
        "schema": "matching-one.p200-n650-context-morphism-map.v1",
        "issue": 200,
        "status": "post_reveal_opportunity_map_no_new_production",
        "input": {
            "batch_csv": repository_path(batch_path),
            "metadata": repository_path(metadata_path),
            "committed_score": repository_path(committed_score_path),
            "samples": metadata["samples"],
            "batches": metadata["batches"],
            "N": 650,
            "source_fibers": source_fibers,
            "sha256": {
                "batch_csv": file_sha256(batch_path),
                "metadata": file_sha256(metadata_path),
                "committed_score": file_sha256(committed_score_path),
                "prediction": file_sha256(prediction_path),
            },
        },
        "identifiability_boundary": {
            "mixed_factor_interaction": (
                "identified by nonzero R_c=J_full-J_local after the exact local incidence subtraction"
            ),
            "geometry_context": (
                "ED/OD compare two static HNF embeddings under shared counters; they are not join-order variables"
            ),
            "path_or_state_memory": (
                "not identified because batches omit intermediate corner states and chronological filtrations"
            ),
            "typed_ambient_H1": (
                "endpoint mixed defect only, with frozen representative-displacement convention; correlated secondary"
            ),
        },
        "local_incidence_subtraction": {
            "performed_in_runner": True,
            "raw_contains_only_residual": True,
            "formula": "R_c=J_full,c-sum_over_65_fibers b1(Inc_c)",
            "exact_p_ref_baseline_per_source_fiber": {
                "black": local_black_density,
                "white": local_white_density,
            },
            "context_reconstruction": color_rows,
        },
        "primary_partition_residual": {
            "state_order": list(PRIMARY),
            "mean": primary_mean,
            "covariance": primary_covariance,
            "common_geometry_S": zero_test(["ES", "OS"], *primary_s),
            "geometry_difference_D": zero_test(["ED", "OD"], *primary_d),
            "orientation_contexts": {
                "state_order": color_names,
                "mean": primary_by_context[0],
                "covariance": primary_by_context[1],
            },
            "second_minus_first_geometry": zero_test(
                ["E_second-E_first", "O_second-O_first"], *cross_geometry
            ),
            "one_orientation_exact_local_reference_scale": {
                "even": local_even_scale,
                "odd": local_odd_scale,
                "scaled_state_order": list(PRIMARY),
                "scaled_mean": [
                    primary_mean[0] / local_even_scale,
                    primary_mean[1] / local_even_scale,
                    primary_mean[2] / local_odd_scale,
                    primary_mean[3] / local_odd_scale,
                ],
            },
            "per_source_fiber_state_mean": [value / source_fibers for value in primary_mean],
            "leave_one_sign_stability": sign_stability,
            "batch_parity_sensitivity_ES_OS": {
                "even_batches": _group_mean(rows, even_batches, (0, 2)),
                "odd_batches": _group_mean(rows, odd_batches, (0, 2)),
                "role": "post_reveal robustness only",
            },
        },
        "typed_ambient_H1": {
            "scaling": "global rank O(1), not divided by the 65 source fibers",
            "state_order": ["ambient_ES", "ambient_ED", "ambient_OS", "ambient_OD"],
            "mean": ambient_mean,
            "covariance": ambient_covariance,
            "common_geometry_S": zero_test(["ambient_ES", "ambient_OS"], *ambient_s),
            "geometry_difference_D": zero_test(["ambient_ED", "ambient_OD"], *ambient_d),
            "orientation_color_contexts": {
                "state_order": color_names,
                "mean": ambient_by_context[0],
                "covariance": ambient_by_context[1],
            },
            "interpretation": (
                "the surviving common even rank defect is convention-labelled endpoint topology, not chronology"
            ),
        },
        "morphism_parameter_freeze": {
            "classification": "post_reveal_conditional_future_challenge",
            "common_connected_residual_density_per_source_fiber": {
                "black": {"estimate": common_density[0], "se": common_density_se[0]},
                "white": {"estimate": common_density[1], "se": common_density_se[1]},
            },
            "same_N650_new_HNF_geometry_prediction": {
                "E": primary_mean[0],
                "O": primary_mean[2],
                "fit_covariance_EO": primary_s[1],
                "condition": "same C2xC5 observable, p_ref, source-fiber count, and graph typing",
                "model_discrepancy": "not estimated; ED/OD are only a two-geometry compatibility check",
            },
            "conditional_N1300_prediction": {
                "source_fibers": target_source_fibers,
                "unrun_geometry_pair": [
                    {
                        "name": "source_11_plus_3i",
                        "source_gaussian": [11, 3],
                        "N260_after_1_plus_i": [8, 14],
                        "N650_after_2_minus_i": [25, -5],
                        "N1300_after_3_plus_i": [30, 20],
                        "final_period_matrix": [[30, -20], [20, 30]],
                    },
                    {
                        "name": "source_7_plus_9i",
                        "source_gaussian": [7, 9],
                        "N260_after_1_plus_i": [-2, 16],
                        "N650_after_2_minus_i": [23, 11],
                        "N1300_after_3_plus_i": [12, 34],
                        "final_period_matrix": [[12, -34], [34, 12]],
                    },
                ],
                "state_order": list(PRIMARY),
                "mean": [target_black + target_white, 0.0, target_black - target_white, 0.0],
                "black_connected_residual": {
                    "mean": target_black,
                    "parameter_fit_se_only": target_source_fibers * common_density_se[0],
                },
                "white_connected_residual": {
                    "mean": target_white,
                    "parameter_fit_se_only": target_source_fibers * common_density_se[1],
                },
                "assumptions": [
                    "connected residual is extensive in the number of C2xC5 source fibers",
                    "geometry-difference ED and OD remain zero",
                    "p_ref and typed graph semantics are unchanged",
                ],
                "model_discrepancy": "unestimated from a single N; this is a freeze, not an established scale law",
            },
        },
        "scientific_card": card,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--prediction", type=Path, required=True)
    parser.add_argument("--score", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--card", type=Path, required=True)
    args = parser.parse_args()
    payload = analyze(args.batches, args.metadata, args.prediction, args.score)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.card.write_text("\n".join(payload["scientific_card"]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
