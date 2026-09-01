#!/usr/bin/env python3
"""Existing-data forward-identifiability audit for Issue #275.

This script consumes the complete aligned-delete-one covariance from the
canonical norm-4 K1/K2 decomposition.  It generates no Monte Carlo samples.
The two primitive coordinates are the direction-normalized activation
responses

    d_i(N) = Delta_4 F_i(p_bar),  i=1,2,

and the scaled coordinates are u_i(N)=N^(13/8)d_i(N).  On each dyadic
lineage N -> 2N -> 4N, the audit compares

    semisimple: u_i(N)=a_i+b_i N^(-delta), kappa=2^(-delta),
    Jordan:     u_i(N)=a_i+c_i log_2(N/N0).

All lineage/activation amplitudes remain free.  Only the common transfer
coordinate kappa is profiled.  The output records the exact design ranks and
the non-uniform kappa -> 1 collision that prevents a uniform model-selection
claim without a separately normalized restricted-trace transport relation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import mpmath as mp
import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "analysis" / "p275_forward_identifiability_manifest.yaml"
OUTPUT_SCHEMA = "matching-one.p275-forward-identifiability.v1"
MANIFEST_SCHEMA = "matching-one.p275-forward-identifiability.manifest.v1"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite_float(value: float) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("non-finite output value")
    return value


def vector_payload(vector: np.ndarray) -> list[float]:
    return [finite_float(value) for value in vector]


def matrix_payload(matrix: np.ndarray) -> list[list[float]]:
    return [vector_payload(row) for row in matrix]


def load_manifest(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("unexpected P275 forward-identifiability manifest")
    return payload


def resolve_repo_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def chi_square_survival(value: float, degrees: int) -> float:
    if value < 0 or degrees <= 0:
        raise ValueError("invalid chi-square request")
    mp.mp.dps = max(mp.mp.dps, 60)
    shape = mp.mpf(degrees) / 2
    probability = mp.gammainc(shape, mp.mpf(str(value)) / 2, mp.inf) / mp.gamma(shape)
    return float(probability)


def covariance_inverse(
    covariance: np.ndarray, relative_cutoff: float, absolute_floor: float
) -> tuple[np.ndarray, dict[str, Any]]:
    symmetric = (covariance + covariance.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    scale = float(max(np.max(np.abs(eigenvalues)), 0.0))
    cutoff = max(scale * relative_cutoff, absolute_floor)
    if float(np.min(eigenvalues)) < -cutoff:
        raise ValueError("covariance has a negative eigenvalue beyond cutoff")
    kept = eigenvalues > cutoff
    if not np.any(kept):
        raise ValueError("covariance has zero numerical rank")
    inverse = (eigenvectors[:, kept] / eigenvalues[kept]) @ eigenvectors[:, kept].T
    return inverse, {
        "eigenvalues_ascending": vector_payload(eigenvalues),
        "relative_cutoff": relative_cutoff,
        "absolute_floor": absolute_floor,
        "applied_cutoff": cutoff,
        "rank": int(np.count_nonzero(kept)),
        "discarded_modes": int(len(eigenvalues) - np.count_nonzero(kept)),
    }


def validate_input(payload: Mapping[str, Any], manifest: Mapping[str, Any], input_path: Path) -> None:
    contract = manifest["input"]
    if payload.get("schema") != contract["schema"]:
        raise ValueError("input schema mismatch")
    if sha256_file(input_path) != contract["sha256"]:
        raise ValueError("input SHA256 mismatch")
    if payload.get("provenance", {}).get("manifest_sha256") != contract["source_manifest_sha256"]:
        raise ValueError("source two-activation manifest hash mismatch")

    required_sizes = sorted({int(n) for chain in manifest["lineages"] for n in chain["sizes"]})
    if sorted(int(n) for n in payload.get("size_order", [])) != required_sizes:
        raise ValueError("input size set differs from frozen lineages")

    decision = payload.get("decision_covariance", {})
    order = decision.get("metric_order_with_N")
    estimate = decision.get("estimate_vector")
    covariance = decision.get("jackknife_covariance")
    if not isinstance(order, list) or not isinstance(estimate, list) or not isinstance(covariance, list):
        raise ValueError("input decision covariance is incomplete")
    if len(order) != len(estimate) or len(covariance) != len(order):
        raise ValueError("input decision dimensions disagree")
    if any(not isinstance(row, list) or len(row) != len(order) for row in covariance):
        raise ValueError("input covariance is not square")

    keys = [(int(row["N"]), str(row["metric"])) for row in order]
    if len(set(keys)) != len(keys):
        raise ValueError("duplicate input decision coordinate")
    for n in required_sizes:
        for metric in manifest["observable"]["metrics"]:
            if (n, metric) not in set(keys):
                raise ValueError(f"missing N={n} metric={metric}")


def select_scaled_activation_block(
    payload: Mapping[str, Any], manifest: Mapping[str, Any]
) -> dict[str, Any]:
    decision = payload["decision_covariance"]
    full_order = decision["metric_order_with_N"]
    full_estimate = np.asarray(decision["estimate_vector"], dtype=float)
    full_covariance = np.asarray(decision["jackknife_covariance"], dtype=float)
    index = {
        (int(row["N"]), str(row["metric"])): position
        for position, row in enumerate(full_order)
    }
    exponent = float(manifest["scaling"]["leading_N_exponent"])

    selected_indices: list[int] = []
    selected_order: list[dict[str, Any]] = []
    for lineage in manifest["lineages"]:
        sizes = [int(n) for n in lineage["sizes"]]
        for activation, metric in enumerate(manifest["observable"]["metrics"], start=1):
            for generation, n in enumerate(sizes):
                selected_indices.append(index[(n, metric)])
                selected_order.append({
                    "lineage": str(lineage["id"]),
                    "activation": f"F{activation}",
                    "generation": generation,
                    "N": n,
                    "metric": metric,
                })

    raw_estimate = full_estimate[selected_indices]
    raw_covariance = full_covariance[np.ix_(selected_indices, selected_indices)]
    scale = np.asarray([float(row["N"]) ** exponent for row in selected_order])
    scaled_estimate = scale * raw_estimate
    scaled_covariance = scale[:, None] * raw_covariance * scale[None, :]
    return {
        "selected_order": selected_order,
        "raw_estimate": raw_estimate,
        "raw_covariance": raw_covariance,
        "scale": scale,
        "scaled_estimate": scaled_estimate,
        "scaled_covariance": scaled_covariance,
    }


def semisimple_block(kappa: float) -> np.ndarray:
    return np.asarray([[1.0, 1.0], [1.0, kappa], [1.0, kappa * kappa]])


def jordan_block() -> np.ndarray:
    return np.asarray([[1.0, 0.0], [1.0, 1.0], [1.0, 2.0]])


def residual_operator(block_count: int, kappa: float) -> np.ndarray:
    row = np.asarray([kappa, -(1.0 + kappa), 1.0])
    operator = np.zeros((block_count, 3 * block_count))
    for block in range(block_count):
        operator[block, 3 * block : 3 * block + 3] = row
    return operator


def score_at_kappa(
    kappa: float,
    estimate: np.ndarray,
    covariance: np.ndarray,
    numerics: Mapping[str, Any],
    fitted_parameters: int,
    residual_order: Sequence[str],
) -> dict[str, Any]:
    chi_square, residual, residual_covariance, spectrum = chi_square_at_kappa(
        kappa, estimate, covariance, numerics
    )
    degrees = int(spectrum["rank"] - fitted_parameters)
    if degrees <= 0:
        raise ValueError("profile consumes every residual covariance mode")
    p_value = chi_square_survival(chi_square, degrees)
    standard_errors = np.sqrt(np.maximum(np.diag(residual_covariance), 0.0))
    alpha = float(numerics["decision_alpha"])
    return {
        "kappa": finite_float(kappa),
        "residual_order": list(residual_order),
        "residual": vector_payload(residual),
        "residual_standard_error": vector_payload(standard_errors),
        "residual_z": vector_payload(residual / standard_errors),
        "residual_covariance": matrix_payload(residual_covariance),
        "covariance_spectrum": spectrum,
        "mahalanobis_chi_square": finite_float(chi_square),
        "degrees_of_freedom": degrees,
        "chi_square_survival_p": finite_float(p_value),
        "decision_alpha": alpha,
        "excluded_at_alpha": bool(p_value < alpha),
    }


def chi_square_at_kappa(
    kappa: float,
    estimate: np.ndarray,
    covariance: np.ndarray,
    numerics: Mapping[str, Any],
) -> tuple[float, np.ndarray, np.ndarray, dict[str, Any]]:
    """Return the GLS residual score without evaluating a survival function.

    This is the hot path for the kappa profiles.  Exact p-values and the full
    diagnostic payload are evaluated only once, at the selected candidates.
    """
    if len(estimate) % 3:
        raise ValueError("scaled vector does not contain complete three-generation blocks")
    operator = residual_operator(len(estimate) // 3, kappa)
    residual = operator @ estimate
    residual_covariance = operator @ covariance @ operator.T
    inverse, spectrum = covariance_inverse(
        residual_covariance,
        float(numerics["eigen_relative_cutoff"]),
        float(numerics["eigen_absolute_floor"]),
    )
    chi_square = float(residual @ inverse @ residual)
    return chi_square, residual, residual_covariance, spectrum


def golden_minimum(
    function: Callable[[float], float], lower: float, upper: float, tolerance: float
) -> tuple[float, float]:
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    left = upper - ratio * (upper - lower)
    right = lower + ratio * (upper - lower)
    f_left, f_right = function(left), function(right)
    while upper - lower > tolerance:
        if f_left <= f_right:
            upper, right, f_right = right, left, f_left
            left = upper - ratio * (upper - lower)
            f_left = function(left)
        else:
            lower, left, f_left = left, right, f_right
            right = lower + ratio * (upper - lower)
            f_right = function(right)
    point = (lower + upper) / 2.0
    return point, function(point)


def bounded_profile(
    function: Callable[[float], float], lower: float, upper: float,
    grid_points: int, tolerance: float
) -> dict[str, Any]:
    grid = np.linspace(lower, upper, grid_points)
    values = np.asarray([function(float(value)) for value in grid])
    candidates: list[tuple[float, float]] = [
        (float(grid[0]), float(values[0])),
        (float(grid[-1]), float(values[-1])),
    ]
    local_indices = [
        index for index in range(1, len(grid) - 1)
        if values[index] <= values[index - 1] and values[index] <= values[index + 1]
    ]
    for index in local_indices:
        candidates.append(golden_minimum(
            function, float(grid[index - 1]), float(grid[index + 1]), tolerance
        ))
    point, value = min(candidates, key=lambda item: item[1])
    boundary_tolerance = max(tolerance, (upper - lower) / max(grid_points - 1, 1))
    return {
        "support": [lower, upper],
        "best_kappa": finite_float(point),
        "best_chi_square": finite_float(value),
        "best_is_lower_boundary": bool(abs(point - lower) <= boundary_tolerance),
        "best_is_upper_boundary": bool(abs(point - upper) <= boundary_tolerance),
        "grid_points": grid_points,
        "local_minima_refined": len(local_indices),
        "refinement_tolerance": tolerance,
    }


def design_audit(manifest: Mapping[str, Any]) -> dict[str, Any]:
    blocks = 2 * len(manifest["lineages"])
    kappa = float(manifest["models"]["fixed_semisimple_q2"]["kappa"])
    semisimple = semisimple_block(kappa)
    jordan = jordan_block()
    full_semisimple = np.kron(np.eye(blocks), semisimple)
    full_jordan = np.kron(np.eye(blocks), jordan)
    combined = np.concatenate([full_semisimple, full_jordan], axis=1)
    rank_s = int(np.linalg.matrix_rank(full_semisimple))
    rank_j = int(np.linalg.matrix_rank(full_jordan))
    rank_combined = int(np.linalg.matrix_rank(combined))
    intersection = rank_s + rank_j - rank_combined

    topological_change = np.asarray([
        [0.5, -0.5, 0.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.0, 1.0],
    ])
    return {
        "single_geometry_raw_basis": {
            "latent_order": ["A_top", "E_top", "B_bulk"],
            "raw_order": ["F1", "F2", "B_bulk"],
            "basis_change_T": matrix_payload(topological_change),
            "determinant_T": 0.5,
            "semisimple_direct_sum_rank": 3,
            "Jordan_Etop_Bbulk_shear_rank": 3,
            "semisimple_image": "R^3",
            "Jordan_image": "R^3",
            "image_intersection": "R^3",
            "conclusion": "single_geometry_exactly_nonidentifiable_after_all_amplitudes_are_retained",
        },
        "dyadic_three_generation_blocks": {
            "block_order": [
                f"{lineage['id']}.F{activation}"
                for lineage in manifest["lineages"]
                for activation in (1, 2)
            ],
            "semisimple_block_X_kappa": matrix_payload(semisimple),
            "Jordan_block_X_J": matrix_payload(jordan),
            "semisimple_block_left_null": [kappa, -(1.0 + kappa), 1.0],
            "Jordan_block_left_null": [1.0, -2.0, 1.0],
            "per_block_rank_semisimple": 2,
            "per_block_rank_Jordan": 2,
            "per_block_intersection_basis": [[1.0, 1.0, 1.0]],
            "per_block_intersection_dimension_for_fixed_kappa_not_1": 1,
            "total_blocks": blocks,
            "total_observation_dimension": 3 * blocks,
            "total_rank_semisimple_fixed_kappa": rank_s,
            "total_rank_Jordan": rank_j,
            "combined_rank": rank_combined,
            "total_image_intersection_dimension": intersection,
            "fixed_kappa_exact_conclusion": "distinct_except_for_each_blocks_pure_leading_constant_line",
        },
        "nonuniform_collision": {
            "parameterization": "kappa=exp(-epsilon), b=-c/epsilon, a_semisimple=a+c/epsilon",
            "limit": "a_semisimple+b*kappa^g -> a+c*g for g=0,1,2 as epsilon->0+",
            "consequence": "the_closure_of_the_unbounded_semisimple_image_contains_the_full_Jordan_image",
            "uniform_separation_requires": [
                "a_physical_spectral_gap_1-kappa_at_least_epsilon",
                "or_an_amplitude_bound",
                "or_the_declared_restricted_trace_modulus_phase_transport_relation",
            ],
        },
    }


def build_report(
    payload: Mapping[str, Any], manifest: Mapping[str, Any],
    input_path: Path, manifest_path: Path
) -> dict[str, Any]:
    validate_input(payload, manifest, input_path)
    block = select_scaled_activation_block(payload, manifest)
    estimate = block["scaled_estimate"]
    covariance = block["scaled_covariance"]
    numerics = manifest["numerics"]
    residual_order = [
        f"{lineage['id']}.F{activation}"
        for lineage in manifest["lineages"]
        for activation in (1, 2)
    ]

    fixed_q2 = score_at_kappa(
        float(manifest["models"]["fixed_semisimple_q2"]["kappa"]),
        estimate, covariance, numerics, fitted_parameters=0,
        residual_order=residual_order,
    )
    h8_contract = manifest["models"]["fixed_h8_alias_transplant"]
    h8_numerics = {
        **numerics,
        "decision_alpha": float(h8_contract["decision_alpha"]),
    }
    fixed_h8_alias = score_at_kappa(
        float(h8_contract["kappa"]),
        estimate, covariance, h8_numerics, fitted_parameters=0,
        residual_order=residual_order,
    )
    jordan = score_at_kappa(
        1.0, estimate, covariance, numerics, fitted_parameters=0,
        residual_order=residual_order,
    )

    def chi_square_at(kappa: float) -> float:
        return chi_square_at_kappa(kappa, estimate, covariance, numerics)[0]

    free_contract = manifest["models"]["free_semisimple_kappa"]
    free_profile = bounded_profile(
        chi_square_at,
        float(free_contract["support"][0]),
        float(free_contract["support"][1]),
        int(free_contract["grid_points"]),
        float(free_contract["refinement_tolerance"]),
    )
    free_score = score_at_kappa(
        float(free_profile["best_kappa"]), estimate, covariance, numerics,
        fitted_parameters=1, residual_order=residual_order,
    )
    free_kappa = float(free_profile["best_kappa"])
    free_profile["relative_exponent_delta_minus_log2_kappa"] = (
        finite_float(-math.log(free_kappa, 2.0)) if free_kappa > 0 else None
    )
    free_profile["score_with_one_profiled_parameter"] = free_score

    physical_contract = manifest["models"]["physical_decaying_semisimple"]
    physical_profile = bounded_profile(
        chi_square_at,
        float(physical_contract["closed_profile_support"][0]),
        float(physical_contract["closed_profile_support"][1]),
        int(physical_contract["grid_points"]),
        float(physical_contract["refinement_tolerance"]),
    )
    physical_score = score_at_kappa(
        float(physical_profile["best_kappa"]), estimate, covariance, numerics,
        fitted_parameters=1, residual_order=residual_order,
    )
    physical_profile.update({
        "physical_open_support": physical_contract["physical_open_support"],
        "score_with_one_profiled_parameter": physical_score,
        "boundary_calibration": "chi_square_df3_is_diagnostic_only_because_the_infimum_is_on_the_kappa_1_collision_boundary",
        "open_interval_result": "no_interior_minimum_the_infimum_is_kappa_to_1_from_below",
        "fixed_boundary_Jordan_score_reference": jordan,
    })

    by_coordinate = []
    for position, row in enumerate(block["selected_order"]):
        by_coordinate.append({
            **row,
            "raw_delta4_F": finite_float(block["raw_estimate"][position]),
            "raw_standard_error": finite_float(math.sqrt(max(
                block["raw_covariance"][position, position], 0.0
            ))),
            "N_power_13_over_8_scaled_u": finite_float(estimate[position]),
            "scaled_standard_error": finite_float(math.sqrt(max(
                covariance[position, position], 0.0
            ))),
        })

    return {
        "schema": OUTPUT_SCHEMA,
        "status": "completed_existing_data_reuse_zero_new_MC",
        "conclusion": {
            "classification": "PARTIALLY_IDENTIFIABLE",
            "uniform_identifiability": False,
            "summary": "fixed_or_spectrally_separated_semisimple_images_are_testable_but_the_unbounded_kappa_to_1_semisimple_closure_contains_the_Jordan_image",
            "not_a_field_identification": True,
        },
        "input": {
            "path": str(input_path.relative_to(ROOT)),
            "schema": payload["schema"],
            "sha256": sha256_file(input_path),
            "source_manifest_sha256": payload["provenance"]["manifest_sha256"],
            "manifest_path": str(manifest_path.relative_to(ROOT)),
            "manifest_sha256": sha256_file(manifest_path),
            "dependency_groups": payload["dependency_groups"],
            "covariance": "complete_selected_submatrix_of_the_pinned_cross_size_aligned_delete_one_covariance",
        },
        "observable": {
            "raw_coordinates": ["angular_delta_F1", "angular_delta_F2"],
            "definition": "d_i(N)=[F_i(first,p_bar)-F_i(second,p_bar)]/Delta_cos_4_theta",
            "units": "dimensionless_probability_response_per_unit_Delta_cos4",
            "root": "same_pooled_moving_matching_root_as_the_input_archive",
            "scaled_coordinates": "u_i(N)=N^(13/8)d_i(N)",
            "redundant_coordinates_not_counted_as_new_rank": ["A_top=d1+d2", "E_top=d2-d1", "delta_p_i=-d_i/Mbar_prime"],
            "selected_order": block["selected_order"],
            "estimate_by_coordinate": by_coordinate,
            "raw_covariance": matrix_payload(block["raw_covariance"]),
            "scaled_covariance": matrix_payload(covariance),
        },
        "model_contracts": {
            "semisimple": {
                "equation": "u_li(g)=a_li+b_li*kappa^g_for_g_0_1_2",
                "allowed_amplitudes": "independent_a_li_and_b_li_for_each_lineage_l_and_activation_i",
                "shared_parameter": "one_common_kappa_across_all_four_blocks",
                "fixed_q2": "kappa=1/2_equivalently_relative_N_minus_1_mode",
                "fixed_H8_alias_transplant": "kappa=2^(-11/8)_from_transporting_a_canonical_weight8_mode_relative_to_the_weight21/4_global_Q4_candidate",
                "physical_decaying_range": "0<kappa<1",
            },
            "Jordan": {
                "equation": "u_li(g)=a_li+c_li*g_for_g_0_1_2",
                "allowed_amplitudes": "independent_a_li_and_c_li_for_each_lineage_l_and_activation_i",
                "reference_scale_gauge": "changing_N0_shifts_a_li_by_c_li_times_a_constant_without_changing_c_li",
                "nilpotent_gauge": "only_mu_times_bottom_overlap_is_identifiable_under_mu_to_s_mu_and_bottom_to_bottom_over_s",
            },
        },
        "design_audit": design_audit(manifest),
        "scores": {
            "fixed_semisimple_q2_kappa_0p5": fixed_q2,
            "fixed_H8_alias_radial_transplant_kappa_2_pow_minus_11_over_8": {
                **fixed_h8_alias,
                "hypothesis_provenance": {
                    "status": h8_contract["angular_source_status"],
                    "commit": h8_contract["angular_source_commit"],
                    "path": h8_contract["angular_source_path"],
                    "angular_decision": h8_contract["angular_source_decision"],
                    "interpretation_commit": h8_contract["interpretation_commit"],
                    "interpretation_path": h8_contract["interpretation_path"],
                    "post_reveal_scalar_alias": h8_contract["post_reveal_scalar_alias"],
                },
                "relative_N_exponent": h8_contract["relative_N_exponent"],
                "assumptions": [
                    "transport_the_primitive_real_C3_H8_alias_into_the_global_K1_K2_residual",
                    "use_canonical_weight8_radial_scaling",
                    "take_the_global_leading_candidate_to_have_dimension21_over_4",
                    "allow_independent_leading_and_subleading_amplitudes_in_all_four_blocks",
                ],
                "claim_boundary": "excludes_only_the_naive_radial_transplant_envelope_not_the_branch_only_primitive_C3_phase_result_or_every_H8_mixed_mechanism",
            },
            "Jordan_kappa_1": jordan,
            "free_kappa_unconstrained": free_profile,
            "physical_decaying_semisimple_0_lt_kappa_lt_1": physical_profile,
        },
        "unique_missing_input": {
            "name": "restricted_trace_modulus_or_phase_transport",
            "requirement": "one_semantics_matched_transfer_relation_for_the_same_B_source_original_q_E_and_pooled_root_physical_normalizer_on_rank0_rank2_restricted_traces",
            "preferred_route": "derive_the_independent_singlet_F_B_tau_and_Jordan_top_F_tilde_4_tau_vectors_then_score_them_on_existing_rho_child_or_P43_P57_assets",
            "fallback_acquisition": "one_phase_calibrated_second_physical_rotation_of_the_same_B_column",
            "why_needed": "it_fixes_the_bulk_column_eigenvalue_or_character_and_removes_the_kappa_to_1_collision_freedom",
        },
        "boundaries": [
            "all_scores_are_post_reveal_existing_archive_diagnostics",
            "no_new_Monte_Carlo_configuration_or_counter_is_generated",
            "q2_exclusion_is_only_the_fixed_kappa_1_over_2_parameterization",
            "the_fixed_H8_alias_radial_transplant_is_a_post_reveal_sector_informed_diagnostic",
            "its_exclusion_does_not_reject_the_primitive_real_C3_H8_phase_result_or_every_H8_mixture",
            "Jordan_non_exclusion_is_not_Jordan_operator_identification",
            "the_free_kappa_best_fit_is_not_a_free_exponent_claim",
            "priority_is_attention_not_permission_or_a_task_lock",
        ],
    }


def markdown_report(report: Mapping[str, Any]) -> str:
    scores = report["scores"]
    q2 = scores["fixed_semisimple_q2_kappa_0p5"]
    h8 = scores["fixed_H8_alias_radial_transplant_kappa_2_pow_minus_11_over_8"]
    jordan = scores["Jordan_kappa_1"]
    free = scores["free_kappa_unconstrained"]
    physical = scores["physical_decaying_semisimple_0_lt_kappa_lt_1"]
    design = report["design_audit"]
    missing = report["unique_missing_input"]

    lines = [
        "# Issue #275 forward identifiability on existing K1/K2 production",
        "",
        "## Decision",
        "",
        "The two mechanism classes are **PARTIALLY_IDENTIFIABLE, but not uniformly identifiable**.",
        "A fixed or spectrally separated semisimple second mode has a different three-generation",
        "image from a Jordan shear.  However, after all amplitudes are retained and no amplitude",
        "bound or spectral gap is invented, the `kappa -> 1` closure of the semisimple image",
        "contains the full Jordan image.  The present result therefore does not identify a Jordan",
        "operator.",
        "",
        "No Monte Carlo sample is generated.  The score consumes the complete selected submatrix",
        "of the pinned cross-size aligned-delete-one covariance.",
        "",
        "## Existing-data scores",
        "",
        "| model | kappa | chi-square / df | survival p | reading |",
        "|:--|--:|--:|--:|:--|",
        f"| fixed ordinary semisimple `q2` | {q2['kappa']:.6g} | {q2['mahalanobis_chi_square']:.6f} / {q2['degrees_of_freedom']} | {q2['chi_square_survival_p']:.6g} | {'excluded at .05' if q2['excluded_at_alpha'] else 'not excluded'} |",
        f"| fixed primitive-C3 H8 radial transplant | {h8['kappa']:.6g} | {h8['mahalanobis_chi_square']:.6f} / {h8['degrees_of_freedom']} | {h8['chi_square_survival_p']:.6g} | {'excluded at .01' if h8['excluded_at_alpha'] else 'not excluded at .01'} |",
        f"| Jordan affine log | {jordan['kappa']:.6g} | {jordan['mahalanobis_chi_square']:.6f} / {jordan['degrees_of_freedom']} | {jordan['chi_square_survival_p']:.6g} | {'excluded at .05' if jordan['excluded_at_alpha'] else 'not excluded'} |",
        f"| free kappa, unconstrained | {free['best_kappa']:.9f} | {free['score_with_one_profiled_parameter']['mahalanobis_chi_square']:.6f} / {free['score_with_one_profiled_parameter']['degrees_of_freedom']} | {free['score_with_one_profiled_parameter']['chi_square_survival_p']:.6g} | descriptive optimum |",
        f"| physical decaying semisimple | `0<kappa<1` | infimum {physical['best_chi_square']:.6f} at `kappa -> 1-` | diagnostic df3 p={physical['score_with_one_profiled_parameter']['chi_square_survival_p']:.6g} | collision boundary, no interior winner |",
        "",
        f"The unconstrained optimum is `kappa={free['best_kappa']:.9f}`, equivalent to",
        f"`delta=-log2(kappa)={free['relative_exponent_delta_minus_log2_kappa']:.9f}`.  Its negative",
        "relative exponent means that the second mode grows relative to the proposed leading",
        "`N^-13/8` response; it is not a more-irrelevant bulk singlet.",
        "",
        "The frozen binary branch-only paired physical-rotation gate at `0b9e89c9` selects",
        "the H8 line over H4 for its primitive real-C3 observer.  A post-reveal H0",
        "signed-real line also survives (`p=.250468`), so the gate identifies an H8/even",
        "branch rather than a unique local H8 field.  A deliberately generous radial transplant of",
        "that alias into the global K1/K2 residual takes a canonical weight-8 correction",
        "relative to the weight-21/4 Q4 candidate, hence `kappa=2^(-11/8)`.  Even with",
        "independent leading and subleading amplitudes in all four blocks, this envelope is",
        f"excluded (`chi2={h8['mahalanobis_chi_square']:.6f}/{h8['degrees_of_freedom']}`;",
        f"nominal `p={h8['chi_square_survival_p']:.6g}`).  The H8 member of the surviving",
        "primitive-sector branch therefore cannot be copied into the global residual as one unmixed canonical",
        "weight-8 radial mode.  This does not reject that branch result or H8-containing",
        "observer-mixing mechanisms.",
        "",
        "## Raw coordinates and units",
        "",
        "For activation `i=1,2`,",
        "",
        "```text",
        "d_i(N) = [F_i(first,p_bar)-F_i(second,p_bar)] / Delta cos(4 theta),",
        "u_i(N) = N^(13/8) d_i(N).",
        "```",
        "",
        "`d_i` is a dimensionless probability response per unit `Delta cos(4 theta)` at the",
        "same pooled moving matching root as the source archive.  `A_top=d1+d2`,",
        "`E_top=d2-d1`, and the linearized root shifts are deterministic transforms and are not",
        "counted as extra evidence coordinates.",
        "",
        "## Exact forward maps",
        "",
        "For each lineage and activation, the fixed-`kappa` semisimple design and Jordan design are",
        "",
        "```text",
        "X_S(kappa) = [[1,1], [1,kappa], [1,kappa^2]],",
        "X_J        = [[1,0], [1,1],     [1,2]].",
        "```",
        "",
        "Thus, from the first two generations,",
        "",
        "```text",
        "semisimple: d_hat_i(4N) = (4N)^(-13/8) *",
        "             [(1+kappa)(2N)^(13/8)d_i(2N) - kappa N^(13/8)d_i(N)],",
        "Jordan:     d_hat_i(4N) = (4N)^(-13/8) *",
        "             [2(2N)^(13/8)d_i(2N) - N^(13/8)d_i(N)].",
        "```",
        "",
        f"There are {design['dyadic_three_generation_blocks']['total_blocks']} blocks.  At fixed",
        f"`kappa != 1`, both complete designs have rank {design['dyadic_three_generation_blocks']['total_rank_Jordan']};",
        f"their combined rank is {design['dyadic_three_generation_blocks']['combined_rank']} and",
        f"their image intersection has dimension {design['dyadic_three_generation_blocks']['total_image_intersection_dimension']}.",
        "Per block the intersection is exactly the pure leading constant line `(1,1,1)`.",
        "",
        "At a single geometry the change of basis from `(A_top,E_top,B_bulk)` to",
        "`(F1,F2,B_bulk)` has determinant `1/2`.  Both the semisimple direct sum and an",
        "`E_top/B_bulk` Jordan shear are rank three and have image `R^3`; adding an arbitrary",
        "bulk column without a transport law is therefore exactly nonidentifying.",
        "",
        "## Why separation is non-uniform",
        "",
        "Write `kappa=exp(-epsilon)`, `b=-c/epsilon`, and",
        "`a_semisimple=a+c/epsilon`.  For generation `g=0,1,2`,",
        "",
        "```text",
        "a_semisimple + b kappa^g  ->  a + c g.",
        "```",
        "",
        "Therefore the closure of the unbounded semisimple image contains the Jordan plane.",
        "The fixed `q2` rejection removes only `kappa=1/2`; it does not eliminate every distinct",
        "bulk singlet whose transfer eigenvalue is allowed to approach the topological one.",
        "",
        "## Unique missing physical input",
        "",
        f"The missing input is **{missing['name']}**: {missing['requirement'].replace('_', ' ')}.",
        "The preferred route is a theory-derived restricted-trace vector for the independent",
        "singlet and Jordan top component, scored on existing rho-child or P43/P57 assets.  If",
        "that relation cannot be derived, the only new acquisition justified by this audit is one",
        "phase-calibrated second physical rotation of the same `B` column—not another untyped",
        "topological coordinate.",
        "",
        "## Provenance and boundary",
        "",
        f"- Input: `{report['input']['path']}` (`{report['input']['sha256']}`).",
        f"- Manifest: `{report['input']['manifest_path']}` (`{report['input']['manifest_sha256']}`).",
        "- All four residual rows are correlated views of the registered dependency groups.",
        "- Scores are post-reveal existing-data diagnostics, not prospective validation.",
        "- Jordan compatibility is not continuum-field, lattice-overlap, or normalizer identification.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--output-md", type=Path)
    args = parser.parse_args()

    manifest_path = args.manifest.resolve()
    manifest = load_manifest(manifest_path)
    input_path = resolve_repo_path(str(manifest["input"]["path"])).resolve()
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    report = build_report(payload, manifest, input_path, manifest_path)

    output_json = (args.output_json or resolve_repo_path(manifest["outputs"]["json"])).resolve()
    output_md = (args.output_md or resolve_repo_path(manifest["outputs"]["markdown"])).resolve()
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    output_md.write_text(markdown_report(report), encoding="utf-8")


if __name__ == "__main__":
    main()
