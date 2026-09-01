#!/usr/bin/env python3
"""Score the P250 old-radius-four plus fresh-degree-five joint operator."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path
import subprocess
import sys
from typing import Mapping, Sequence

import numpy as np
import scipy
from scipy.stats import chi2, f

from score_z5_projective_leg_joint_annihilation import (
    MONOMIALS_2,
    candidate_matrix,
    complex_payload,
    hand_hankel,
    jackknife_covariance,
    means,
    pair,
    pinned_bytes,
    projective_q,
    read_batches,
    reflect,
    residual_labels,
    rotate_power,
    sha256_bytes,
    transformed_basis,
)


ROWS_3 = ((3, 0), (2, 1), (1, 2), (0, 3))
RANK = 5
EIGEN_RELATIVE_CUTOFF = 1e-10


def validate_stream(rows: Sequence[dict], spec: Mapping[str, object]) -> dict:
    expected_batches = int(spec["batches"])
    expected_samples = int(spec["samples"])
    expected_per_batch = int(spec["samples_per_batch"])
    if len(rows) != expected_batches:
        raise ValueError(f"{spec['name']}: expected {expected_batches} batches, found {len(rows)}")
    ids = [int(row["batch"]) for row in rows]
    if ids != list(range(expected_batches)):
        raise ValueError(f"{spec['name']}: batch ids are not contiguous 0..B-1")
    batch_sizes = {int(row["samples"]) for row in rows}
    if batch_sizes != {expected_per_batch}:
        raise ValueError(f"{spec['name']}: unexpected batch sizes {sorted(batch_sizes)}")
    if sum(int(row["samples"]) for row in rows) != expected_samples:
        raise ValueError(f"{spec['name']}: sample total changed")
    for index, row in enumerate(rows):
        if int(row["replica_first"]) != index * expected_per_batch:
            raise ValueError(f"{spec['name']}: replica intervals changed at batch {index}")
    return {
        "batch_ids": [ids[0], ids[-1]],
        "batches": expected_batches,
        "samples": expected_samples,
        "samples_per_batch": expected_per_batch,
        "field_hashes_unique": len({row["field_sha256"] for row in rows}),
        "translation_hashes_unique": len({row["translation_sha256"] for row in rows}),
    }


def pinned_json(spec: Mapping[str, object], label: str) -> tuple[dict, dict]:
    payload = pinned_bytes(Path(spec["path"]), spec["commit"])
    observed = sha256_bytes(payload)
    if observed != spec["sha256"]:
        raise ValueError(f"{label} hash changed")
    return json.loads(payload), {**spec, "observed_sha256": observed}


def pinned_text(spec: Mapping[str, object], label: str) -> tuple[str, dict]:
    payload = pinned_bytes(Path(spec["path"]), spec["commit"])
    observed = sha256_bytes(payload)
    if observed != spec["sha256"]:
        raise ValueError(f"{label} hash changed")
    return payload.decode("utf-8"), {**spec, "observed_sha256": observed}


def validate_input_semantics(
    manifest: Mapping[str, object],
    radius5_reference: Mapping[str, object],
) -> dict:
    expected = manifest["shared_semantics"]["expected"]
    old_response, old_response_audit = pinned_json(
        manifest["old_stream"]["response"], "old response",
    )
    old_gate, old_gate_audit = pinned_json(
        manifest["old_stream"]["exact_gate"], "old exact gate",
    )
    fresh_response, fresh_response_audit = pinned_json(
        manifest["fresh_stream"]["response"], "fresh response",
    )
    fresh_gate, fresh_gate_audit = pinned_json(
        manifest["fresh_stream"]["exact_gate"], "fresh exact gate",
    )
    old_runner, old_runner_audit = pinned_text(
        manifest["shared_semantics"]["old_runner"], "old runner",
    )
    fresh_runner, fresh_runner_audit = pinned_text(
        manifest["shared_semantics"]["fresh_runner"], "fresh runner",
    )
    shared_runner, shared_runner_audit = pinned_text(
        manifest["shared_semantics"]["shared_geometry_runner"], "shared geometry runner",
    )
    if (
        old_response["manifest_runner_commit"]
        != manifest["shared_semantics"]["old_runner"]["declared_manifest_runner_commit"]
    ):
        raise ValueError("old response runner provenance changed")
    if (
        fresh_response["manifest_runner_commit"]
        != manifest["shared_semantics"]["fresh_runner"]["declared_manifest_runner_commit"]
    ):
        raise ValueError("fresh response runner provenance changed")

    parent_order = int(expected["parent_order"])
    child_order = int(expected["child_order"])
    p_fixed = float(expected["p"])
    hands = list(expected["hands"])
    charges = list(expected["charges"])
    if old_response["observable"]["parent_order"] != parent_order:
        raise ValueError("old response parent order changed")
    if old_response["observable"]["child_order"] != child_order:
        raise ValueError("old response is not N505")
    if old_gate["cross_scale_gate"]["child_order"] != child_order:
        raise ValueError("old exact gate is not N505")
    if old_gate["cross_scale_gate"]["parent_order"] != parent_order:
        raise ValueError("old exact gate parent order changed")
    if f"CHILD_ORDER = {child_order}" not in shared_runner:
        raise ValueError("shared geometry runner no longer fixes N505")
    for token in ("CHILD_ORDER", "contexts", "range(CHILD_ORDER)"):
        if token not in fresh_runner:
            raise ValueError(f"fresh runner lost N505 evidence token {token}")
    if old_response["run"]["p"] != p_fixed or fresh_response["run"]["p"] != p_fixed:
        raise ValueError("old/fresh p differs from the frozen common value")
    if old_response["observable"]["hands"] != hands or fresh_response["observable"]["hands"] != hands:
        raise ValueError("old/fresh hand order changed")
    if old_response["observable"]["charges"] != charges or fresh_response["observable"]["charges"] != charges:
        raise ValueError("old/fresh charge order changed")
    if fresh_gate["hands"] != hands or fresh_gate["charges"] != charges:
        raise ValueError("fresh exact gate hand/charge scope changed")
    if not old_gate["passed"] or not fresh_gate["passed"]:
        raise ValueError("one input exact gate is not passing")
    if fresh_gate["radius"] != 5 or fresh_gate["shell_points"] != 20:
        raise ValueError("fresh exact gate is not the radius-five shell")
    if fresh_response["observable"]["shell"] != fresh_gate["shell"]:
        raise ValueError("fresh response/exact-gate shell mismatch")

    gauge = manifest["shared_semantics"]["gauge"]
    if old_response["observable"]["gauge"] != gauge["old_response_label"]:
        raise ValueError("old response gauge label changed")
    rotation_gate = old_gate["rotation_fiber_gate"]
    if not rotation_gate["passed"] or rotation_gate["gauge_equation"] != gauge["equation"]:
        raise ValueError("old affine-fiber gauge gate changed")
    for token in ("rotation_gauges", "gauge_charged_rows"):
        if token not in old_runner or token not in fresh_runner:
            raise ValueError(f"old/fresh runners no longer share {token}")

    translation = manifest["shared_semantics"]["translation"]
    for source, token, label in (
        (old_runner, translation["old_salt_token"], "old"),
        (fresh_runner, translation["fresh_salt_token"], "fresh"),
    ):
        if token not in source or "% PARENT_GEOMETRY.n" not in source:
            raise ValueError(f"{label} translation convention changed")
        if "translation_digest.update" not in source:
            raise ValueError(f"{label} translation audit digest disappeared")
    if translation["old_salt_token"] == translation["fresh_salt_token"]:
        raise ValueError("old and fresh translation salts must remain distinct")

    aliases = manifest["radius5_alias_metadata"]
    candidate_by_name = {row["name"]: row for row in manifest["candidate_maps"]}
    expected_alias_values = set(radius5_reference["parameter_free_maps"])
    if set(aliases.values()) != expected_alias_values:
        raise ValueError("radius-five alias set no longer matches the published score")
    for name, alias in aliases.items():
        candidate = candidate_by_name[name]
        if not candidate["coefficient_conjugation"]:
            raise ValueError(f"radius-five alias {alias} lost coefficient conjugation")
        if alias == "identity_conjugation":
            valid = not candidate["alexander_reflection"] and candidate["rotation_power"] == 0
        else:
            valid = (
                candidate["alexander_reflection"]
                and alias == f"Alexander_R{candidate['rotation_power']}_conjugation"
            )
        if not valid or candidate.get("radius5_alias") != alias:
            raise ValueError(f"radius-five alias metadata changed for {name}")
    return {
        "passed": True,
        "parent_order": parent_order,
        "child_order": child_order,
        "p": p_fixed,
        "hands": hands,
        "charges": charges,
        "gauge": {
            "old_exact_gate_passed": True,
            "shared_rotation_gauge_functions_pinned": True,
            "equation": gauge["equation"],
        },
        "translation": {
            "same_uniform_parent_translation_convention": True,
            "independent_seed_and_salt": True,
            "old_salt_token": translation["old_salt_token"],
            "fresh_salt_token": translation["fresh_salt_token"],
        },
        "radius5_aliases": aliases,
        "artifacts": {
            "old_response": old_response_audit,
            "old_exact_gate": old_gate_audit,
            "fresh_response": fresh_response_audit,
            "fresh_exact_gate": fresh_gate_audit,
            "old_runner": old_runner_audit,
            "fresh_runner": fresh_runner_audit,
            "shared_geometry_runner": shared_runner_audit,
        },
        "fresh_N505_evidence": (
            "the pinned fresh runner imports CHILD_ORDER and contexts from the pinned "
            "shared geometry runner, where CHILD_ORDER=505; its response pins that runner commit"
        ),
    }


def transformed_rows(candidate: Mapping[str, object]) -> tuple[tuple[int, int], ...]:
    return tuple(
        rotate_power(
            reflect(point) if bool(candidate["alexander_reflection"]) else point,
            int(candidate["rotation_power"]),
        )
        for point in ROWS_3
    )


def mixed_moment(
    old_values: Mapping[str, float],
    fresh_values: Mapping[str, float],
    point: tuple[int, int],
    hand: str,
    charge: int,
) -> complex:
    degree = abs(point[0]) + abs(point[1])
    if degree <= 4:
        return pair(old_values, point, hand, charge)
    if degree == 5:
        return pair(fresh_values, point, hand, charge)
    raise ValueError(f"point {point} is outside the frozen radius-five domain")


def extension_matrix(
    old_values: Mapping[str, float],
    fresh_values: Mapping[str, float],
    hand: str,
    basis: Sequence[tuple[int, int]],
    left_rows: Sequence[tuple[int, int]],
) -> np.ndarray:
    return np.asarray(
        [
            [
                mixed_moment(
                    old_values,
                    fresh_values,
                    (left[0] + right[0], left[1] + right[1]),
                    hand,
                    charge,
                )
                for right in basis
            ]
            for charge in (1, 2)
            for left in left_rows
        ],
        dtype=complex,
    )


def extension_row_labels(
    hand: str,
    left_rows: Sequence[tuple[int, int]],
    *,
    conjugated: bool,
) -> list[str]:
    labels: list[str] = []
    for charge in (1, 2):
        for left in left_rows:
            prefix = (
                f"extension:{hand}:r{charge}:u3({left[0]},{left[1]})"
                + (":conjugated" if conjugated else "")
            )
            labels.extend((prefix + ":re", prefix + ":im"))
    return labels


def candidate_blocks(
    old_values: Mapping[str, float],
    fresh_values: Mapping[str, float],
    candidate: Mapping[str, object],
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    old_matrix, _ = candidate_matrix(old_values, candidate)
    target_basis = transformed_basis(
        alexander_reflection=bool(candidate["alexander_reflection"]),
        rotation_power=int(candidate["rotation_power"]),
    )
    target_rows = transformed_rows(candidate)
    conjugated = bool(candidate["coefficient_conjugation"])
    plus_extension = extension_matrix(
        old_values, fresh_values, "plus", MONOMIALS_2, ROWS_3,
    )
    if conjugated:
        plus_extension = plus_extension.conjugate()
    minus_extension = extension_matrix(
        old_values, fresh_values, "minus", target_basis, target_rows,
    )
    extension = np.vstack((plus_extension, minus_extension))
    labels = (
        extension_row_labels("plus", ROWS_3, conjugated=conjugated)
        + extension_row_labels("minus", target_rows, conjugated=False)
    )
    if old_matrix.shape != (24, 6) or extension.shape != (16, 6):
        raise ValueError("augmented candidate matrices changed shape")
    return old_matrix, extension, labels


def augmented_residual_from_blocks(
    old_matrix: np.ndarray,
    extension: np.ndarray,
    pivot: Mapping[str, object],
) -> tuple[np.ndarray, np.ndarray]:
    """Return the locally complete 35-complex rank-five residual.

    The five pivot equations and five pivot columns are frozen from the old
    24x6 candidate matrix.  Hence the first 19 complex coordinates exactly
    replay the old joint-annihilation chart and all 16 extension rows are
    honest additional equations for the same projective coefficient vector.
    """
    rows = tuple(int(value) for value in pivot["rows"])
    columns = tuple(int(value) for value in pivot["columns"])
    if len(rows) != RANK or len(columns) != RANK:
        raise ValueError("frozen pivot is not rank five")
    if any(row >= old_matrix.shape[0] for row in rows):
        raise ValueError("frozen pivot must use only old rows")
    coefficients = projective_q(old_matrix, pivot)
    other_rows = [index for index in range(old_matrix.shape[0]) if index not in rows]
    residual = np.concatenate((old_matrix[other_rows] @ coefficients, extension @ coefficients))
    expected = old_matrix.shape[0] - RANK + extension.shape[0]
    if residual.shape != (expected,):
        raise ValueError("augmented residual dimension changed")
    return residual, coefficients


def realify(vector: Sequence[complex]) -> np.ndarray:
    values = np.asarray(vector, dtype=complex).reshape(-1)
    return np.column_stack((values.real, values.imag)).reshape(-1)


def candidate_residual(
    old_values: Mapping[str, float],
    fresh_values: Mapping[str, float],
    candidate: Mapping[str, object],
    pivot: Mapping[str, object],
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    old_matrix, extension, extension_labels = candidate_blocks(old_values, fresh_values, candidate)
    complex_residual, coefficients = augmented_residual_from_blocks(old_matrix, extension, pivot)
    old_row_labels = candidate_matrix(old_values, candidate)[1]
    labels = residual_labels(old_row_labels, pivot) + extension_labels
    if len(labels) != 70:
        raise ValueError("augmented residual labels changed")
    return realify(complex_residual), coefficients, labels


def hand_blocks(
    old_values: Mapping[str, float],
    fresh_values: Mapping[str, float],
    candidate: Mapping[str, object],
    hand: str,
) -> tuple[np.ndarray, np.ndarray, Sequence[tuple[int, int]], Sequence[tuple[int, int]], bool]:
    if hand == "plus":
        basis = MONOMIALS_2
        left_rows = ROWS_3
        conjugated = bool(candidate["coefficient_conjugation"])
    elif hand == "minus":
        basis = transformed_basis(
            alexander_reflection=bool(candidate["alexander_reflection"]),
            rotation_power=int(candidate["rotation_power"]),
        )
        left_rows = transformed_rows(candidate)
        conjugated = False
    else:
        raise ValueError(f"unknown hand {hand}")
    old_matrix = hand_hankel(old_values, hand, basis)
    extension = extension_matrix(old_values, fresh_values, hand, basis, left_rows)
    if conjugated:
        old_matrix = old_matrix.conjugate()
        extension = extension.conjugate()
    if old_matrix.shape != (12, 6) or extension.shape != (8, 6):
        raise ValueError("hand-specific augmented matrices changed shape")
    return old_matrix, extension, basis, left_rows, conjugated


def hand_residual(
    old_values: Mapping[str, float],
    fresh_values: Mapping[str, float],
    candidate: Mapping[str, object],
    hand: str,
    pivot: Mapping[str, object],
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    old_matrix, extension, basis, left_rows, conjugated = hand_blocks(
        old_values, fresh_values, candidate, hand,
    )
    complex_residual, coefficients = augmented_residual_from_blocks(old_matrix, extension, pivot)
    # A hand block has seven old nonpivot rows and eight extension rows.
    if complex_residual.shape != (15,):
        raise ValueError("hand-specific residual is not 15 complex coordinates")
    old_labels = []
    pivot_rows = set(int(value) for value in pivot["rows"])
    row_index = 0
    for charge in (1, 2):
        for left in basis:
            if row_index not in pivot_rows:
                prefix = (
                    f"old:{hand}:r{charge}:u({left[0]},{left[1]})"
                    + (":conjugated" if conjugated else "")
                )
                old_labels.extend((prefix + ":re", prefix + ":im"))
            row_index += 1
    labels = old_labels + extension_row_labels(hand, left_rows, conjugated=conjugated)
    return realify(complex_residual), coefficients, labels


def covariance_score_from_components(
    point: Sequence[float],
    old_deleted: Sequence[Sequence[float]],
    fresh_deleted: Sequence[Sequence[float]],
    *,
    eigen_relative_cutoff: float,
) -> dict:
    old_covariance = jackknife_covariance(old_deleted)
    fresh_covariance = jackknife_covariance(fresh_deleted)
    covariance = old_covariance + fresh_covariance
    point_array = np.asarray(point, dtype=float)
    scales = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    if np.any(scales <= 0.0):
        raise ValueError("augmented residual has a zero total covariance scale")
    correlation = covariance / scales[:, None] / scales[None, :]
    correlation = 0.5 * (correlation + correlation.T)
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    keep = eigenvalues > eigen_relative_cutoff * eigenvalues[-1]
    if not np.any(keep):
        raise ValueError("augmented covariance has no resolved mode")
    normalized = point_array / scales
    projection = eigenvectors[:, keep].T @ normalized
    statistic = float(np.sum(projection**2 / eigenvalues[keep]))
    degrees = int(np.count_nonzero(keep))
    batches = min(len(old_deleted), len(fresh_deleted))
    if degrees >= batches:
        hotelling_f = math.inf
        denominator_df = batches - degrees
        finite_p = 0.0
    else:
        hotelling_f = (batches - degrees) * statistic / (degrees * (batches - 1))
        denominator_df = batches - degrees
        finite_p = float(f.sf(hotelling_f, degrees, denominator_df))
    return {
        "residual": point_array.tolist(),
        "old_influence_covariance": old_covariance.tolist(),
        "fresh_influence_covariance": fresh_covariance.tolist(),
        "total_covariance": covariance.tolist(),
        "covariance_identity_max_abs_error": float(
            np.max(np.abs(covariance - old_covariance - fresh_covariance))
        ),
        "correlation_eigenvalues": eigenvalues.tolist(),
        "eigen_relative_cutoff": eigen_relative_cutoff,
        "resolved_covariance_modes": degrees,
        "discarded_covariance_modes": len(point_array) - degrees,
        "asymptotic_chi_square": statistic,
        "asymptotic_degrees_of_freedom": degrees,
        "asymptotic_survival_p": float(chi2.sf(statistic, degrees)),
        "finite_batch_Hotelling_F": hotelling_f,
        "finite_batch_numerator_df": degrees,
        "finite_batch_denominator_df": denominator_df,
        "finite_batch_survival_p": finite_p,
        "finite_batch_calibration": (
            "project-compatible min-stream-B Hotelling approximation for a sum of two "
            "independent jackknife covariance estimates; not an exact two-sample pivot law"
        ),
    }


def replay_audit(
    point: np.ndarray,
    score: Mapping[str, object],
    reference_block: Mapping[str, object],
    old_width: int,
) -> dict:
    expected_score = reference_block["score"]
    expected_point = np.asarray(expected_score["residual"], dtype=float)
    observed_old_covariance = np.asarray(score["old_influence_covariance"], dtype=float)[:old_width, :old_width]
    expected_covariance = np.asarray(expected_score["covariance"], dtype=float)
    fresh_old_covariance = np.asarray(score["fresh_influence_covariance"], dtype=float)[:old_width, :old_width]
    point_difference = float(np.max(np.abs(point[:old_width] - expected_point)))
    covariance_difference = float(np.max(np.abs(observed_old_covariance - expected_covariance)))
    fresh_old_max = float(np.max(np.abs(fresh_old_covariance)))
    passed = point_difference <= 1e-14 and covariance_difference <= 1e-14 and fresh_old_max <= 1e-25
    return {
        "passed": passed,
        "old_residual_width": old_width,
        "old_point_max_abs_difference": point_difference,
        "old_covariance_max_abs_difference": covariance_difference,
        "fresh_contribution_to_old_block_max_abs": fresh_old_max,
    }


def score_hand(
    old_full: Mapping[str, float],
    fresh_full: Mapping[str, float],
    old_deleted_values: Sequence[Mapping[str, float]],
    fresh_deleted_values: Sequence[Mapping[str, float]],
    candidate: Mapping[str, object],
    hand: str,
    reference_block: Mapping[str, object],
    alpha: float,
    cutoff: float,
) -> dict:
    pivot = reference_block["pivot"]
    point, coefficients, labels = hand_residual(old_full, fresh_full, candidate, hand, pivot)
    old_deleted = [
        hand_residual(row, fresh_full, candidate, hand, pivot)[0]
        for row in old_deleted_values
    ]
    fresh_deleted = [
        hand_residual(old_full, row, candidate, hand, pivot)[0]
        for row in fresh_deleted_values
    ]
    score = covariance_score_from_components(
        point, old_deleted, fresh_deleted, eigen_relative_cutoff=cutoff,
    )
    replay = replay_audit(point, score, reference_block, old_width=14)
    if not replay["passed"]:
        raise ValueError(f"{hand} old-block replay failed: {replay}")
    return {
        "matrix_shape": [20, 6],
        "old_matrix_shape": [12, 6],
        "extension_matrix_shape": [8, 6],
        "frozen_old_pivot": pivot,
        "projective_q": complex_payload(coefficients),
        "residual_order": labels,
        "replay": replay,
        "score": score,
        "decision": "survives" if score["finite_batch_survival_p"] >= alpha else "rejected",
    }


def score_candidate(
    old_full: Mapping[str, float],
    fresh_full: Mapping[str, float],
    old_deleted_values: Sequence[Mapping[str, float]],
    fresh_deleted_values: Sequence[Mapping[str, float]],
    candidate: Mapping[str, object],
    reference_candidate: Mapping[str, object],
    alpha: float,
    cutoff: float,
) -> tuple[dict, np.ndarray, np.ndarray]:
    pivot = reference_candidate["joint_rank5"]["pivot"]
    point, coefficients, labels = candidate_residual(old_full, fresh_full, candidate, pivot)
    old_deleted = np.asarray([
        candidate_residual(row, fresh_full, candidate, pivot)[0]
        for row in old_deleted_values
    ])
    fresh_deleted = np.asarray([
        candidate_residual(old_full, row, candidate, pivot)[0]
        for row in fresh_deleted_values
    ])
    score = covariance_score_from_components(
        point, old_deleted, fresh_deleted, eigen_relative_cutoff=cutoff,
    )
    replay = replay_audit(point, score, reference_candidate["joint_rank5"], old_width=38)
    if not replay["passed"]:
        raise ValueError(f"{candidate['name']} old joint replay failed: {replay}")
    finite_p = float(score["finite_batch_survival_p"])
    output = {
        "scope": candidate["scope"],
        "radius5_alias": candidate.get("radius5_alias"),
        "alexander_reflection": bool(candidate["alexander_reflection"]),
        "rotation_power": int(candidate["rotation_power"]),
        "coefficient_conjugation": bool(candidate["coefficient_conjugation"]),
        "matrix_shape": [40, 6],
        "old_matrix_shape": [24, 6],
        "extension_matrix_shape": [16, 6],
        "rank_null": "rank(candidate_mapped_old_plus_degree5_operator)<=5",
        "frozen_old_pivot": pivot,
        "projective_q": complex_payload(coefficients),
        "residual_order": labels,
        "replay": replay,
        "joint_score": score,
        "decision_alpha": alpha,
        "decision": "survives" if finite_p >= alpha else "rejected",
    }
    return output, old_deleted, fresh_deleted


def score(
    old_rows: Sequence[dict],
    fresh_rows: Sequence[dict],
    manifest: Mapping[str, object],
    reference: Mapping[str, object],
) -> dict:
    old_full = means(old_rows)
    fresh_full = means(fresh_rows)
    old_deleted_values = [means(old_rows, index) for index in range(len(old_rows))]
    fresh_deleted_values = [means(fresh_rows, index) for index in range(len(fresh_rows))]
    alpha = float(manifest["decision_alpha"])
    cutoff = float(manifest["covariance"]["correlation_eigen_relative_cutoff"])

    candidates = {row["name"]: row for row in manifest["candidate_maps"]}
    results: dict[str, dict] = {}
    old_influences: dict[str, np.ndarray] = {}
    fresh_influences: dict[str, np.ndarray] = {}
    for candidate in manifest["candidate_maps"]:
        name = candidate["name"]
        if name not in reference["candidate_maps"]:
            raise ValueError(f"candidate {name} is absent from the pinned joint score")
        result, old_deleted, fresh_deleted = score_candidate(
            old_full, fresh_full, old_deleted_values, fresh_deleted_values,
            candidate, reference["candidate_maps"][name], alpha, cutoff,
        )
        results[name] = result
        old_influences[name] = old_deleted
        fresh_influences[name] = fresh_deleted

    # Hand-specific blocks are shared objects, not copied into every candidate.
    plus_reference_name = "orientation_preserving_R0_linear"
    plus_candidate = candidates[plus_reference_name]
    plus_support = score_hand(
        old_full, fresh_full, old_deleted_values, fresh_deleted_values,
        plus_candidate, "plus",
        reference["candidate_maps"][plus_reference_name]["plus_rank5_support"],
        alpha, cutoff,
    )
    hand_supports = {"plus_source": plus_support}
    minus_support_ids = {}
    for candidate in manifest["candidate_maps"]:
        geometry_id = (
            f"minus_{'alexander' if candidate['alexander_reflection'] else 'orientation'}"
            f"_R{candidate['rotation_power']}"
        )
        minus_support_ids[candidate["name"]] = geometry_id
        if geometry_id not in hand_supports:
            hand_supports[geometry_id] = score_hand(
                old_full, fresh_full, old_deleted_values, fresh_deleted_values,
                candidate, "minus",
                reference["candidate_maps"][candidate["name"]]["minus_rank5_support"],
                alpha, cutoff,
            )
        results[candidate["name"]]["hand_specific_support_refs"] = {
            "plus": "plus_source",
            "minus": geometry_id,
        }

    primary = list(manifest["primary_retrospective_candidates"])
    primary_old = np.concatenate([old_influences[name] for name in primary], axis=1)
    primary_fresh = np.concatenate([fresh_influences[name] for name in primary], axis=1)
    cross_old = jackknife_covariance(primary_old)
    cross_fresh = jackknife_covariance(primary_fresh)
    cross_total = cross_old + cross_fresh
    slices = {
        name: [70 * index, 70 * (index + 1)]
        for index, name in enumerate(primary)
    }
    primary_survivors = [
        name for name in primary if results[name]["decision"] == "survives"
    ]
    secondary = [
        candidate["name"] for candidate in manifest["candidate_maps"]
        if candidate["name"] not in primary
    ]
    secondary_survivors = [name for name in secondary if results[name]["decision"] == "survives"]
    hand_complete = [
        name for name in primary
        if all(
            hand_supports[results[name]["hand_specific_support_refs"][hand]]["decision"] == "survives"
            for hand in ("plus", "minus")
        )
    ]
    if primary_survivors:
        decision = "declared_fixed_bridge_family_survives"
    elif hand_complete:
        decision = "declared_fixed_bridges_rejected_general_5plus5_still_viable"
    else:
        decision = "declared_fixed_bridges_rejected_and_full_hand_rank5_support_incomplete"
    return {
        "schema": "matching-one/z5-projective-leg-augmented-joint-score/v2",
        "status": "existing_data_old_plus_degree5_joint_operator",
        "operator": {
            "old_matrix_shape": [24, 6],
            "extension_matrix_shape": [16, 6],
            "augmented_matrix_shape": [40, 6],
            "rank_null": "rank<=5",
            "residual": "19 old plus 16 extension complex Schur coordinates, realified to 70",
            "pivot_policy": "candidate-specific pivot rows and columns frozen from 99d23a7; never reselected with fresh data or inside delete-one recomputation",
        },
        "candidate_maps": results,
        "hand_specific_support_blocks": hand_supports,
        "primary_retrospective_candidates": primary,
        "primary_retrospective_survivors": primary_survivors,
        "secondary_candidates": secondary,
        "secondary_survivors": secondary_survivors,
        "primary_candidates_with_both_hand_support": hand_complete,
        "primary_cross_candidate_covariance": {
            "candidate_order": primary,
            "candidate_slices": slices,
            "dimension": int(cross_total.shape[0]),
            "old_influence_covariance": cross_old.tolist(),
            "fresh_influence_covariance": cross_fresh.tolist(),
            "total_covariance": cross_total.tolist(),
            "scope": "saved for correlated contrasts; no p-values are multiplied or added",
        },
        "decision_alpha": alpha,
        "decision": decision,
        "claim_boundary": [
            "The fixed-map family itself predates the radius-five reveal, but this augmented joint gate was designed after the marginal R2/R3 conflict was observed; it is a primary retrospective zero-sample synthesis, not a prospective confirmation.",
            "The other eleven maps were frozen in the later joint-annihilation manifest and are reported as secondary retrospective existing-data views.",
            "A candidate is rejected by one 70-real-coordinate joint Wald gate; old and fresh p-values are never combined as independent votes.",
            "The old and fresh streams have independent seeds.  Within each stream one delete-one removes every hand, charge and displacement together; equal batch indices across streams are never paired.",
            "Hand-specific augmented support is reported because failure there cannot be interpreted as bridge-only failure.",
            "The finite-batch Hotelling calibration is the existing project approximation for summed independent jackknife covariances, not an exact two-sample pivot law.",
            "Survival means compatibility with a declared truncated rank-at-most-five chart, not exact rank five, a closed transfer algebra or a continuum field identity.",
        ],
    }


def compact_score_covariance(score_row: dict, prefix: str, arrays: dict[str, np.ndarray]) -> None:
    point_key = prefix + "__point"
    old_key = prefix + "__old_covariance"
    fresh_key = prefix + "__fresh_covariance"
    eigen_key = prefix + "__correlation_eigenvalues"
    arrays[point_key] = np.asarray(score_row.pop("residual"), dtype=float)
    arrays[old_key] = np.asarray(score_row.pop("old_influence_covariance"), dtype=float)
    arrays[fresh_key] = np.asarray(score_row.pop("fresh_influence_covariance"), dtype=float)
    # Total covariance is intentionally not stored: it is exactly old + fresh.
    score_row.pop("total_covariance")
    arrays[eigen_key] = np.asarray(score_row.pop("correlation_eigenvalues"), dtype=float)
    score_row["covariance_payload_refs"] = {
        "point": point_key,
        "old": old_key,
        "fresh": fresh_key,
        "total": f"derive `{old_key} + {fresh_key}`",
        "correlation_eigenvalues": eigen_key,
    }


def write_compact_covariance_payload(result: dict, path: Path) -> dict:
    arrays: dict[str, np.ndarray] = {}
    for name, candidate in result["candidate_maps"].items():
        compact_score_covariance(candidate["joint_score"], f"joint__{name}", arrays)
    for support_id, support in result["hand_specific_support_blocks"].items():
        compact_score_covariance(support["score"], f"hand__{support_id}", arrays)
    cross = result["primary_cross_candidate_covariance"]
    arrays["primary_cross__old_covariance"] = np.asarray(
        cross.pop("old_influence_covariance"), dtype=float,
    )
    arrays["primary_cross__fresh_covariance"] = np.asarray(
        cross.pop("fresh_influence_covariance"), dtype=float,
    )
    cross.pop("total_covariance")
    cross["covariance_payload_refs"] = {
        "old": "primary_cross__old_covariance",
        "fresh": "primary_cross__fresh_covariance",
        "total": "derive `primary_cross__old_covariance + primary_cross__fresh_covariance`",
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    payload = path.read_bytes()
    return {
        "path": path.as_posix(),
        "sha256": sha256_bytes(payload),
        "format": "numpy savez_compressed NPZ",
        "array_count": len(arrays),
        "arrays": {
            key: {"shape": list(value.shape), "dtype": str(value.dtype)}
            for key, value in arrays.items()
        },
        "derivation": "every total covariance is the exact array sum old + fresh",
        "deduplication": (
            "one canonical plus-hand block and eight minus-hand geometry blocks are stored; "
            "candidate summaries reference them rather than copying 32 hand blocks"
        ),
    }


def render(result: Mapping[str, object]) -> str:
    lines = [
        "# P250 augmented old-plus-degree-five joint operator",
        "",
        "| candidate | scope | plus p | minus p | joint p | joint df | decision |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for name, row in result["candidate_maps"].items():
        refs = row["hand_specific_support_refs"]
        plus_p = result["hand_specific_support_blocks"][refs["plus"]]["score"]["finite_batch_survival_p"]
        minus_p = result["hand_specific_support_blocks"][refs["minus"]]["score"]["finite_batch_survival_p"]
        joint = row["joint_score"]
        lines.append(
            f"| {name} | {row['scope']} | {plus_p:.6g} | {minus_p:.6g} | "
            f"{joint['finite_batch_survival_p']:.6g} | {joint['resolved_covariance_modes']} | {row['decision']} |"
        )
    lines += [
        "",
        f"Primary retrospective survivors: `{result['primary_retrospective_survivors']}`.",
        f"Secondary survivors: `{result['secondary_survivors']}`.",
        f"Primary candidates with both hand-specific rank-five gates surviving: `{result['primary_candidates_with_both_hand_support']}`.",
        f"Decision: `{result['decision']}`.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--covariance-output", type=Path, required=True)
    args = parser.parse_args()
    manifest_payload = args.manifest.read_bytes()
    manifest = json.loads(manifest_payload)

    streams = {}
    audits = {}
    for key in ("old_stream", "fresh_stream"):
        spec = manifest[key]
        payload = pinned_bytes(Path(spec["path"]), spec["commit"])
        if sha256_bytes(payload) != spec["sha256"]:
            raise ValueError(f"{key} hash changed")
        rows = read_batches(payload)
        streams[key] = rows
        audits[key] = validate_stream(rows, spec)

    reference_spec = manifest["reference_joint_score"]
    reference_payload = pinned_bytes(Path(reference_spec["path"]), reference_spec["commit"])
    if sha256_bytes(reference_payload) != reference_spec["sha256"]:
        raise ValueError("pinned 99d23a7 joint score hash changed")
    reference = json.loads(reference_payload)

    radius5_spec = manifest["reference_radius5_score"]
    radius5_payload = pinned_bytes(Path(radius5_spec["path"]), radius5_spec["commit"])
    if sha256_bytes(radius5_payload) != radius5_spec["sha256"]:
        raise ValueError("pinned 11130ae radius-five score hash changed")
    radius5_reference = json.loads(radius5_payload)
    semantic_audit = validate_input_semantics(manifest, radius5_reference)

    result = score(streams["old_stream"], streams["fresh_stream"], manifest, reference)
    result["manifest"] = {"path": args.manifest.as_posix(), "sha256": sha256_bytes(manifest_payload)}
    result["inputs"] = {
        "old_stream": {**manifest["old_stream"], "audit": audits["old_stream"]},
        "fresh_stream": {**manifest["fresh_stream"], "audit": audits["fresh_stream"]},
        "cross_stream_independence": {
            "seed_distinct": manifest["old_stream"]["seed"] != manifest["fresh_stream"]["seed"],
            "batch_indices_are_not_paired": True,
        },
        "reference_joint_score": reference_spec,
        "reference_radius5_score": {
            **radius5_spec,
            "published_surviving_maps": radius5_reference["surviving_maps"],
            "published_decision": radius5_reference["decision"],
        },
    }
    result["input_semantic_audit"] = semantic_audit
    result["runtime"] = {
        "python": sys.version.split()[0],
        "machine": platform.machine(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "command": (
            f"python3 {Path(__file__).as_posix()} {args.manifest.as_posix()} "
            f"--output {args.output.as_posix()} --markdown {args.markdown.as_posix()} "
            f"--covariance-output {args.covariance_output.as_posix()}"
        ),
        "server_used": False,
    }
    result["covariance_payload"] = write_compact_covariance_payload(
        result, args.covariance_output,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
