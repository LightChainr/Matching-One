#!/usr/bin/env python3
"""Score candidate-constrained P250 two-hand joint annihilation nulls."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
from itertools import combinations
import json
import math
from pathlib import Path
import subprocess
from typing import Mapping, Sequence

import numpy as np
from scipy.stats import chi2, f


MONOMIALS_2 = ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2))
RANK = 5
EIGEN_RELATIVE_CUTOFF = 1e-10


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def pinned_bytes(path: Path, commit: str) -> bytes:
    """Read a working-tree artifact or its immutable Git object fallback."""
    if path.is_file():
        return path.read_bytes()
    completed = subprocess.run(
        ["git", "show", f"{commit}:{path.as_posix()}"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return completed.stdout


def read_batches(payload: bytes) -> list[dict]:
    text = payload.decode("utf-8")
    return [
        {
            key: int(value) if key in {"batch", "replica_first", "samples"}
            else value if key in {"field_sha256", "translation_sha256"}
            else float(value)
            for key, value in row.items()
        }
        for row in csv.DictReader(io.StringIO(text))
    ]


def validate_batches(rows: Sequence[dict], manifest: Mapping[str, object]) -> dict:
    expected_batches = int(manifest["batches"])
    expected_samples = int(manifest["samples"])
    if len(rows) != expected_batches:
        raise ValueError(f"expected {expected_batches} batches, found {len(rows)}")
    observed_ids = [int(row["batch"]) for row in rows]
    if observed_ids != list(range(expected_batches)):
        raise ValueError("batch ids are not the frozen contiguous order")
    observed_samples = sum(int(row["samples"]) for row in rows)
    if observed_samples != expected_samples:
        raise ValueError(f"expected {expected_samples} samples, found {observed_samples}")
    samples_per_batch = {int(row["samples"]) for row in rows}
    if len(samples_per_batch) != 1:
        raise ValueError("joint-annihilation scorer requires equal batch sizes")
    per_batch = samples_per_batch.pop()
    for index, row in enumerate(rows):
        if int(row["replica_first"]) != index * per_batch:
            raise ValueError("replica intervals do not match the frozen aligned batches")
    return {
        "batch_ids": [observed_ids[0], observed_ids[-1]],
        "batches": len(rows),
        "samples": observed_samples,
        "samples_per_batch": per_batch,
        "aligned_delete_one_contract": "one row deletes plus, minus, both charges, and every displacement together",
    }


def means(rows: Sequence[dict], excluded: int | None = None) -> dict[str, float]:
    kept = [row for index, row in enumerate(rows) if index != excluded]
    samples = sum(row["samples"] for row in kept)
    fields = [key for key in kept[0] if key.startswith(("ap", "am"))]
    return {key: sum(row[key] for row in kept) / samples for key in fields}


def label(a: int, b: int) -> str:
    return f"a{'p' if a >= 0 else 'm'}{abs(a)}_b{'p' if b >= 0 else 'm'}{abs(b)}"


def pair(values: Mapping[str, float], point: tuple[int, int], hand: str, charge: int) -> complex:
    prefix = f"{label(*point)}_r{charge}_{hand}_"
    return complex(values[prefix + "re"], values[prefix + "im"])


def reflect(point: tuple[int, int]) -> tuple[int, int]:
    return point[0], -point[1]


def rotate(point: tuple[int, int]) -> tuple[int, int]:
    return -point[1], point[0]


def rotate_power(point: tuple[int, int], power: int) -> tuple[int, int]:
    for _ in range(power % 4):
        point = rotate(point)
    return point


def transformed_basis(*, alexander_reflection: bool, rotation_power: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        rotate_power(reflect(point) if alexander_reflection else point, rotation_power)
        for point in MONOMIALS_2
    )


def hand_hankel(values: Mapping[str, float], hand: str, basis: Sequence[tuple[int, int]]) -> np.ndarray:
    return np.asarray([
        [pair(values, (left[0] + right[0], left[1] + right[1]), hand, charge) for right in basis]
        for charge in (1, 2) for left in basis
    ], dtype=complex)


def hand_row_labels(hand: str, basis: Sequence[tuple[int, int]], *, conjugated: bool) -> list[dict]:
    return [
        {
            "hand": hand,
            "charge": charge,
            "left_monomial": list(left),
            "moment_conjugated": conjugated,
        }
        for charge in (1, 2) for left in basis
    ]


def joint_stack(plus: np.ndarray, minus: np.ndarray, *, coefficient_conjugation: bool) -> np.ndarray:
    """Eliminate q projectively while preserving a conjugating bridge exactly.

    If q_minus=conj(q_plus), set r=q_minus.  Then H_plus q_plus=0 is
    equivalent to conj(H_plus) r=0, so both equations remain one complex
    rank-five null rather than a realified direction comparison.
    """
    source = plus.conjugate() if coefficient_conjugation else plus
    return np.vstack((source, minus))


def candidate_matrix(values: Mapping[str, float], candidate: Mapping[str, object]) -> tuple[np.ndarray, list[dict]]:
    basis0 = MONOMIALS_2
    target_basis = transformed_basis(
        alexander_reflection=bool(candidate["alexander_reflection"]),
        rotation_power=int(candidate["rotation_power"]),
    )
    conjugate = bool(candidate["coefficient_conjugation"])
    plus = hand_hankel(values, "plus", basis0)
    minus = hand_hankel(values, "minus", target_basis)
    matrix = joint_stack(plus, minus, coefficient_conjugation=conjugate)
    labels = (
        hand_row_labels("plus", basis0, conjugated=conjugate)
        + hand_row_labels("minus", target_basis, conjugated=False)
    )
    return matrix, labels


def maximum_volume_pivot(matrix: np.ndarray, rank: int) -> dict:
    best: tuple[float, tuple[int, ...], tuple[int, ...]] | None = None
    for columns in combinations(range(matrix.shape[1]), rank):
        for rows in combinations(range(matrix.shape[0]), rank):
            volume = float(abs(np.linalg.det(matrix[np.ix_(rows, columns)])))
            candidate = (volume, tuple(rows), tuple(columns))
            if best is None or candidate[0] > best[0]:
                best = candidate
    if best is None or not best[0] > 0.0:
        raise ValueError("no invertible rank-five pivot chart")
    pivot = matrix[np.ix_(best[1], best[2])]
    return {
        "rows": best[1],
        "columns": best[2],
        "abs_determinant": best[0],
        "condition_number": float(np.linalg.cond(pivot)),
    }


def schur_complement(matrix: np.ndarray, pivot: Mapping[str, object]) -> np.ndarray:
    rows = tuple(pivot["rows"])
    columns = tuple(pivot["columns"])
    other_rows = tuple(index for index in range(matrix.shape[0]) if index not in rows)
    other_columns = tuple(index for index in range(matrix.shape[1]) if index not in columns)
    p = matrix[np.ix_(rows, columns)]
    q = matrix[np.ix_(rows, other_columns)]
    r = matrix[np.ix_(other_rows, columns)]
    s = matrix[np.ix_(other_rows, other_columns)]
    return s - r @ np.linalg.solve(p, q)


def projective_q(matrix: np.ndarray, pivot: Mapping[str, object]) -> np.ndarray:
    rows = tuple(pivot["rows"])
    columns = tuple(pivot["columns"])
    other_columns = tuple(index for index in range(matrix.shape[1]) if index not in columns)
    if len(other_columns) != 1:
        raise ValueError("rank-five six-column chart must have one normalization column")
    normalized_column = other_columns[0]
    p = matrix[np.ix_(rows, columns)]
    q_column = matrix[np.ix_(rows, (normalized_column,))].reshape(-1)
    coefficients = np.zeros(matrix.shape[1], dtype=complex)
    coefficients[normalized_column] = 1.0
    coefficients[list(columns)] = -np.linalg.solve(p, q_column)
    return coefficients


def complex_payload(values: Sequence[complex]) -> list[dict]:
    return [
        {"re": float(value.real), "im": float(value.imag), "abs": float(abs(value))}
        for value in values
    ]


def realify(matrix: np.ndarray) -> list[float]:
    output = []
    for value in matrix.ravel():
        output.extend((float(value.real), float(value.imag)))
    return output


def jackknife_covariance(rows: Sequence[Sequence[float]]) -> np.ndarray:
    values = np.asarray(rows, dtype=float)
    center = values.mean(axis=0)
    factor = (len(values) - 1) / len(values)
    centered = values - center
    return factor * centered.T @ centered


def covariance_score(point: Sequence[float], deleted: Sequence[Sequence[float]]) -> dict:
    covariance = jackknife_covariance(deleted)
    point_array = np.asarray(point, dtype=float)
    scales = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    if np.any(scales <= 0.0):
        raise ValueError("joint-null residual has a zero covariance diagonal")
    correlation = covariance / scales[:, None] / scales[None, :]
    correlation = 0.5 * (correlation + correlation.T)
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    largest = float(eigenvalues[-1])
    kept = eigenvalues > EIGEN_RELATIVE_CUTOFF * largest
    if not np.any(kept):
        raise ValueError("joint-null covariance has no resolved mode")
    normalized = point_array / scales
    projections = eigenvectors[:, kept].T @ normalized
    statistic = float(np.sum(projections**2 / eigenvalues[kept]))
    degrees = int(np.count_nonzero(kept))
    batches = len(deleted)
    if degrees >= batches:
        hotelling_f = math.inf
        hotelling_p = 0.0
        denominator_df = batches - degrees
    else:
        hotelling_f = (batches - degrees) * statistic / (degrees * (batches - 1))
        denominator_df = batches - degrees
        hotelling_p = float(f.sf(hotelling_f, degrees, denominator_df))
    return {
        "residual": list(point),
        "covariance": covariance.tolist(),
        "correlation_eigenvalues": eigenvalues.tolist(),
        "eigen_relative_cutoff": EIGEN_RELATIVE_CUTOFF,
        "resolved_covariance_modes": degrees,
        "discarded_covariance_modes": len(point) - degrees,
        "asymptotic_chi_square": statistic,
        "asymptotic_degrees_of_freedom": degrees,
        "asymptotic_survival_p": float(chi2.sf(statistic, degrees)),
        "finite_batch_hotelling_F": hotelling_f,
        "finite_batch_numerator_df": degrees,
        "finite_batch_denominator_df": denominator_df,
        "finite_batch_survival_p": hotelling_p,
    }


def residual_labels(row_labels: Sequence[dict], pivot: Mapping[str, object]) -> list[str]:
    pivot_rows = set(pivot["rows"])
    labels = []
    for index, row in enumerate(row_labels):
        if index in pivot_rows:
            continue
        a, b = row["left_monomial"]
        prefix = (
            f"{row['hand']}:r{row['charge']}:u({a},{b})"
            + (":conjugated" if row["moment_conjugated"] else "")
        )
        labels.extend((prefix + ":re", prefix + ":im"))
    return labels


def rank_five_score(
    matrix: np.ndarray,
    deleted_matrices: Sequence[np.ndarray],
    row_labels: Sequence[dict],
) -> tuple[dict, list[list[float]]]:
    pivot = maximum_volume_pivot(matrix, RANK)
    point = realify(schur_complement(matrix, pivot))
    deleted = [realify(schur_complement(row, pivot)) for row in deleted_matrices]
    score = covariance_score(point, deleted)
    coefficients = projective_q(matrix, pivot)
    other_columns = [index for index in range(matrix.shape[1]) if index not in pivot["columns"]]
    result = {
        "null": "rank(candidate_mapped_joint_Hankel)<=5",
        "matrix_shape": list(matrix.shape),
        "pivot": {
            "rows": list(pivot["rows"]),
            "columns": list(pivot["columns"]),
            "abs_determinant": pivot["abs_determinant"],
            "condition_number": pivot["condition_number"],
        },
        "normalization_column": other_columns[0],
        "joint_projective_q": complex_payload(coefficients),
        "full_matrix_residual_max_abs": float(np.max(np.abs(matrix @ coefficients))),
        "singular_values_descriptive_only": [float(value) for value in np.linalg.svd(matrix, compute_uv=False)],
        "Schur_shape": [matrix.shape[0] - RANK, matrix.shape[1] - RANK],
        "residual_order": residual_labels(row_labels, pivot),
        "score": score,
    }
    return result, deleted


def hand_support_score(
    matrix: np.ndarray,
    deleted_matrices: Sequence[np.ndarray],
    hand: str,
    basis: Sequence[tuple[int, int]],
    *,
    conjugated: bool,
) -> dict:
    labels = hand_row_labels(hand, basis, conjugated=conjugated)
    result, _ = rank_five_score(matrix, deleted_matrices, labels)
    return result


def compare_identity_replay(observed: Mapping[str, object], reference: Mapping[str, object]) -> dict:
    expected = reference["groups"]["shared_block"]["rank_nulls"]["5"]
    observed_score = observed["score"]
    expected_score = expected["score"]
    residual_difference = float(np.max(np.abs(
        np.asarray(observed_score["residual"]) - np.asarray(expected_score["residual"])
    )))
    covariance_difference = float(np.max(np.abs(
        np.asarray(observed_score["covariance"]) - np.asarray(expected_score["covariance"])
    )))
    p_difference = abs(
        float(observed_score["finite_batch_survival_p"])
        - float(expected_score["finite_batch_survival_p"])
    )
    statistic_difference = abs(
        float(observed_score["asymptotic_chi_square"])
        - float(expected_score["asymptotic_chi_square"])
    )
    pivot_matches = observed["pivot"] == expected["pivot"]
    passed = (
        pivot_matches
        and residual_difference <= 1e-14
        and covariance_difference <= 1e-14
        and p_difference <= 1e-15
        and statistic_difference <= 1e-13
        and int(observed_score["resolved_covariance_modes"])
        == int(expected_score["resolved_covariance_modes"])
    )
    return {
        "passed": passed,
        "reference_null": expected["null"],
        "pivot_exact_match": pivot_matches,
        "residual_max_abs_difference": residual_difference,
        "covariance_max_abs_difference": covariance_difference,
        "finite_batch_p_abs_difference": p_difference,
        "asymptotic_chi_square_abs_difference": statistic_difference,
        "expected_finite_batch_p": expected_score["finite_batch_survival_p"],
        "observed_finite_batch_p": observed_score["finite_batch_survival_p"],
        "expected_resolved_modes": expected_score["resolved_covariance_modes"],
        "observed_resolved_modes": observed_score["resolved_covariance_modes"],
    }


def candidate_payload(
    values: Mapping[str, float],
    deleted_values: Sequence[Mapping[str, float]],
    candidate: Mapping[str, object],
    alpha: float,
) -> tuple[dict, list[list[float]]]:
    matrix, labels = candidate_matrix(values, candidate)
    deleted_matrices = [candidate_matrix(row, candidate)[0] for row in deleted_values]
    joint, deleted_joint = rank_five_score(matrix, deleted_matrices, labels)

    basis0 = MONOMIALS_2
    target_basis = transformed_basis(
        alexander_reflection=bool(candidate["alexander_reflection"]),
        rotation_power=int(candidate["rotation_power"]),
    )
    conjugate = bool(candidate["coefficient_conjugation"])
    plus = hand_hankel(values, "plus", basis0)
    deleted_plus = [hand_hankel(row, "plus", basis0) for row in deleted_values]
    if conjugate:
        plus = plus.conjugate()
        deleted_plus = [row.conjugate() for row in deleted_plus]
    minus = hand_hankel(values, "minus", target_basis)
    deleted_minus = [hand_hankel(row, "minus", target_basis) for row in deleted_values]
    plus_support = hand_support_score(
        plus, deleted_plus, "plus", basis0, conjugated=conjugate,
    )
    minus_support = hand_support_score(
        minus, deleted_minus, "minus", target_basis, conjugated=False,
    )

    joint_p = float(joint["score"]["finite_batch_survival_p"])
    output = {
        "alexander_reflection": bool(candidate["alexander_reflection"]),
        "rotation_power": int(candidate["rotation_power"]),
        "coefficient_conjugation": conjugate,
        "role": candidate["role"],
        "target_basis": [list(point) for point in target_basis],
        "plus_rank5_support": plus_support,
        "minus_rank5_support": minus_support,
        "joint_rank5": joint,
        "decision_alpha": alpha,
        "decision": "survives" if joint_p >= alpha else "rejected",
    }
    return output, deleted_joint


def score(rows: Sequence[dict], manifest: Mapping[str, object], reference: Mapping[str, object]) -> dict:
    alpha = float(manifest["decision_alpha"])
    values = means(rows)
    deleted_values = [means(rows, index) for index in range(len(rows))]
    candidates = {row["name"]: row for row in manifest["candidate_maps"]}
    identity_name = manifest["identity_linear_replay"]["candidate"]
    if identity_name not in candidates:
        raise ValueError("identity-linear replay candidate is absent")

    results: dict[str, dict] = {}
    deleted_by_candidate: dict[str, list[list[float]]] = {}
    identity, identity_deleted = candidate_payload(
        values, deleted_values, candidates[identity_name], alpha,
    )
    replay = compare_identity_replay(identity["joint_rank5"], reference)
    if not replay["passed"]:
        raise ValueError(f"identity-linear replay failed: {replay}")
    results[identity_name] = identity
    deleted_by_candidate[identity_name] = identity_deleted

    for name, candidate in candidates.items():
        if name == identity_name:
            continue
        result, deleted = candidate_payload(values, deleted_values, candidate, alpha)
        results[name] = result
        deleted_by_candidate[name] = deleted

    core_order = list(manifest["core_cross_candidate_covariance"])
    core_slices = {}
    cursor = 0
    for name in core_order:
        width = len(results[name]["joint_rank5"]["score"]["residual"])
        core_slices[name] = [cursor, cursor + width]
        cursor += width
    core_deleted = [
        np.concatenate([deleted_by_candidate[name][batch] for name in core_order]).tolist()
        for batch in range(len(rows))
    ]
    core_covariance = jackknife_covariance(core_deleted)

    post_radius5_primary = manifest["post_radius5_primary_candidate"]
    primary_p = results[post_radius5_primary]["joint_rank5"]["score"]["finite_batch_survival_p"]
    primary_family = list(manifest["primary_family"])
    family_survivors = [
        name for name in primary_family
        if results[name]["joint_rank5"]["score"]["finite_batch_survival_p"] >= alpha
    ]
    all_survivors = [
        name for name in candidates
        if results[name]["joint_rank5"]["score"]["finite_batch_survival_p"] >= alpha
    ]
    return {
        "schema": "matching-one/z5-projective-leg-joint-annihilation-score/v1",
        "status": "existing_data_candidate_constrained_joint_null",
        "identity_linear_replay": replay,
        "candidate_maps": results,
        "post_radius5_primary_candidate": post_radius5_primary,
        "post_radius5_primary_decision": (
            "joint_null_survives" if primary_p >= alpha else "joint_null_rejected"
        ),
        "primary_family": primary_family,
        "primary_family_survivors": family_survivors,
        "all_candidate_survivors": all_survivors,
        "core_cross_candidate_covariance": {
            "candidate_order": core_order,
            "candidate_slices": core_slices,
            "dimension": int(core_covariance.shape[0]),
            "covariance": core_covariance.tolist(),
            "scope": "saved for correlated candidate contrasts; no combined vote or p-value is formed",
        },
        "decision_alpha": alpha,
        "decision": (
            "radius5_selected_R2_joint_null_survives"
            if primary_p >= alpha
            else "radius5_selected_R2_joint_null_rejected"
        ),
        "claim_boundary": [
            "This score tests a degree-two, radius-four candidate-mapped common annihilator and not exact rank five, a closed transfer algebra, ordered TxTy/TyTx words, or a continuum field identity.",
            "All candidate scores reuse the same 80k/400-batch archive and are correlated views of one dependency group.",
            "The later radius-five result uses an independent fresh stream but its published map score also reuses this old archive, so its p-value is not independent of this score.",
            "Survival retains a declared truncated map; it does not select a map by largest p-value.  The primary family is rejected only if every declared member is rejected.",
            "The finite-batch Hotelling calibration follows the existing nonlinear Schur-jackknife project contract and is not an exact finite-replica pivot conditional law.",
        ],
    }


def render(result: Mapping[str, object]) -> str:
    lines = [
        "# P250 candidate-constrained joint annihilation score",
        "",
        "| candidate | plus rank<=5 p | minus rank<=5 p | joint rank<=5 p | joint df | decision |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for name, row in result["candidate_maps"].items():
        plus_p = row["plus_rank5_support"]["score"]["finite_batch_survival_p"]
        minus_p = row["minus_rank5_support"]["score"]["finite_batch_survival_p"]
        joint = row["joint_rank5"]["score"]
        lines.append(
            f"| {name} | {plus_p:.6g} | {minus_p:.6g} | "
            f"{joint['finite_batch_survival_p']:.6g} | {joint['resolved_covariance_modes']} | {row['decision']} |"
        )
    replay = result["identity_linear_replay"]
    lines += [
        "",
        f"Identity-linear replay exact: `{replay['passed']}`; p difference `{replay['finite_batch_p_abs_difference']:.3g}`.",
        f"Post-radius-five primary: `{result['post_radius5_primary_candidate']}` -> `{result['post_radius5_primary_decision']}`.",
        f"Primary-family survivors: `{result['primary_family_survivors']}`.",
        f"Decision: `{result['decision']}`.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    manifest_payload = args.manifest.read_bytes()
    manifest = json.loads(manifest_payload)

    input_spec = manifest["input_batches"]
    input_payload = pinned_bytes(Path(input_spec["path"]), input_spec["commit"])
    if sha256_bytes(input_payload) != input_spec["sha256"]:
        raise ValueError("P250 bivariate batch hash changed")
    rows = read_batches(input_payload)
    batch_audit = validate_batches(rows, manifest)

    reference_spec = manifest["identity_linear_replay"]["reference_score"]
    reference_payload = pinned_bytes(Path(reference_spec["path"]), reference_spec["commit"])
    if sha256_bytes(reference_payload) != reference_spec["sha256"]:
        raise ValueError("a770ac9 reference score hash changed")
    reference = json.loads(reference_payload)

    result = score(rows, manifest, reference)
    result["manifest"] = {
        "path": args.manifest.as_posix(),
        "sha256": sha256_bytes(manifest_payload),
    }
    result["input_batches"] = {**input_spec, "batch_audit": batch_audit}
    result["reference_score"] = reference_spec
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
