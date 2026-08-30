#!/usr/bin/env python3
"""One-shot P250 joint annihilation plus radius-five extension map score."""

from __future__ import annotations

import argparse
from itertools import combinations
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Mapping, Sequence

import numpy as np
from scipy.stats import chi2, f

from score_z5_projective_leg_bivariate_state import means as old_means, read_batches as read_old
from score_z5_projective_leg_cross_scale import jackknife_covariance
from score_z5_projective_leg_hankel_rank import MONOMIALS_2
from score_z5_projective_leg_radius5_morphism import (
    ROWS_3,
    extension_matrix,
    new_means,
    read_new,
    transformed_rows,
)
from score_z5_projective_leg_annihilator_bridge import transformed_basis


RANK = 5
EIGEN_RELATIVE_CUTOFF = 1e-10


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_batches(rows: Sequence[dict], *, samples: int, batches: int, name: str) -> dict:
    if len(rows) != batches:
        raise ValueError(f"{name}: expected {batches} batches, found {len(rows)}")
    ids = [int(row["batch"]) for row in rows]
    if ids != list(range(batches)):
        raise ValueError(f"{name}: batch ids are not contiguous")
    counts = [int(row["samples"]) for row in rows]
    if sum(counts) != samples or len(set(counts)) != 1:
        raise ValueError(f"{name}: sample total or equal-batch contract failed")
    width = counts[0]
    offsets = [int(row["replica_first"]) for row in rows]
    if offsets != [index * width for index in range(batches)]:
        raise ValueError(f"{name}: replica intervals are not the frozen domain")
    return {
        "samples": sum(counts),
        "batches": batches,
        "samples_per_batch": width,
        "replica_interval": [offsets[0], offsets[-1] + width],
    }


def candidate_matrix(
    old_values: Mapping[str, float],
    fresh_values: Mapping[str, float],
    candidate: Mapping[str, object],
) -> tuple[np.ndarray, list[dict]]:
    """Build the full 40x6 candidate map, including old and extension rows."""
    source_basis = tuple(MONOMIALS_2)
    source_rows = tuple(ROWS_3)
    target_basis = transformed_basis(
        alexander_reflection=bool(candidate["alexander_reflection"]),
        rotation_power=int(candidate["rotation_power"]),
    )
    target_rows = transformed_rows(
        bool(candidate["alexander_reflection"]), int(candidate["rotation_power"]),
    )
    plus = extension_matrix(
        old_values, fresh_values, "plus", source_basis, source_rows,
    ).conjugate()
    minus = extension_matrix(
        old_values, fresh_values, "minus", target_basis, target_rows,
    )
    labels = []
    for hand, basis, rows, conjugated in (
        ("plus", source_basis, source_rows, True),
        ("minus", target_basis, target_rows, False),
    ):
        for charge in (1, 2):
            for left in basis:
                labels.append({
                    "hand": hand,
                    "charge": charge,
                    "left_shift": list(left),
                    "stage": "radius4_annihilation",
                    "conjugated": conjugated,
                })
        for charge in (1, 2):
            for left in rows:
                labels.append({
                    "hand": hand,
                    "charge": charge,
                    "left_shift": list(left),
                    "stage": "radius5_extension",
                    "conjugated": conjugated,
                })
    return np.vstack((plus, minus)), labels


def maximum_volume_pivot(matrix: np.ndarray, rank: int = RANK) -> dict:
    """Exact global minor scan, vectorized in chunks but lexicographically stable."""
    row_sets = np.asarray(list(combinations(range(matrix.shape[0]), rank)), dtype=int)
    best_volume = -1.0
    best_rows: tuple[int, ...] | None = None
    best_columns: tuple[int, ...] | None = None
    chunk = 20000
    for columns in combinations(range(matrix.shape[1]), rank):
        column_index = np.asarray(columns, dtype=int)
        for start in range(0, len(row_sets), chunk):
            selected = row_sets[start : start + chunk]
            blocks = matrix[selected[:, :, None], column_index[None, None, :]]
            volumes = np.abs(np.linalg.det(blocks))
            local = int(np.argmax(volumes))
            volume = float(volumes[local])
            if volume > best_volume:
                best_volume = volume
                best_rows = tuple(int(value) for value in selected[local])
                best_columns = tuple(int(value) for value in columns)
    if best_rows is None or best_columns is None or not best_volume > 0.0:
        raise ValueError("no nonsingular rank-five pivot")
    pivot = matrix[np.ix_(best_rows, best_columns)]
    return {
        "rows": best_rows,
        "columns": best_columns,
        "abs_determinant": best_volume,
        "condition_number": float(np.linalg.cond(pivot)),
    }


def schur_complement(matrix: np.ndarray, pivot: Mapping[str, object]) -> np.ndarray:
    rows = tuple(int(value) for value in pivot["rows"])
    columns = tuple(int(value) for value in pivot["columns"])
    other_rows = tuple(index for index in range(matrix.shape[0]) if index not in rows)
    other_columns = tuple(index for index in range(matrix.shape[1]) if index not in columns)
    p = matrix[np.ix_(rows, columns)]
    q = matrix[np.ix_(rows, other_columns)]
    r = matrix[np.ix_(other_rows, columns)]
    s = matrix[np.ix_(other_rows, other_columns)]
    return s - r @ np.linalg.solve(p, q)


def projective_q(matrix: np.ndarray, pivot: Mapping[str, object]) -> np.ndarray:
    rows = tuple(int(value) for value in pivot["rows"])
    columns = tuple(int(value) for value in pivot["columns"])
    other = tuple(index for index in range(matrix.shape[1]) if index not in columns)
    if len(other) != 1:
        raise ValueError("rank-five chart must leave one normalization column")
    output = np.zeros(matrix.shape[1], dtype=complex)
    output[other[0]] = 1.0
    output[list(columns)] = -np.linalg.solve(
        matrix[np.ix_(rows, columns)], matrix[np.ix_(rows, other)].reshape(-1),
    )
    return output


def realify(values: np.ndarray) -> np.ndarray:
    flat = values.ravel()
    output = np.empty(2 * len(flat), dtype=float)
    output[0::2] = flat.real
    output[1::2] = flat.imag
    return output


def centered_influences(deleted: Sequence[Sequence[float]]) -> np.ndarray:
    values = np.asarray(deleted, dtype=float)
    centered = values - values.mean(axis=0)
    return math.sqrt((len(values) - 1) / len(values)) * centered


def score_from_influences(
    point: Sequence[float], old_influence: np.ndarray, fresh_influence: np.ndarray,
    *, effective_batches: int,
) -> dict:
    old_covariance = old_influence.T @ old_influence
    fresh_covariance = fresh_influence.T @ fresh_influence
    covariance = old_covariance + fresh_covariance
    point_array = np.asarray(point, dtype=float)
    scales = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    if np.any(scales <= 0.0):
        raise ValueError("augmented residual has a zero covariance scale")
    correlation = covariance / scales[:, None] / scales[None, :]
    correlation = 0.5 * (correlation + correlation.T)
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    kept = eigenvalues > EIGEN_RELATIVE_CUTOFF * eigenvalues[-1]
    normalized = point_array / scales
    projection = eigenvectors[:, kept].T @ normalized
    statistic = float(np.sum(projection**2 / eigenvalues[kept]))
    degrees = int(np.count_nonzero(kept))
    if degrees >= effective_batches:
        hotelling = math.inf
        probability = 0.0
    else:
        hotelling = (
            (effective_batches - degrees) * statistic
            / (degrees * (effective_batches - 1))
        )
        probability = float(f.sf(hotelling, degrees, effective_batches - degrees))
    return {
        "residual": point_array.tolist(),
        "old_covariance": old_covariance.tolist(),
        "fresh_covariance": fresh_covariance.tolist(),
        "covariance": covariance.tolist(),
        "covariance_addition_max_abs_error": float(np.max(np.abs(
            covariance - old_covariance - fresh_covariance
        ))),
        "eigen_relative_cutoff": EIGEN_RELATIVE_CUTOFF,
        "resolved_covariance_modes": degrees,
        "discarded_covariance_modes": int(len(point_array) - degrees),
        "asymptotic_chi_square": statistic,
        "asymptotic_degrees_of_freedom": degrees,
        "asymptotic_survival_p": float(chi2.sf(statistic, degrees)),
        "finite_batch_hotelling_F": hotelling,
        "finite_batch_numerator_df": degrees,
        "finite_batch_denominator_df": effective_batches - degrees,
        "finite_batch_survival_p": probability,
        "effective_batches_conservative": effective_batches,
    }


def git_head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"], check=True, text=True,
        stdout=subprocess.PIPE,
    ).stdout.strip()


def compute(
    old_rows: Sequence[dict], fresh_rows: Sequence[dict], manifest: Mapping[str, object],
) -> tuple[dict, dict[str, np.ndarray]]:
    old_full = old_means(old_rows)
    fresh_full = new_means(fresh_rows)
    old_deleted_values = [old_means(old_rows, index) for index in range(len(old_rows))]
    fresh_deleted_values = [new_means(fresh_rows, index) for index in range(len(fresh_rows))]
    alpha = float(manifest["decision_alpha"])
    effective_batches = min(len(old_rows), len(fresh_rows))

    results = {}
    points = []
    old_influences = []
    fresh_influences = []
    for candidate in manifest["candidate_maps"]:
        name = str(candidate["name"])
        matrix, labels = candidate_matrix(old_full, fresh_full, candidate)
        pivot = maximum_volume_pivot(matrix)
        point = realify(schur_complement(matrix, pivot))
        old_deleted = [
            realify(schur_complement(candidate_matrix(values, fresh_full, candidate)[0], pivot))
            for values in old_deleted_values
        ]
        fresh_deleted = [
            realify(schur_complement(candidate_matrix(old_full, values, candidate)[0], pivot))
            for values in fresh_deleted_values
        ]
        old_influence = centered_influences(old_deleted)
        fresh_influence = centered_influences(fresh_deleted)
        statistical = score_from_influences(
            point, old_influence, fresh_influence, effective_batches=effective_batches,
        )
        q = projective_q(matrix, pivot)
        residual_rows = [index for index in range(matrix.shape[0]) if index not in pivot["rows"]]
        residual_labels = []
        for index in residual_rows:
            prefix = labels[index]
            residual_labels.extend((
                {**prefix, "component": "real"},
                {**prefix, "component": "imaginary"},
            ))
        results[name] = {
            "candidate": dict(candidate),
            "null": "rank(augmented_candidate_matrix)<=5",
            "matrix_shape_complex": list(matrix.shape),
            "pivot": {
                "rows": list(pivot["rows"]),
                "columns": list(pivot["columns"]),
                "abs_determinant": pivot["abs_determinant"],
                "condition_number": pivot["condition_number"],
            },
            "Schur_shape_complex": [matrix.shape[0] - RANK, matrix.shape[1] - RANK],
            "projective_q": [
                {"re": float(value.real), "im": float(value.imag), "abs": float(abs(value))}
                for value in q
            ],
            "full_matrix_residual_max_abs": float(np.max(np.abs(matrix @ q))),
            "residual_order": residual_labels,
            "score": statistical,
            "decision": (
                "survives" if statistical["finite_batch_survival_p"] >= alpha else "rejected"
            ),
        }
        points.append(point)
        old_influences.append(old_influence)
        fresh_influences.append(fresh_influence)

    candidate_names = [str(row["name"]) for row in manifest["candidate_maps"]]
    survivors = [name for name in candidate_names if results[name]["decision"] == "survives"]
    alexander_survivors = [name for name in survivors if name.startswith("Alexander_")]
    old_stack = np.concatenate(old_influences, axis=1)
    fresh_stack = np.concatenate(fresh_influences, axis=1)
    cross_covariance = old_stack.T @ old_stack + fresh_stack.T @ fresh_stack
    result = {
        "schema": "matching-one/z5-projective-leg-augmented-map-score/v1",
        "status": "existing_data_augmented_candidate_map_scored_once",
        "freeze_commit": "2bca045",
        "scorer_commit": git_head(),
        "candidate_order": candidate_names,
        "candidate_maps": results,
        "surviving_maps": survivors,
        "alexander_family_survivors": alexander_survivors,
        "decision_alpha": alpha,
        "decision": (
            "fixed_parameter_free_map_survives"
            if survivors else "all_frozen_parameter_free_maps_rejected"
        ),
        "alexander_union_decision": (
            "survives" if alexander_survivors else "rejected"
        ),
        "cross_candidate_covariance": {
            "candidate_order": candidate_names,
            "candidate_slices": {
                name: [70 * index, 70 * (index + 1)]
                for index, name in enumerate(candidate_names)
            },
            "old_covariance": (old_stack.T @ old_stack).tolist(),
            "fresh_covariance": (fresh_stack.T @ fresh_stack).tolist(),
            "covariance": cross_covariance.tolist(),
            "scope": "audit only; no combined candidate vote or largest-p selection",
        },
        "deduplication": {
            "one_gate_per_candidate": True,
            "published_radius4_and_radius5_p_values_recounted": False,
            "R2_R3_old_and_fresh_gates_treated_as_independent": False,
        },
        "claim_boundary": list(manifest["claim_boundary"]),
    }
    influence_payload = {
        "candidate_names": np.asarray(candidate_names),
        "point_residuals": np.asarray(points),
        "old_centered_jackknife_influences": np.asarray(old_influences),
        "fresh_centered_jackknife_influences": np.asarray(fresh_influences),
    }
    return result, influence_payload


def render(result: Mapping[str, object]) -> str:
    lines = [
        "# P250 augmented annihilation-extension map score",
        "",
        "| frozen candidate | chi-square / df | finite-batch p | decision |",
        "|---|---:|---:|---|",
    ]
    for name in result["candidate_order"]:
        row = result["candidate_maps"][name]
        score = row["score"]
        lines.append(
            f"| {name} | {score['asymptotic_chi_square']:.6g} / "
            f"{score['resolved_covariance_modes']} | "
            f"{score['finite_batch_survival_p']:.6g} | {row['decision']} |"
        )
    lines += [
        "",
        f"Survivors: `{result['surviving_maps']}`.",
        f"Alexander union: `{result['alexander_union_decision']}`.",
        f"Decision: `{result['decision']}`.",
        "",
        "Each map is counted once in the augmented operator; earlier radius-four and radius-five p-values are not separate votes.",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--influences", type=Path, required=True)
    args = parser.parse_args()
    for target in (args.output, args.markdown, args.influences):
        if target.exists():
            raise ValueError(f"refusing to overwrite one-shot artifact: {target}")
    manifest = json.loads(args.manifest.read_text())
    old_spec = manifest["old_batches"]
    fresh_spec = manifest["fresh_batches"]
    old_path = Path(old_spec["path"])
    fresh_path = Path(fresh_spec["path"])
    if sha256(old_path) != old_spec["sha256"]:
        raise ValueError("old batch hash changed")
    if sha256(fresh_path) != fresh_spec["sha256"]:
        raise ValueError("fresh batch hash changed")
    old_rows = read_old(old_path)
    fresh_rows = read_new(fresh_path)
    audits = {
        "old": audit_batches(
            old_rows, samples=int(old_spec["samples"]), batches=int(old_spec["batches"]), name="old",
        ),
        "fresh": audit_batches(
            fresh_rows, samples=int(fresh_spec["samples"]), batches=int(fresh_spec["batches"]), name="fresh",
        ),
    }
    result, influences = compute(old_rows, fresh_rows, manifest)
    args.influences.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.influences, **influences)
    result["inputs"] = {
        "old": {**old_spec, "audit": audits["old"]},
        "fresh": {**fresh_spec, "audit": audits["fresh"]},
    }
    result["influences"] = {
        "path": args.influences.as_posix(),
        "sha256": sha256(args.influences),
        "old_shape": list(influences["old_centered_jackknife_influences"].shape),
        "fresh_shape": list(influences["fresh_centered_jackknife_influences"].shape),
        "contract": "For each source, influence.T @ influence is its covariance contribution.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
