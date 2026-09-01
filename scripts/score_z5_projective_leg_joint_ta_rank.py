#!/usr/bin/env python3
"""Score a shared rank-2/rank-3 recurrence for the pinned P250 T/A data.

The input remains branch-only.  It is read byte-for-byte with ``git show`` at
the commit and path frozen in the manifest, then checked against its SHA256.
All eight complex T/A x hand x charge sequences share one recurrence.  The
coefficients are fit on d1..d4, while d5 remains a strict holdout.  Delete-one
jackknife replicates always remove one aligned batch from every sequence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import itertools
import json
import math
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np


SEQUENCE_SPECS = (
    ("T_plus_r1", "T", "plus", 1),
    ("T_plus_r2", "T", "plus", 2),
    ("T_minus_r1", "T", "minus", 1),
    ("T_minus_r2", "T", "minus", 2),
    ("A_plus_r1", "A", "plus", 1),
    ("A_plus_r2", "A", "plus", 2),
    ("A_minus_r1", "A", "minus", 1),
    ("A_minus_r2", "A", "minus", 2),
)


class NotScoreable(RuntimeError):
    """Raised when a frozen contract fails instead of being silently repaired."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        default="analysis/p250_joint_ta_hankel_manifest.json",
        help="Frozen analysis manifest.",
    )
    parser.add_argument(
        "--output-json",
        default="results/p250-projective-leg-joint-ta-rank/latest.json",
    )
    parser.add_argument(
        "--output-md",
        default="results/p250-projective-leg-joint-ta-rank/latest.md",
    )
    return parser.parse_args()


def git(repo: Path, *args: str, text: bool = True) -> str | bytes:
    completed = subprocess.run(
        ["git", *args],
        cwd=repo,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=text,
    )
    return completed.stdout


def load_manifest(path: Path) -> dict[str, Any]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "matching-one.p250-joint-ta-hankel-manifest.v1":
        raise NotScoreable("unexpected manifest schema")
    if manifest.get("analysis_schema") != "matching-one.z5-projective-leg-joint-ta-rank.v1":
        raise NotScoreable("unexpected analysis schema")
    if manifest.get("new_samples") is not False:
        raise NotScoreable("this scorer is restricted to existing-data reuse")
    source = manifest.get("source", {})
    for key in ("commit", "path", "sha256", "aligned_batches", "samples", "dependency_group"):
        if key not in source:
            raise NotScoreable(f"missing source.{key}")
    if manifest.get("semantics", {}).get("sequence_order") != [spec[0] for spec in SEQUENCE_SPECS]:
        raise NotScoreable("manifest sequence order does not match scorer contract")
    model_ranks = [model.get("rank") for model in manifest.get("models", [])]
    if model_ranks != [2, 3]:
        raise NotScoreable("only the frozen rank-2 then rank-3 gate is allowed")
    return manifest


def read_pinned_csv(repo: Path, source: dict[str, Any]) -> tuple[bytes, list[dict[str, str]]]:
    spec = f"{source['commit']}:{source['path']}"
    blob = git(repo, "show", spec, text=False)
    assert isinstance(blob, bytes)
    digest = hashlib.sha256(blob).hexdigest()
    if digest != source["sha256"]:
        raise NotScoreable(
            f"pinned CSV SHA256 mismatch: expected {source['sha256']}, observed {digest}"
        )
    try:
        decoded = blob.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise NotScoreable("pinned CSV is not UTF-8") from exc
    reader = csv.DictReader(io.StringIO(decoded))
    rows = list(reader)
    if not reader.fieldnames:
        raise NotScoreable("pinned CSV has no header")
    required = {"batch", "replica_first", "samples", "field_sha256", "translation_sha256"}
    for distance in range(1, 6):
        for _, gamma, hand, charge in SEQUENCE_SPECS:
            required.add(f"d{distance}_{gamma}{charge}_{hand}_re")
            required.add(f"d{distance}_{gamma}{charge}_{hand}_im")
    missing = sorted(required.difference(reader.fieldnames))
    if missing:
        raise NotScoreable(f"pinned CSV is missing required columns: {missing}")
    return blob, rows


def validate_rows(rows: Sequence[dict[str, str]], source: dict[str, Any]) -> np.ndarray:
    expected_batches = int(source["aligned_batches"])
    if len(rows) != expected_batches:
        raise NotScoreable(f"expected {expected_batches} batches, found {len(rows)}")
    batches = [int(row["batch"]) for row in rows]
    if batches != list(range(expected_batches)):
        raise NotScoreable("batch labels are not the unique ordered range 0..B-1")
    samples = np.asarray([int(row["samples"]) for row in rows], dtype=np.int64)
    if np.any(samples <= 0):
        raise NotScoreable("all aligned batches must contain positive sample counts")
    expected_per_batch = source.get("samples_per_batch")
    if expected_per_batch is not None and np.any(samples != int(expected_per_batch)):
        raise NotScoreable("batch sample counts differ from frozen samples_per_batch")
    if int(samples.sum()) != int(source["samples"]):
        raise NotScoreable("batch sample counts do not sum to frozen source.samples")
    for row in rows:
        for distance in range(1, 6):
            for _, gamma, hand, charge in SEQUENCE_SPECS:
                for part in ("re", "im"):
                    value = float(row[f"d{distance}_{gamma}{charge}_{hand}_{part}"])
                    if not math.isfinite(value):
                        raise NotScoreable("non-finite response entry in pinned CSV")
    return samples


def aggregate_sequences(
    rows: Sequence[dict[str, str]], samples: np.ndarray, omit: int | None = None
) -> dict[str, np.ndarray]:
    keep = np.ones(len(rows), dtype=bool)
    if omit is not None:
        keep[omit] = False
    denominator = int(samples[keep].sum())
    if denominator <= 0:
        raise NotScoreable("delete-one aggregation has no samples")
    sequences: dict[str, np.ndarray] = {}
    for label, gamma, hand, charge in SEQUENCE_SPECS:
        values: list[complex] = []
        for distance in range(1, 6):
            real = sum(
                float(rows[index][f"d{distance}_{gamma}{charge}_{hand}_re"])
                for index in range(len(rows))
                if keep[index]
            )
            imag = sum(
                float(rows[index][f"d{distance}_{gamma}{charge}_{hand}_im"])
                for index in range(len(rows))
                if keep[index]
            )
            values.append(complex(real, imag) / denominator)
        sequences[label] = np.asarray(values, dtype=np.complex128)
    return sequences


def fit_common_recurrence(sequences: dict[str, np.ndarray], rank: int) -> np.ndarray:
    design: list[list[complex]] = []
    targets: list[complex] = []
    # zero-based target indices rank..3 correspond to the frozen d1..d4 fit window.
    for label, *_ in SEQUENCE_SPECS:
        sequence = sequences[label]
        for target_index in range(rank, 4):
            design.append([sequence[target_index - lag] for lag in range(1, rank + 1)])
            targets.append(sequence[target_index])
    matrix = np.asarray(design, dtype=np.complex128)
    vector = np.asarray(targets, dtype=np.complex128)
    if matrix.shape[0] <= rank:
        raise NotScoreable("common recurrence is not overdetermined on d1..d4")
    coefficients, _, matrix_rank, singular_values = np.linalg.lstsq(matrix, vector, rcond=None)
    if int(matrix_rank) != rank:
        raise NotScoreable(f"rank-{rank} recurrence design is singular")
    if not np.all(np.isfinite(singular_values)) or not np.all(np.isfinite(coefficients)):
        raise NotScoreable(f"rank-{rank} recurrence fit is non-finite")
    return coefficients


def residuals_complex(
    sequences: dict[str, np.ndarray], coefficients: np.ndarray, heldout: bool
) -> np.ndarray:
    rank = len(coefficients)
    residuals: list[complex] = []
    target_indices: Iterable[int] = (4,) if heldout else range(rank, 4)
    for target_index in target_indices:
        for label, *_ in SEQUENCE_SPECS:
            sequence = sequences[label]
            prediction = sum(
                coefficients[lag - 1] * sequence[target_index - lag]
                for lag in range(1, rank + 1)
            )
            residuals.append(sequence[target_index] - prediction)
    return np.asarray(residuals, dtype=np.complex128)


def complex_vector_to_real(vector: np.ndarray) -> np.ndarray:
    result = np.empty(2 * len(vector), dtype=np.float64)
    result[0::2] = vector.real
    result[1::2] = vector.imag
    return result


def complex_coefficients_to_real(coefficients: np.ndarray) -> np.ndarray:
    return complex_vector_to_real(coefficients)


def jackknife_covariance(replicates: np.ndarray) -> np.ndarray:
    if replicates.ndim != 2 or replicates.shape[0] < 2:
        raise NotScoreable("jackknife covariance needs at least two aligned deletions")
    centered = replicates - replicates.mean(axis=0, keepdims=True)
    batch_count = replicates.shape[0]
    covariance = (batch_count - 1.0) / batch_count * centered.T @ centered
    return (covariance + covariance.T) / 2.0


def chi_square_survival(chi_square: float, degrees_of_freedom: int) -> float:
    if degrees_of_freedom <= 0 or chi_square < 0 or not math.isfinite(chi_square):
        raise NotScoreable("invalid chi-square arguments")
    x = chi_square / 2.0
    a = degrees_of_freedom / 2.0
    if degrees_of_freedom % 2 == 0:
        term = 1.0
        total = term
        for index in range(1, int(a)):
            term *= x / index
            total += term
        return min(1.0, max(0.0, math.exp(-x) * total))
    # Q(1/2,x)=erfc(sqrt(x)); step upward with
    # Q(a+1,x)=Q(a,x)+x**a exp(-x)/Gamma(a+1).
    value = math.erfc(math.sqrt(x))
    current_a = 0.5
    while current_a < a:
        log_term = current_a * math.log(x) - x - math.lgamma(current_a + 1.0) if x > 0 else -math.inf
        value += math.exp(log_term) if math.isfinite(log_term) else 0.0
        current_a += 1.0
    return min(1.0, max(0.0, value))


def covariance_score(
    residual: np.ndarray, covariance: np.ndarray, nominal_df: int
) -> dict[str, Any]:
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    largest = float(max(0.0, eigenvalues[-1]))
    cutoff = largest * 1.0e-10
    positive = eigenvalues > cutoff
    numerical_rank = int(np.count_nonzero(positive))
    if numerical_rank < nominal_df:
        return {
            "scoreable": False,
            "reason": "covariance_rank_below_nominal_degrees_of_freedom",
            "nominal_degrees_of_freedom": nominal_df,
            "covariance_numerical_rank": numerical_rank,
            "pseudoinverse_relative_cutoff": 1.0e-10,
            "pseudoinverse_absolute_cutoff": cutoff,
            "covariance_eigenvalues": [float(value) for value in eigenvalues],
        }
    kept_values = eigenvalues[positive]
    kept_vectors = eigenvectors[:, positive]
    coordinates = kept_vectors.T @ residual
    chi_square = float(np.sum((coordinates * coordinates) / kept_values))
    projection = kept_vectors @ coordinates
    residual_norm = float(np.linalg.norm(residual))
    null_norm = float(np.linalg.norm(residual - projection))
    null_fraction = null_norm / max(residual_norm, np.finfo(float).tiny)
    if null_fraction > 1.0e-6:
        return {
            "scoreable": False,
            "reason": "full_sample_residual_has_non_negligible_covariance_null_component",
            "nominal_degrees_of_freedom": nominal_df,
            "covariance_numerical_rank": numerical_rank,
            "pseudoinverse_relative_cutoff": 1.0e-10,
            "pseudoinverse_absolute_cutoff": cutoff,
            "nullspace_fraction": null_fraction,
            "covariance_eigenvalues": [float(value) for value in eigenvalues],
        }
    p_value = chi_square_survival(chi_square, nominal_df)
    condition_number = float(kept_values[-1] / kept_values[0])
    return {
        "scoreable": True,
        "chi_square": chi_square,
        "degrees_of_freedom": nominal_df,
        "p_value": p_value,
        "covariance_numerical_rank": numerical_rank,
        "covariance_condition_number": condition_number,
        "pseudoinverse_relative_cutoff": 1.0e-10,
        "pseudoinverse_absolute_cutoff": cutoff,
        "nullspace_fraction": null_fraction,
        "covariance_eigenvalues": [float(value) for value in eigenvalues],
    }


def complex_object(value: complex) -> dict[str, float]:
    return {
        "re": float(value.real),
        "im": float(value.imag),
        "abs": float(abs(value)),
        "phase": float(math.atan2(value.imag, value.real)),
    }


def vector_labels(rank: int, heldout: bool) -> list[str]:
    target_distances = [5] if heldout else list(range(rank + 1, 5))
    labels: list[str] = []
    for target in target_distances:
        for sequence, *_ in SEQUENCE_SPECS:
            labels.extend([f"d{target}:{sequence}:re", f"d{target}:{sequence}:im"])
    return labels


def align_roots(full_roots: np.ndarray, replicate_roots: np.ndarray) -> tuple[np.ndarray, float, float]:
    candidates: list[tuple[float, tuple[int, ...]]] = []
    for permutation in itertools.permutations(range(len(replicate_roots))):
        cost = sum(abs(full_roots[index] - replicate_roots[permutation[index]]) for index in range(len(full_roots)))
        candidates.append((float(cost), permutation))
    candidates.sort(key=lambda item: item[0])
    best_cost, best_permutation = candidates[0]
    second_cost = candidates[1][0] if len(candidates) > 1 else math.inf
    aligned = replicate_roots[list(best_permutation)]
    return aligned, best_cost, second_cost


def root_diagnostics(coefficients: np.ndarray, loo_coefficients: Sequence[np.ndarray]) -> dict[str, Any]:
    polynomial = np.concatenate(([1.0 + 0.0j], -coefficients))
    full_roots = np.roots(polynomial)
    aligned_replicates: list[np.ndarray] = []
    assignment_margins: list[float] = []
    for replicate_coefficients in loo_coefficients:
        replicate_roots = np.roots(np.concatenate(([1.0 + 0.0j], -replicate_coefficients)))
        aligned, best_cost, second_cost = align_roots(full_roots, replicate_roots)
        aligned_replicates.append(aligned)
        assignment_margins.append(second_cost - best_cost)
    aligned_array = np.asarray(aligned_replicates, dtype=np.complex128)
    root_summaries: list[dict[str, Any]] = []
    batch_count = aligned_array.shape[0]
    for index, root in enumerate(full_roots):
        replicate = aligned_array[:, index]
        centered = replicate - replicate.mean()
        variance = (batch_count - 1.0) / batch_count * float(np.sum(np.abs(centered) ** 2))
        phases = np.angle(replicate)
        phase_delta = np.angle(np.exp(1j * (phases - np.angle(root))))
        phase_centered = phase_delta - phase_delta.mean()
        phase_variance = (batch_count - 1.0) / batch_count * float(np.sum(phase_centered**2))
        root_summaries.append(
            {
                "root": complex_object(complex(root)),
                "jackknife_complex_se": math.sqrt(max(0.0, variance)),
                "jackknife_phase_se": math.sqrt(max(0.0, phase_variance)),
            }
        )
    minimum_margin = min(assignment_margins) if assignment_margins else math.inf
    finite = bool(np.all(np.isfinite(full_roots)) and np.all(np.isfinite(aligned_array)))
    labels_stable = finite and minimum_margin > 1.0e-10
    return {
        "interpretation": "descriptive_only" if labels_stable else "root_not_interpretable",
        "labels_stable_across_delete_one": labels_stable,
        "minimum_best_vs_second_assignment_cost_margin": minimum_margin,
        "roots": root_summaries,
    }


def coefficient_covariance(loo_coefficients: Sequence[np.ndarray]) -> np.ndarray:
    matrix = np.asarray([complex_coefficients_to_real(value) for value in loo_coefficients])
    return jackknife_covariance(matrix)


def score_rank(
    full_sequences: dict[str, np.ndarray],
    loo_sequences: Sequence[dict[str, np.ndarray]],
    rank: int,
    alpha: float,
) -> dict[str, Any]:
    coefficients = fit_common_recurrence(full_sequences, rank)
    full_training = complex_vector_to_real(residuals_complex(full_sequences, coefficients, heldout=False))
    full_holdout = complex_vector_to_real(residuals_complex(full_sequences, coefficients, heldout=True))
    loo_coefficients: list[np.ndarray] = []
    loo_training: list[np.ndarray] = []
    loo_holdout: list[np.ndarray] = []
    for sequences in loo_sequences:
        replicate_coefficients = fit_common_recurrence(sequences, rank)
        loo_coefficients.append(replicate_coefficients)
        loo_training.append(
            complex_vector_to_real(residuals_complex(sequences, replicate_coefficients, heldout=False))
        )
        loo_holdout.append(
            complex_vector_to_real(residuals_complex(sequences, replicate_coefficients, heldout=True))
        )
    training_covariance = jackknife_covariance(np.asarray(loo_training))
    holdout_covariance = jackknife_covariance(np.asarray(loo_holdout))
    coefficient_cov = coefficient_covariance(loo_coefficients)
    training_nominal_df = len(full_training) - 2 * rank
    holdout_nominal_df = len(full_holdout)
    training_score = covariance_score(full_training, training_covariance, training_nominal_df)
    holdout_score = covariance_score(full_holdout, holdout_covariance, holdout_nominal_df)
    closes = bool(
        training_score.get("scoreable")
        and holdout_score.get("scoreable")
        and training_score["p_value"] >= alpha
        and holdout_score["p_value"] >= alpha
    )
    return {
        "rank": rank,
        "recurrence": "s(d) = sum_j coefficient[lag_j] * s(d-j)",
        "fit_distances": [1, 2, 3, 4],
        "training_targets": list(range(rank + 1, 5)),
        "strict_holdout_distance": 5,
        "coefficients": [complex_object(complex(value)) for value in coefficients],
        "coefficient_vector_labels": [
            label for lag in range(1, rank + 1) for label in (f"lag_{lag}:re", f"lag_{lag}:im")
        ],
        "coefficient_jackknife_covariance": coefficient_cov.tolist(),
        "roots": root_diagnostics(coefficients, loo_coefficients),
        "training": {
            "residual_vector_labels": vector_labels(rank, heldout=False),
            "residual_vector": full_training.tolist(),
            "jackknife_covariance": training_covariance.tolist(),
            "score": training_score,
        },
        "holdout_d5": {
            "residual_vector_labels": vector_labels(rank, heldout=True),
            "residual_vector": full_holdout.tolist(),
            "jackknife_covariance": holdout_covariance.tolist(),
            "score": holdout_score,
        },
        "closes_at_alpha": closes,
    }


def format_number(value: Any, digits: int = 6) -> str:
    if value is None:
        return "n/a"
    if isinstance(value, bool):
        return "yes" if value else "no"
    return f"{float(value):.{digits}g}"


def score_cell(score: dict[str, Any]) -> str:
    if not score.get("scoreable"):
        return f"not scoreable ({score.get('reason', 'unspecified')})"
    return (
        f"chi2={format_number(score['chi_square'])}, df={score['degrees_of_freedom']}, "
        f"p={format_number(score['p_value'])}, cov-rank={score['covariance_numerical_rank']}"
    )


def render_markdown(result: dict[str, Any], manifest_path: str) -> str:
    source = result["source"]
    lines = [
        "# P250 joint T/A common-annihilator rank gate",
        "",
        f"**Decision:** `{result['decision']['outcome']}`.",
        "",
        "This is a zero-new-sample, exploratory reuse of one pinned branch-only N505 CSV. "
        "All T/A, hand, and charge coordinates share the same 160 aligned batches and the "
        f"single dependency group `{source['dependency_group']}`; they are not independent evidence.",
        "",
        "## Rank gate",
        "",
        "| shared recurrence | d1..d4 training | strict d5 holdout | closes at alpha=0.01 | roots |",
        "|---|---|---|---|---|",
    ]
    for model in result["models"]:
        lines.append(
            f"| rank {model['rank']} | {score_cell(model['training']['score'])} | "
            f"{score_cell(model['holdout_d5']['score'])} | "
            f"{'yes' if model['closes_at_alpha'] else 'no'} | "
            f"{model['roots']['interpretation']} |"
        )
    lines.extend(
        [
            "",
            "The recurrence is `s(d)=a1 s(d-1)+...+aK s(d-K)`, with one complex coefficient "
            "vector shared by all eight sequences. The fit uses only targets within d1..d4. "
            "The d5 residual is computed after the coefficients are frozen.",
            "",
            "## Coefficients and descriptive roots",
            "",
        ]
    )
    for model in result["models"]:
        coefficients = ", ".join(
            f"a{index + 1}={format_number(item['re'])}{float(item['im']):+.6g}i"
            for index, item in enumerate(model["coefficients"])
        )
        roots = ", ".join(
            f"{format_number(item['root']['re'])}{float(item['root']['im']):+.6g}i "
            f"(complex jackknife SE {format_number(item['jackknife_complex_se'])})"
            for item in model["roots"]["roots"]
        )
        lines.extend(
            [
                f"- Rank {model['rank']}: {coefficients}.",
                f"  Roots ({model['roots']['interpretation']}): {roots}.",
            ]
        )
    lines.extend(
        [
            "",
            "## Statistical contract",
            "",
            "Each delete-one replicate removes the same batch from every T/A x hand x charge "
            "coordinate, refits the common coefficients, and recomputes both training and d5 "
            "residuals. The JSON stores the complete jackknife covariance matrices. The "
            "covariance-aware quadratic forms use a fixed relative eigenvalue cutoff of 1e-10.",
            "",
            "The exact coordinate change is `X=T+A`, `Y=T-A`; therefore common recurrence "
            "closure is basis-equivalent in T/A and X/Y. T/A is used only as the stored numerical basis.",
            "",
            "## Provenance and boundaries",
            "",
            f"- Manifest: `{manifest_path}`.",
            f"- CSV: `{source['commit']}:{source['path']}`.",
            f"- CSV SHA256: `{source['sha256']}`.",
            f"- Samples/batches: {source['samples']} / {source['aligned_batches']} aligned batches.",
            f"- Analysis worktree HEAD: `{result['provenance']['analysis_worktree_head']}`.",
            "- The input remains branch-only; this analysis does not merge or copy its CSV into the Draft branch.",
            "- Rank 4 is intentionally not attempted because five distances would leave no honest holdout.",
            "- A rank-3 pass is finite-sequence closure, not a three-field, C4, or Z5 fusion theorem.",
            "- A rank-3 failure excludes the stated shared scalar rank-at-most-3 recurrence on these "
            "  five distances, but does not establish noncommuting directional transfer operators.",
            "- Root labels are descriptive; an ambiguous delete-one assignment is emitted as "
            "  `root_not_interpretable` and is never used in the rank decision. Even a matchable "
            "  label is not a precise root estimate when its reported jackknife SE is large.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    repo = Path(git(Path.cwd(), "rev-parse", "--show-toplevel").strip())
    manifest_path = Path(args.manifest)
    if not manifest_path.is_absolute():
        manifest_path = repo / manifest_path
    manifest = load_manifest(manifest_path)
    source = manifest["source"]
    blob, rows = read_pinned_csv(repo, source)
    samples = validate_rows(rows, source)
    full_sequences = aggregate_sequences(rows, samples)
    # The only resampling unit is the aligned batch.  No nested or per-channel
    # deletion is permitted by this contract.
    loo_sequences = [aggregate_sequences(rows, samples, omit=index) for index in range(len(rows))]
    alpha = float(manifest["decision"]["alpha"])
    models = [score_rank(full_sequences, loo_sequences, rank, alpha) for rank in (2, 3)]
    for model in models:
        for component in ("training", "holdout_d5"):
            if not model[component]["score"].get("scoreable"):
                reason = model[component]["score"].get("reason", "unspecified")
                raise NotScoreable(f"rank-{model['rank']} {component}: {reason}")
    rank2, rank3 = models
    if rank2["closes_at_alpha"]:
        outcome = manifest["decision"]["rank2_pass"]
    elif rank3["closes_at_alpha"]:
        outcome = manifest["decision"]["rank2_fail_rank3_pass"]
    else:
        outcome = manifest["decision"]["rank3_fail"]
    command = (
        f"python3 scripts/{Path(__file__).name} --manifest {Path(args.manifest).as_posix()} "
        f"--output-json {Path(args.output_json).as_posix()} --output-md {Path(args.output_md).as_posix()}"
    )
    result: dict[str, Any] = {
        "schema": manifest["analysis_schema"],
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "scoreable",
        "analysis_class": manifest["status"],
        "new_samples": False,
        "source": {
            **source,
            "observed_sha256": hashlib.sha256(blob).hexdigest(),
        },
        "semantics": manifest["semantics"],
        "inference": {
            "unit": "aligned_batch",
            "delete_one_replicates": len(rows),
            "channels_share_dependency": True,
            "full_covariance_saved": True,
            "pseudoinverse_relative_cutoff": 1.0e-10,
        },
        "sequence_means": {
            label: [complex_object(complex(value)) for value in full_sequences[label]]
            for label, *_ in SEQUENCE_SPECS
        },
        "models": models,
        "decision": {
            "alpha": alpha,
            "outcome": outcome,
            "rank2_closes": rank2["closes_at_alpha"],
            "rank3_closes": rank3["closes_at_alpha"],
            "chronology": "retrospective_after_related_P250_results",
            "evidence_boundary": "exploratory_existing_data_reuse_not_independent_confirmation",
        },
        "stop_boundaries": manifest["stop_boundaries"],
        "provenance": {
            "manifest": str(manifest_path.relative_to(repo)),
            "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
            "analysis_worktree_head": str(git(repo, "rev-parse", "HEAD")).strip(),
            "python": sys.version.splitlines()[0],
            "numpy": np.__version__,
            "command": command,
        },
    }
    output_json = Path(args.output_json)
    output_md = Path(args.output_md)
    if not output_json.is_absolute():
        output_json = repo / output_json
    if not output_md.is_absolute():
        output_md = repo / output_md
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    output_md.write_text(render_markdown(result, str(manifest_path.relative_to(repo))), encoding="utf-8")
    print(json.dumps({"outcome": outcome, "json": str(output_json), "markdown": str(output_md)}))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (NotScoreable, subprocess.CalledProcessError, ValueError, KeyError) as exc:
        print(f"not_scoreable: {exc}", file=sys.stderr)
        raise SystemExit(2)
