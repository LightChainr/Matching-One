#!/usr/bin/env python3
"""Score three typed observer models on the frozen N112 rho-child batches.

The inputs are read directly with ``git show`` from one pinned commit.  This
script reconstructs the nine batch-mean coordinates and their complete
covariance of the mean.  It generates no random samples and does not change
the archived frozen score.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import platform
import subprocess
import time
from pathlib import Path
from typing import Any, Mapping

import mpmath as mp
import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "analysis" / "p275_rho_observer_models_manifest.yaml"
DEFAULT_JSON = ROOT / "results" / "p275-rho-observer-models" / "latest.json"
DEFAULT_MD = ROOT / "results" / "p275-rho-observer-models" / "latest.md"
MANIFEST_SCHEMA = "matching-one.p275-rho-observer-models.manifest.v1"
OUTPUT_SCHEMA = "matching-one.p275-rho-observer-models.v1"


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def finite_float(value: float) -> float:
    output = float(value)
    if not math.isfinite(output):
        raise ValueError("non-finite result")
    return output


def vector_payload(vector: np.ndarray) -> list[float]:
    return [finite_float(value) for value in vector]


def matrix_payload(matrix: np.ndarray) -> list[list[float]]:
    return [vector_payload(row) for row in matrix]


def load_manifest(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("unexpected rho-observer manifest schema")
    return payload


def git_output(arguments: list[str]) -> bytes:
    process = subprocess.run(
        ["git", *arguments],
        cwd=ROOT,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if process.returncode != 0:
        message = process.stderr.decode("utf-8", errors="replace").strip()
        raise RuntimeError(f"git {' '.join(arguments)} failed: {message}")
    return process.stdout


def read_pinned_object(commit: str, contract: Mapping[str, Any]) -> tuple[bytes, str]:
    path = str(contract["path"])
    object_name = f"{commit}:{path}"
    blob = git_output(["rev-parse", object_name]).decode("ascii").strip()
    if blob != contract["git_blob_sha1"]:
        raise ValueError(f"Git blob mismatch for {path}: {blob}")
    payload = git_output(["show", object_name])
    digest = sha256_bytes(payload)
    if digest != contract["sha256"]:
        raise ValueError(f"SHA256 mismatch for {path}: {digest}")
    return payload, blob


def load_inputs(manifest: Mapping[str, Any]) -> dict[str, Any]:
    source = manifest["source"]
    commit = str(source["commit"])
    resolved_commit = git_output(["rev-parse", commit]).decode("ascii").strip()
    if resolved_commit != commit:
        raise ValueError(f"commit does not resolve exactly: {resolved_commit}")

    batch_bytes, batch_blob = read_pinned_object(commit, source["batches"])
    score_bytes, score_blob = read_pinned_object(commit, source["score"])
    score = json.loads(score_bytes.decode("utf-8"))
    if score.get("batch_sha256") != source["batches"]["sha256"]:
        raise ValueError("score envelope does not reference the pinned batches")
    if score.get("schema") != "matching-one/p267-rho-child-etop-c3-score/v1":
        raise ValueError("unexpected frozen score schema")
    return {
        "commit": commit,
        "batch_bytes": batch_bytes,
        "score": score,
        "provenance": {
            "commit": commit,
            "batches": {
                "path": source["batches"]["path"],
                "git_blob_sha1": batch_blob,
                "sha256": source["batches"]["sha256"],
            },
            "score": {
                "path": source["score"]["path"],
                "git_blob_sha1": score_blob,
                "sha256": source["score"]["sha256"],
            },
            "dependency_group": source["dependency_group"],
            "acquisition": source["acquisition"],
        },
    }


def reconstruct_batches(batch_bytes: bytes, manifest: Mapping[str, Any]) -> dict[str, Any]:
    source = manifest["source"]
    children = list(manifest["observable"]["geometry"])
    stream = io.StringIO(batch_bytes.decode("utf-8"), newline="")
    rows = list(csv.DictReader(stream))
    if len(rows) != int(source["expected_batches"]):
        raise ValueError("unexpected batch count")

    vectors: list[list[float]] = []
    total_samples = 0
    replica_intervals: list[tuple[int, int]] = []
    common_field_hashes: list[str] = []
    for position, row in enumerate(rows):
        if int(row["batch"]) != position:
            raise ValueError("batch order is not contiguous")
        samples = int(row["samples"])
        if samples != int(source["expected_samples_per_batch"]):
            raise ValueError("batch sample count differs from the frozen contract")
        replica_first = int(row["replica_first"])
        replica_intervals.append((replica_first, replica_first + samples))
        if position and replica_intervals[-2][1] != replica_first:
            raise ValueError("replica intervals are not contiguous")
        total_samples += samples
        common_field_hashes.append(row["common_field_sha256"])

        vector: list[float] = []
        for child in children:
            invalid = int(row[f"{child}_invalid"])
            if invalid != 0:
                raise ValueError(f"nonzero topology invariant failures for {child}")
            rank_total = sum(int(row[f"{child}_rank{rank}"]) for rank in range(3))
            if rank_total != samples:
                raise ValueError(f"rank counts do not sum to samples for {child}")
            etop = (
                int(row[f"{child}_rank0"]) + int(row[f"{child}_rank2"])
            ) / samples
            h4_re = float(row[f"{child}_H4_re"]) / samples
            h4_im = float(row[f"{child}_H4_im"]) / samples
            vector.extend([etop, h4_re, h4_im])
        vectors.append(vector)

    if total_samples != int(source["expected_total_samples"]):
        raise ValueError("total samples differ from the frozen contract")
    if len(set(common_field_hashes)) != len(common_field_hashes):
        raise ValueError("duplicate common-field batch hash")

    batch_matrix = np.asarray(vectors, dtype=float)
    mean = np.mean(batch_matrix, axis=0)
    centered = batch_matrix - mean
    count = len(batch_matrix)
    covariance = centered.T @ centered / (count * (count - 1))
    covariance = (covariance + covariance.T) / 2.0
    eigenvalues = np.linalg.eigvalsh(covariance)
    if float(np.min(eigenvalues)) <= 0:
        raise ValueError("reconstructed covariance is not positive definite")
    return {
        "batch_matrix": batch_matrix,
        "mean": mean,
        "covariance": covariance,
        "eigenvalues": eigenvalues,
        "condition_number": float(np.linalg.cond(covariance)),
        "batch_count": count,
        "total_samples": total_samples,
        "replica_interval": [replica_intervals[0][0], replica_intervals[-1][1]],
        "unique_common_field_hashes": len(set(common_field_hashes)),
    }


def continuum_vector(score: Mapping[str, Any], manifest: Mapping[str, Any]) -> np.ndarray:
    expected_children = list(manifest["observable"]["geometry"])
    if score.get("child_order") != expected_children:
        raise ValueError("score child order differs from the manifest")
    rows = score.get("continuum")
    if not isinstance(rows, list) or len(rows) != len(expected_children):
        raise ValueError("frozen score continuum is incomplete")
    vector: list[float] = []
    for expected, row in zip(expected_children, rows):
        if row.get("child") != expected:
            raise ValueError("continuum child order mismatch")
        vector.extend(
            [float(row["Etop"]), float(row["H4_re_im"][0]), float(row["H4_re_im"][1])]
        )
    return np.asarray(vector, dtype=float)


def build_designs(reference: np.ndarray) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    designs = {
        "normalizer_only": np.zeros((9, 3)),
        "rank1_mass_only": np.zeros((9, 3)),
        "independent_real_rescalings": np.zeros((9, 6)),
    }
    constraints = {
        "normalizer_only": np.zeros((6, 9)),
        "rank1_mass_only": np.zeros((6, 9)),
        "independent_real_rescalings": np.zeros((3, 9)),
    }
    for child in range(3):
        start = 3 * child
        etop, h4_re, h4_im = reference[start : start + 3]
        h4 = np.asarray([h4_re, h4_im])
        rank1_mass = 1.0 - etop
        if etop <= 0 or rank1_mass <= 0 or float(np.linalg.norm(h4)) == 0:
            raise ValueError("degenerate continuum reference")

        designs["normalizer_only"][start : start + 3, child] = reference[
            start : start + 3
        ]
        designs["rank1_mass_only"][start : start + 3, child] = np.r_[
            1.0, -h4 / rank1_mass
        ]
        designs["independent_real_rescalings"][start, 2 * child] = 1.0
        designs["independent_real_rescalings"][
            start + 1 : start + 3, 2 * child + 1
        ] = h4

        for component in range(2):
            constraints["normalizer_only"][2 * child + component, start] = -h4[
                component
            ]
            constraints["normalizer_only"][
                2 * child + component, start + 1 + component
            ] = etop
            constraints["rank1_mass_only"][2 * child + component, start] = h4[
                component
            ]
            constraints["rank1_mass_only"][
                2 * child + component, start + 1 + component
            ] = rank1_mass
        constraints["independent_real_rescalings"][
            child, start + 1 : start + 3
        ] = [-h4_im, h4_re]
    return designs, constraints


def chi_square_survival(value: float, degrees: int) -> float:
    mp.mp.dps = max(mp.mp.dps, 60)
    shape = mp.mpf(degrees) / 2
    output = mp.gammainc(shape, mp.mpf(str(value)) / 2, mp.inf) / mp.gamma(shape)
    return finite_float(output)


def f_survival(value: float, numerator_df: int, denominator_df: int) -> float:
    mp.mp.dps = max(mp.mp.dps, 60)
    x = mp.mpf(str(value))
    d1 = mp.mpf(numerator_df)
    d2 = mp.mpf(denominator_df)
    z = d2 / (d2 + d1 * x)
    output = mp.betainc(d2 / 2, d1 / 2, 0, z, regularized=True)
    return finite_float(output)


def score_model(
    mean: np.ndarray,
    reference: np.ndarray,
    covariance: np.ndarray,
    design: np.ndarray,
    constraint: np.ndarray,
    batches: int,
    alpha: float,
) -> dict[str, Any]:
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("model design is rank deficient")
    if float(np.max(np.abs(constraint @ design))) > 2e-14:
        raise ValueError("declared constraints do not annihilate the design")

    delta = mean - reference
    cholesky = np.linalg.cholesky(covariance)
    whitened_design = np.linalg.solve(cholesky, design)
    whitened_delta = np.linalg.solve(cholesky, delta)
    coefficients = np.linalg.lstsq(whitened_design, whitened_delta, rcond=None)[0]
    whitened_residual = whitened_delta - whitened_design @ coefficients
    statistic = float(whitened_residual @ whitened_residual)
    degrees = int(len(mean) - design.shape[1])

    constraint_residual = constraint @ delta
    constraint_covariance = constraint @ covariance @ constraint.T
    constraint_eigenvalues = np.linalg.eigvalsh(
        (constraint_covariance + constraint_covariance.T) / 2.0
    )
    if float(np.min(constraint_eigenvalues)) <= 0:
        raise ValueError("constraint covariance is not positive definite")
    independent = float(
        constraint_residual
        @ np.linalg.solve(constraint_covariance, constraint_residual)
    )
    if not np.isclose(statistic, independent, rtol=3e-12, atol=1e-10):
        raise ValueError("GLS and exact-constraint scores disagree")

    denominator_df = batches - degrees
    if denominator_df <= 0:
        raise ValueError("too few batches for Hotelling reference")
    hotelling_f = statistic * denominator_df / (degrees * (batches - 1))
    hotelling_p = f_survival(hotelling_f, degrees, denominator_df)
    nominal_p = chi_square_survival(statistic, degrees)
    return {
        "nuisance_dimension": int(design.shape[1]),
        "constraint_dimension": degrees,
        "design_rank": int(np.linalg.matrix_rank(design)),
        "mahalanobis_T2": finite_float(statistic),
        "nominal_chi_square_survival_p": nominal_p,
        "Hotelling_F": finite_float(hotelling_f),
        "Hotelling_degrees_of_freedom": [degrees, denominator_df],
        "Hotelling_p_gaussian_batch_reference": hotelling_p,
        "decision_alpha": alpha,
        "excluded_at_alpha_under_gaussian_batch_reference": bool(hotelling_p < alpha),
        "coefficients": vector_payload(coefficients),
        "prediction": vector_payload(reference + design @ coefficients),
        "residual": vector_payload(delta - design @ coefficients),
        "linear_constraints": matrix_payload(constraint),
        "constraint_residual": vector_payload(constraint_residual),
        "constraint_covariance": matrix_payload(constraint_covariance),
        "constraint_covariance_eigenvalues": vector_payload(constraint_eigenvalues),
        "independent_constraint_T2": finite_float(independent),
    }


def replay_frozen_score(
    mean: np.ndarray,
    covariance: np.ndarray,
    reference: np.ndarray,
    score: Mapping[str, Any],
) -> dict[str, Any]:
    delta_e = mean[::3] - reference[::3]
    phase = np.exp(-2j * np.pi * np.arange(3) / 3) / 3
    dft = np.vstack([phase.real, phase.imag])
    primary = dft @ delta_e
    e_indices = [0, 3, 6]
    primary_covariance = dft @ covariance[np.ix_(e_indices, e_indices)] @ dft.T
    primary_t2 = float(primary @ np.linalg.solve(primary_covariance, primary))
    archived_primary = score["primary_nontrivial_Etop_r1"]
    if not np.allclose(primary, archived_primary["value_re_im"], rtol=0, atol=3e-16):
        raise ValueError("reconstructed primary value differs from frozen score")
    if not np.allclose(
        primary_covariance, archived_primary["covariance_2x2"], rtol=1e-10, atol=1e-22
    ):
        raise ValueError("reconstructed primary covariance differs from frozen score")
    if abs(primary_t2 - float(archived_primary["chi_square"])) > 1e-8:
        raise ValueError("reconstructed primary T2 differs from frozen score")

    vectors: list[np.ndarray] = []
    for kind, character in (("E", 0), ("E", 1), ("H", 0), ("H", 1)):
        vector = np.zeros(9, dtype=complex)
        for child in range(3):
            weight = 1 / 3 if character == 0 else phase[child]
            if kind == "E":
                vector[3 * child] = weight
            else:
                vector[3 * child + 1] = weight
                vector[3 * child + 2] = 1j * weight
        vectors.append(vector)
    e0, e1, h0, h1 = vectors
    residual = mean - reference
    plugin = (h0 @ residual) * (e1 @ residual) - (h1 @ residual) * (e0 @ residual)
    bias = h0 @ covariance @ e1 - h1 @ covariance @ e0
    corrected = plugin - bias
    archived_determinant = score["secondary_observer_ray_determinant"]["value_re_im"]
    if not np.allclose(
        [corrected.real, corrected.imag], archived_determinant, rtol=0, atol=1e-16
    ):
        raise ValueError("reconstructed determinant point differs from frozen score")
    return {
        "primary_Etop_r1": {
            "value_re_im": vector_payload(primary),
            "covariance": matrix_payload(primary_covariance),
            "T2": finite_float(primary_t2),
            "archived_T2": finite_float(archived_primary["chi_square"]),
        },
        "secondary_observer_ray_point": {
            "plugin": [finite_float(plugin.real), finite_float(plugin.imag)],
            "mean_covariance_bias": [finite_float(bias.real), finite_float(bias.imag)],
            "bias_corrected": [finite_float(corrected.real), finite_float(corrected.imag)],
            "archived_bias_corrected": [finite_float(value) for value in archived_determinant],
            "nonlinear_jackknife_covariance_recomputed": False,
        },
    }


def analyze(manifest_path: Path) -> dict[str, Any]:
    started = time.perf_counter()
    manifest = load_manifest(manifest_path)
    loaded = load_inputs(manifest)
    batches = reconstruct_batches(loaded["batch_bytes"], manifest)
    reference = continuum_vector(loaded["score"], manifest)
    designs, constraints = build_designs(reference)
    alpha = float(manifest["statistics"]["decision_alpha"])
    model_results = {
        name: score_model(
            batches["mean"],
            reference,
            batches["covariance"],
            designs[name],
            constraints[name],
            batches["batch_count"],
            alpha,
        )
        for name in ("normalizer_only", "rank1_mass_only", "independent_real_rescalings")
    }
    replay = replay_frozen_score(
        batches["mean"], batches["covariance"], reference, loaded["score"]
    )
    decision = (
        "ALL_THREE_DECLARED_PRODUCTION_PARAMETERIZATIONS_EXCLUDED"
        if all(
            row["excluded_at_alpha_under_gaussian_batch_reference"]
            for row in model_results.values()
        )
        else "AT_LEAST_ONE_DECLARED_PARAMETERIZATION_NOT_EXCLUDED"
    )
    return {
        "schema": OUTPUT_SCHEMA,
        "status": "post_reveal_existing_production_reanalysis",
        "issue": int(manifest["issue"]),
        "decision": decision,
        "provenance": loaded["provenance"],
        "observable_contract": manifest["observable"],
        "reconstruction": {
            "batch_count": batches["batch_count"],
            "total_samples": batches["total_samples"],
            "replica_interval": batches["replica_interval"],
            "unique_common_field_hashes": batches["unique_common_field_hashes"],
            "coordinate_order": manifest["observable"]["coordinate_order"],
            "mean": vector_payload(batches["mean"]),
            "continuum_reference": vector_payload(reference),
            "covariance_of_mean": matrix_payload(batches["covariance"]),
            "covariance_eigenvalues": vector_payload(batches["eigenvalues"]),
            "covariance_condition_number": finite_float(batches["condition_number"]),
        },
        "frozen_score_replay": replay,
        "models": model_results,
        "statistics_boundary": manifest["statistics"],
        "dependency_interpretation": {
            "evidence_units": 1,
            "models_are_correlated_reuses": True,
            "p_values_must_not_be_combined": True,
        },
        "claim_boundary": manifest["claim_boundary"],
        "runtime": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "wall_seconds": finite_float(time.perf_counter() - started),
            "new_random_samples": 0,
        },
    }


def render_markdown(payload: Mapping[str, Any]) -> str:
    models = payload["models"]
    lines = [
        "# P275 rho-child observer-model elimination",
        "",
        f"Decision: `{payload['decision']}`.",
        "",
        "This is a post-reveal reanalysis of one frozen N112 square-bond dependency block. "
        "It generated no new Monte Carlo samples.",
        "",
        "## Production reconstruction",
        "",
        f"- Commit: `{payload['provenance']['commit']}`",
        f"- Dependency group: `{payload['provenance']['dependency_group']}`",
        f"- Batches / samples: `{payload['reconstruction']['batch_count']}` / "
        f"`{payload['reconstruction']['total_samples']:,}`",
        f"- Covariance condition number: `{payload['reconstruction']['covariance_condition_number']:.8g}`",
        "- The nine-dimensional mean and full covariance of the mean were reconstructed "
        "directly from the pinned `batches.csv`; the archived Etop-r1 and determinant "
        "points were replayed against the pinned `score.json`.",
        "",
        "## Declared model scores",
        "",
        "| parameterization | T2 / constraints | Hotelling F reference p | decision at 0.01 |",
        "|---|---:|---:|---|",
    ]
    for name in ("normalizer_only", "rank1_mass_only", "independent_real_rescalings"):
        row = models[name]
        lines.append(
            f"| `{name}` | {row['mahalanobis_T2']:.6f} / "
            f"{row['constraint_dimension']} | {row['Hotelling_p_gaussian_batch_reference']:.6g} | "
            f"{'excluded' if row['excluded_at_alpha_under_gaussian_batch_reference'] else 'not excluded'} |"
        )
    lines.extend(
        [
            "",
            "The fixed probability identities rejected here are: common-denominator-only "
            "scaling, rank-1-total-mass-only reweighting with fixed internal winding composition, "
            "and the broader class in which E is free but primitive H4 may only undergo a real "
            "rescaling on each child.",
            "",
            "## Interpretation boundary",
            "",
            "All three scores reuse the same 100 aligned batches and form one evidence unit. "
            "The Hotelling values are finite-batch references under iid Gaussian batch means; "
            "the hypotheses were defined after reveal and their p-values are not prospective or additive.",
            "",
            "These exclusions are parameterization-specific. `E_top` is an Alexander-even topology "
            "coordinate, not an identified energy operator; primitive `H4` is a direction-weighted "
            "observer, not a local spin-4 identification. The result does not decide square-site "
            "original U, Q4/Jordan identity, H4 versus H8, or an asymptotic exponent.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--output-md", type=Path, default=DEFAULT_MD)
    arguments = parser.parse_args()

    payload = analyze(arguments.manifest)
    arguments.output_json.parent.mkdir(parents=True, exist_ok=True)
    arguments.output_md.parent.mkdir(parents=True, exist_ok=True)
    arguments.output_json.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    arguments.output_md.write_text(render_markdown(payload), encoding="utf-8")
    print(
        json.dumps(
            {
                "schema": payload["schema"],
                "decision": payload["decision"],
                "models": {
                    name: {
                        "T2": row["mahalanobis_T2"],
                        "df": row["constraint_dimension"],
                        "Hotelling_p": row["Hotelling_p_gaussian_batch_reference"],
                    }
                    for name, row in payload["models"].items()
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
