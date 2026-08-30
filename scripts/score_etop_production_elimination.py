#!/usr/bin/env python3
"""Score frozen E_top model classes on existing production covariance.

The input is the complete cross-size decision vector and aligned-delete-one
covariance emitted by ``analyze_two_activation_h4.py``.  This scorer generates
no Monte Carlo samples.  It applies the exact integer change of basis

    A_top = delta_F1 + delta_F2,   E_top = delta_F2 - delta_F1

and measures the distance from five explicitly bounded model images to the
resulting production estimate.
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
DEFAULT_MANIFEST = ROOT / "analysis" / "etop_production_elimination_manifest.yaml"
OUTPUT_SCHEMA = "matching-one.etop-production-elimination.v1"
TRANSFORM_COORDINATES = ("A_top", "E_top", "F1", "F2")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def _json_number(value: float) -> float:
    if not math.isfinite(float(value)):
        raise ValueError("non-finite value cannot be serialized")
    return float(value)


def _matrix_payload(matrix: np.ndarray) -> list[list[float]]:
    return [[_json_number(value) for value in row] for row in matrix]


def load_manifest(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("manifest must be a mapping")
    if payload.get("schema") != "matching-one.etop-production-elimination.manifest.v1":
        raise ValueError("unexpected E_top manifest schema")
    return payload


def validate_input(payload: Mapping[str, Any], manifest: Mapping[str, Any]) -> None:
    contract = manifest["input"]
    if payload.get("schema") != contract["schema"]:
        raise ValueError("pinned two-activation schema mismatch")
    if payload.get("size_order") != contract["required_sizes"]:
        raise ValueError("pinned production size order mismatch")
    if payload.get("provenance", {}).get("manifest_sha256") != contract["source_manifest_sha256"]:
        raise ValueError("two-activation source manifest hash mismatch")

    decision = payload.get("decision_covariance", {})
    order = decision.get("metric_order_with_N")
    estimates = decision.get("estimate_vector")
    covariance = decision.get("jackknife_covariance")
    if not isinstance(order, list) or not isinstance(estimates, list) or not isinstance(covariance, list):
        raise ValueError("two-activation decision covariance is incomplete")
    if len(order) != len(estimates) or len(covariance) != len(order):
        raise ValueError("two-activation decision dimensions disagree")
    if any(not isinstance(row, list) or len(row) != len(order) for row in covariance):
        raise ValueError("two-activation covariance is not square")

    index: dict[tuple[int, str], int] = {}
    for position, item in enumerate(order):
        key = (int(item["N"]), str(item["metric"]))
        if key in index:
            raise ValueError(f"duplicate decision coordinate {key}")
        index[key] = position
    for n in contract["required_sizes"]:
        for metric in contract["required_metrics"]:
            if (int(n), metric) not in index:
                raise ValueError(f"missing decision coordinate N={n} {metric}")


def build_exact_transform(
    order: Sequence[Mapping[str, Any]], sizes: Sequence[int]
) -> tuple[np.ndarray, list[dict[str, Any]], dict[tuple[int, str], int]]:
    index = {
        (int(item["N"]), str(item["metric"])): position
        for position, item in enumerate(order)
    }
    rows: list[np.ndarray] = []
    output_order: list[dict[str, Any]] = []
    coefficients = {
        "A_top": (1, 1),
        "E_top": (-1, 1),
        "F1": (1, 0),
        "F2": (0, 1),
    }
    for n in sizes:
        i1 = index[(int(n), "angular_delta_F1")]
        i2 = index[(int(n), "angular_delta_F2")]
        for coordinate in TRANSFORM_COORDINATES:
            row = np.zeros(len(order), dtype=float)
            left, right = coefficients[coordinate]
            row[i1] = left
            row[i2] = right
            rows.append(row)
            output_order.append({"N": int(n), "coordinate": coordinate})
    return np.vstack(rows), output_order, index


def identity_audit(
    estimate: np.ndarray,
    covariance: np.ndarray,
    index: Mapping[tuple[int, str], int],
    sizes: Sequence[int],
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    mean_residuals: dict[str, float] = {}
    covariance_residuals: dict[str, float] = {}
    for n in sizes:
        i1 = index[(int(n), "angular_delta_F1")]
        i2 = index[(int(n), "angular_delta_F2")]
        im = index[(int(n), "angular_delta_M")]
        mean_residuals[str(n)] = _json_number(estimate[im] - estimate[i1] - estimate[i2])
        covariance_residuals[str(n)] = _json_number(
            np.max(np.abs(covariance[im, :] - covariance[i1, :] - covariance[i2, :]))
        )
    max_mean = max(abs(value) for value in mean_residuals.values())
    max_covariance = max(abs(value) for value in covariance_residuals.values())
    tolerance = manifest["exact_transform"]["identity_tolerance"]
    if max_mean > float(tolerance["estimate_absolute"]):
        raise ValueError("delta_M estimate identity exceeds frozen tolerance")
    if max_covariance > float(tolerance["covariance_absolute"]):
        raise ValueError("delta_M covariance identity exceeds frozen tolerance")
    return {
        "identity": "angular_delta_M = angular_delta_F1 + angular_delta_F2",
        "mean_residual_by_N": mean_residuals,
        "covariance_row_max_residual_by_N": covariance_residuals,
        "max_abs_mean_residual": max_mean,
        "max_abs_covariance_residual": max_covariance,
        "tolerance": tolerance,
        "passed": True,
    }


def coordinate_rows(
    output_order: Sequence[Mapping[str, Any]], coordinate: str
) -> np.ndarray:
    rows = np.zeros((sum(item["coordinate"] == coordinate for item in output_order), len(output_order)))
    destination = 0
    for source, item in enumerate(output_order):
        if item["coordinate"] == coordinate:
            rows[destination, source] = 1.0
            destination += 1
    return rows


def covariance_inverse(
    covariance: np.ndarray, relative_cutoff: float, absolute_floor: float
) -> tuple[np.ndarray, dict[str, Any]]:
    symmetric = (covariance + covariance.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    scale = float(max(np.max(np.abs(eigenvalues)), 0.0))
    cutoff = max(scale * relative_cutoff, absolute_floor)
    if float(np.min(eigenvalues)) < -cutoff:
        raise ValueError("covariance has a negative eigenvalue beyond the frozen cutoff")
    kept = eigenvalues > cutoff
    rank = int(np.count_nonzero(kept))
    if rank == 0:
        raise ValueError("covariance has zero numerical rank")
    inverse = (eigenvectors[:, kept] / eigenvalues[kept]) @ eigenvectors[:, kept].T
    return inverse, {
        "eigenvalues_ascending": [_json_number(value) for value in eigenvalues],
        "relative_cutoff": relative_cutoff,
        "absolute_floor": absolute_floor,
        "applied_cutoff": cutoff,
        "kept": [bool(value) for value in kept],
        "rank": rank,
        "discarded_modes": int(len(eigenvalues) - rank),
    }


def chi_square_survival(value: float, degrees: int) -> float:
    if value < 0 or degrees <= 0:
        raise ValueError("invalid chi-square arguments")
    mp.mp.dps = max(mp.mp.dps, 50)
    shape = mp.mpf(degrees) / 2
    return float(mp.gammainc(shape, mp.mpf(str(value)) / 2, mp.inf) / mp.gamma(shape))


def chi_square_critical(alpha: float, degrees: int) -> float:
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie in (0,1)")
    lower, upper = 0.0, float(max(2 * degrees, 1))
    while chi_square_survival(upper, degrees) > alpha:
        upper *= 2.0
    for _ in range(100):
        midpoint = (lower + upper) / 2.0
        if chi_square_survival(midpoint, degrees) > alpha:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def score_residual(
    residual: np.ndarray,
    covariance: np.ndarray,
    nuisance_parameters: int,
    numerics: Mapping[str, Any],
) -> dict[str, Any]:
    inverse, spectrum = covariance_inverse(
        covariance,
        float(numerics["eigen_relative_cutoff"]),
        float(numerics["eigen_absolute_floor"]),
    )
    dual = inverse @ residual
    chi_square = float(residual @ dual)
    degrees = spectrum["rank"] - nuisance_parameters
    if degrees <= 0:
        raise ValueError("model consumes every covariance mode")
    decision_alpha = float(numerics["decision_alpha"])
    strong_alpha = float(numerics["strong_reference_alpha"])
    critical = chi_square_critical(decision_alpha, degrees)
    strong_critical = chi_square_critical(strong_alpha, degrees)
    return {
        "residual_vector": [_json_number(value) for value in residual],
        "residual_covariance": _matrix_payload(covariance),
        "spectrum": spectrum,
        "dual_witness_Vplus_r": [_json_number(value) for value in dual],
        "mahalanobis_chi_square": chi_square,
        "degrees_of_freedom": degrees,
        "chi_square_survival_p": chi_square_survival(chi_square, degrees),
        "decision": {
            "alpha": decision_alpha,
            "critical_chi_square": critical,
            "separation_margin": chi_square - critical,
            "model_image_excluded": bool(chi_square > critical),
        },
        "strong_reference": {
            "alpha": strong_alpha,
            "critical_chi_square": strong_critical,
            "separation_margin": chi_square - strong_critical,
            "model_image_excluded": bool(chi_square > strong_critical),
        },
    }


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
    function: Callable[[float], float], lower: float, upper: float, grid_points: int, tolerance: float
) -> tuple[float, float, dict[str, Any]]:
    grid = np.linspace(lower, upper, grid_points)
    values = np.asarray([function(float(point)) for point in grid])
    candidates: list[tuple[float, float]] = [
        (float(grid[0]), float(values[0])),
        (float(grid[-1]), float(values[-1])),
    ]
    local_indices = [
        index
        for index in range(1, len(grid) - 1)
        if values[index] <= values[index - 1] and values[index] <= values[index + 1]
    ]
    for index in local_indices:
        candidates.append(
            golden_minimum(function, float(grid[index - 1]), float(grid[index + 1]), tolerance)
        )
    point, value = min(candidates, key=lambda item: item[1])
    return point, value, {
        "support": [lower, upper],
        "grid_points": grid_points,
        "grid_minimum": {
            "parameter": float(grid[int(np.argmin(values))]),
            "chi_square": float(np.min(values)),
        },
        "local_minima_refined": len(local_indices),
        "refinement_tolerance": tolerance,
        "best_is_boundary": bool(abs(point - lower) <= tolerance or abs(point - upper) <= tolerance),
    }


def _model_contracts(manifest: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    contracts = {str(item["id"]): item for item in manifest["models"]}
    required = {
        "M0_PURE_ALEXANDER_ODD",
        "M1_SECOND_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO",
        "M2_FIRST_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO",
        "M3_COMMON_PROJECTIVE_RANK_PLANE_LINE",
        "M4_SINGLE_FIXED_H4_POWER",
    }
    if set(contracts) != required:
        raise ValueError("frozen E_top model family drift")
    return contracts


def build_report(
    payload: Mapping[str, Any], manifest: Mapping[str, Any], *, input_path: Path, manifest_path: Path
) -> dict[str, Any]:
    validate_input(payload, manifest)
    decision = payload["decision_covariance"]
    order = decision["metric_order_with_N"]
    sizes = [int(value) for value in payload["size_order"]]
    estimate = np.asarray(decision["estimate_vector"], dtype=float)
    covariance = np.asarray(decision["jackknife_covariance"], dtype=float)
    transform, output_order, index = build_exact_transform(order, sizes)
    transformed_estimate = transform @ estimate
    transformed_covariance = transform @ covariance @ transform.T
    audit = identity_audit(estimate, covariance, index, sizes, manifest)
    _transformed_inverse, transformed_spectrum = covariance_inverse(
        transformed_covariance,
        float(manifest["numerics"]["eigen_relative_cutoff"]),
        float(manifest["numerics"]["eigen_absolute_floor"]),
    )
    selectors = {
        coordinate: coordinate_rows(output_order, coordinate)
        for coordinate in TRANSFORM_COORDINATES
    }
    vectors = {name: selector @ transformed_estimate for name, selector in selectors.items()}
    numerics = manifest["numerics"]
    contracts = _model_contracts(manifest)

    models: dict[str, Any] = {}
    fixed = (
        ("M0_PURE_ALEXANDER_ODD", "E_top"),
        ("M1_SECOND_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO", "F2"),
        ("M2_FIRST_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO", "F1"),
    )
    for model_id, coordinate in fixed:
        selector = selectors[coordinate]
        score = score_residual(
            vectors[coordinate], selector @ transformed_covariance @ selector.T, 0, numerics
        )
        models[model_id] = {
            "contract": contracts[model_id],
            "fitted_parameters": {},
            "score": score,
        }

    e_selector, a_selector = selectors["E_top"], selectors["A_top"]
    e_vector, a_vector = vectors["E_top"], vectors["A_top"]
    model3 = contracts["M3_COMMON_PROJECTIVE_RANK_PLANE_LINE"]

    def projective_score(parameter: float) -> float:
        residual_transform = e_selector - parameter * a_selector
        residual = residual_transform @ transformed_estimate
        residual_covariance = residual_transform @ transformed_covariance @ residual_transform.T
        inverse, _spectrum = covariance_inverse(
            residual_covariance,
            float(numerics["eigen_relative_cutoff"]),
            float(numerics["eigen_absolute_floor"]),
        )
        return float(residual @ inverse @ residual)

    lambda_best, _minimum, search = bounded_profile(
        projective_score,
        float(model3["support"][0]),
        float(model3["support"][1]),
        int(model3["grid_points"]),
        float(model3["refinement_tolerance"]),
    )
    projective_transform = e_selector - lambda_best * a_selector
    projective = score_residual(
        projective_transform @ transformed_estimate,
        projective_transform @ transformed_covariance @ projective_transform.T,
        1,
        numerics,
    )
    models[model3["id"]] = {
        "contract": model3,
        "fitted_parameters": {"lambda": lambda_best},
        "search_audit": search,
        "score": projective,
    }

    model4 = contracts["M4_SINGLE_FIXED_H4_POWER"]
    exponent = float(model4["exponent"]["numerator"]) / float(model4["exponent"]["denominator"])
    shape = np.asarray([float(n) ** (-exponent) for n in sizes])
    e_covariance = e_selector @ transformed_covariance @ e_selector.T
    e_inverse, amplitude_spectrum = covariance_inverse(
        e_covariance,
        float(numerics["eigen_relative_cutoff"]),
        float(numerics["eigen_absolute_floor"]),
    )
    denominator = float(shape @ e_inverse @ shape)
    if denominator <= 0:
        raise ValueError("fixed-power amplitude is not identifiable")
    unconstrained_amplitude = float(shape @ e_inverse @ e_vector / denominator)
    lower, upper = (float(value) for value in model4["support"])
    amplitude = min(max(unconstrained_amplitude, lower), upper)
    fixed_power = score_residual(e_vector - amplitude * shape, e_covariance, 1, numerics)
    models[model4["id"]] = {
        "contract": model4,
        "fitted_parameters": {"c": amplitude},
        "fit_audit": {
            "unconstrained_amplitude": unconstrained_amplitude,
            "support": [lower, upper],
            "best_is_boundary": bool(amplitude != unconstrained_amplitude),
            "shape_vector_N_minus_13_over_8": [_json_number(value) for value in shape],
            "amplitude_fit_spectrum": amplitude_spectrum,
        },
        "score": fixed_power,
    }

    by_n = {}
    group_by_size = {
        int(size): group["id"]
        for group in payload["dependency_groups"]
        for size in group["sizes"]
    }
    for position, n in enumerate(sizes):
        base = 4 * position
        by_n[str(n)] = {
            "N": n,
            "dependency_group": group_by_size[n],
            "A_top": _json_number(transformed_estimate[base]),
            "E_top": _json_number(transformed_estimate[base + 1]),
            "F1": _json_number(transformed_estimate[base + 2]),
            "F2": _json_number(transformed_estimate[base + 3]),
            "E_top_standard_error": _json_number(
                math.sqrt(max(0.0, transformed_covariance[base + 1, base + 1]))
            ),
            "source_provenance": payload["by_N"][str(n)]["provenance"],
        }

    return {
        "schema": OUTPUT_SCHEMA,
        "status": manifest["status"],
        "input": {
            "path": str(input_path.relative_to(ROOT)) if input_path.is_relative_to(ROOT) else str(input_path),
            "schema": payload["schema"],
            "sha256": sha256_file(input_path),
            "pinned_artifact_commit": manifest["input"]["artifact_commit"],
            "pinned_sha256": manifest["input"]["sha256"],
            "pinned_source_manifest_sha256": manifest["input"]["source_manifest_sha256"],
        },
        "scorer": {
            "path": "scripts/score_etop_production_elimination.py",
            "sha256": sha256_file(Path(__file__)),
            "manifest_path": str(manifest_path.relative_to(ROOT)) if manifest_path.is_relative_to(ROOT) else str(manifest_path),
            "manifest_sha256": sha256_file(manifest_path),
            "numpy_version": np.__version__,
            "mpmath_version": mp.__version__,
        },
        "exact_transform": {
            "source_coordinate_order": order,
            "output_coordinate_order": output_order,
            "integer_matrix": [[int(value) for value in row] for row in transform],
            "definitions": manifest["exact_transform"]["definitions"],
            "semantic_boundary": manifest["exact_transform"]["semantic_boundary"],
            "identity_audit": audit,
        },
        "transformed_production_block": {
            "estimate_vector": [_json_number(value) for value in transformed_estimate],
            "covariance": _matrix_payload(transformed_covariance),
            "dimension": len(transformed_estimate),
            "spectrum": transformed_spectrum,
        },
        "by_N": by_n,
        "dependency_groups": payload["dependency_groups"],
        "models": models,
        "primary_conclusion": {
            "pure_Alexander_odd_excluded": models["M0_PURE_ALEXANDER_ODD"]["score"]["decision"]["model_image_excluded"],
            "both_activation_directional_H4_responses_nonzero": bool(
                models["M1_SECOND_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO"]["score"]["decision"]["model_image_excluded"]
                and models["M2_FIRST_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO"]["score"]["decision"]["model_image_excluded"]
            ),
            "common_projective_line_excluded_only_on_declared_domain": models[
                "M3_COMMON_PROJECTIVE_RANK_PLANE_LINE"
            ]["score"]["decision"]["model_image_excluded"],
            "fixed_13_over_8_one_amplitude_excluded_only_without_corrections": models[
                "M4_SINGLE_FIXED_H4_POWER"
            ]["score"]["decision"]["model_image_excluded"],
        },
        "certificate_boundary": manifest["claim_boundary"],
    }


def render_markdown(report: Mapping[str, Any]) -> str:
    models = report["models"]
    order = [
        "M0_PURE_ALEXANDER_ODD",
        "M1_SECOND_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO",
        "M2_FIRST_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO",
        "M3_COMMON_PROJECTIVE_RANK_PLANE_LINE",
        "M4_SINGLE_FIXED_H4_POWER",
    ]
    lines = [
        "# Production E_top model elimination",
        "",
        "The existing ten-size threshold-rank production block rejects a pure",
        "Alexander-odd state response. Both activation-resolved directional H4",
        "response components are nonzero on the declared production block,",
        "and neither one common `E=lambda A` line nor one uncorrected",
        "`E=c N^(-13/8)` amplitude describes the complete declared archive set.",
        "",
        "No Monte Carlo samples are generated here. The exact state transform is",
        "",
        "```text",
        "A_top = delta_F1 + delta_F2",
        "E_top = delta_F2 - delta_F1",
        "```",
        "",
        "and every score uses the pinned full cross-size aligned-delete-one covariance.",
        "",
        "| model | fitted parameter | chi-square / df | survival p | alpha=.05 |",
        "|:--|:--|--:|--:|:--|",
    ]
    for model_id in order:
        row = models[model_id]
        parameters = row["fitted_parameters"]
        parameter_text = ", ".join(f"{key}={value:.9g}" for key, value in parameters.items()) or "none"
        score = row["score"]
        decision = "excluded" if score["decision"]["model_image_excluded"] else "not excluded"
        lines.append(
            f"| `{model_id}` | {parameter_text} | {score['mahalanobis_chi_square']:.6f} / "
            f"{score['degrees_of_freedom']} | {score['chi_square_survival_p']:.6g} | {decision} |"
        )
    lines.extend(
        [
            "",
            "## Direction-normalized E_top production coordinates",
            "",
            "| N | E_top | SE | dependency group |",
            "|--:|--:|--:|:--|",
        ]
    )
    for n in sorted(int(value) for value in report["by_N"]):
        row = report["by_N"][str(n)]
        lines.append(
            f"| {n} | {row['E_top']:+.9e} | {row['E_top_standard_error']:.3e} | "
            f"`{row['dependency_group']}` |"
        )
    audit = report["exact_transform"]["identity_audit"]
    m0 = models["M0_PURE_ALEXANDER_ODD"]["score"]
    lines.extend(
        [
            "",
            "## Certificate form and boundary",
            "",
            f"The pure-odd distance is `{m0['mahalanobis_chi_square']:.6f}`. At the",
            f"stored strong reference alpha `{m0['strong_reference']['alpha']:.1e}`, the",
            f"chi-square critical value is `{m0['strong_reference']['critical_chi_square']:.6f}`",
            f"and the separation margin is `{m0['strong_reference']['separation_margin']:.6f}`.",
            "",
            f"The redundant matching-coordinate audit has maximum mean residual "
            f"`{audit['max_abs_mean_residual']:.3e}` and covariance-row residual "
            f"`{audit['max_abs_covariance_residual']:.3e}`.",
            "",
            "This artifact is a hash-bound exact linear transform followed by a",
            "floating jackknife Mahalanobis confidence-set separation. It is not an",
            "exact probability bound, an interval-reconstructed LDL certificate, or",
            "an SOS certificate. M3 excludes only one common projective line across",
            "the declared sizes/geometries; M4 excludes only one uncorrected fixed-power",
            "amplitude. Neither result excludes every H4 or multi-field mechanism.",
            "The F1/F2 rows are activation-resolved directional responses; their",
            "nonzero values do not assert the mere existence of K1 or K2, which was",
            "already part of the input construction.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--output-md", type=Path)
    args = parser.parse_args()

    manifest_path = args.manifest.resolve()
    manifest = load_manifest(manifest_path)
    input_path = (args.input or (ROOT / manifest["input"]["path"])).resolve()
    if sha256_file(input_path) != manifest["input"]["sha256"]:
        raise SystemExit("pinned two-activation input SHA256 mismatch")
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    report = build_report(payload, manifest, input_path=input_path, manifest_path=manifest_path)

    output_json = (args.output_json or (ROOT / manifest["outputs"]["json"])).resolve()
    output_md = (args.output_md or (ROOT / manifest["outputs"]["markdown"])).resolve()
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    output_md.write_text(render_markdown(report), encoding="utf-8")


if __name__ == "__main__":
    main()
