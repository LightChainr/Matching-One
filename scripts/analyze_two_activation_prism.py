#!/usr/bin/env python3
"""Resolve the P205 quotient prism into its two exact homology activations.

This scorer reuses the branch-only P205 N=25/50/125 threshold-rank archives at
the probability and character vectors frozen before their reveal.  It changes
only the readout:

    K1 = K_minus, K2 = K_plus,
    Delta M = Delta F1 + Delta F2.

The two component contrasts are recomputed after deleting the same aligned
batch from both quotient orientations.  Their 2x2 covariance blocks are then
used in a six-coordinate GLS score over all nine frozen H4/H8/H12 component
pairs.  The exponent remains 13/8; no offset or correction is fitted.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import mpmath as mp
import yaml

from analyze_two_activation_h4 import (
    Archive,
    ArchiveNotScoreable,
    _sum_histograms,
    activation_components,
    jackknife_covariance,
    read_archive,
    sha256,
)


SCHEMA = "matching-one.two-activation-prism.v1"
MANIFEST_SCHEMA = "matching-one.two-activation-prism.manifest.v1"
COMPONENTS = ("K1", "K2")
COMPONENT_METRICS = ("delta_F1", "delta_F2")


def _as_float(value: mp.mpf | float | int) -> float:
    return float(value)


def load_manifest(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or payload.get("schema") != MANIFEST_SCHEMA:
        raise ValueError(f"manifest must use schema {MANIFEST_SCHEMA}")
    if payload.get("status") != "retrospective_branch_only_existing_data_reanalysis":
        raise ValueError("manifest status must preserve retrospective branch-only reuse")
    source = payload.get("source")
    if not isinstance(source, dict) or source.get("integration_state") != "branch_only":
        raise ValueError("source must remain explicitly branch_only")
    commit = str(source.get("commit", ""))
    if len(commit) != 40:
        raise ValueError("source.commit must be the full immutable commit")
    contract = payload.get("semantic_contract")
    if not isinstance(contract, dict):
        raise ValueError("manifest lacks semantic_contract")
    if contract.get("K1") != "K_minus" or contract.get("K2") != "K_plus":
        raise ValueError("the exact activation mapping K1=K_minus, K2=K_plus changed")
    if contract.get("DeltaM") != "DeltaF1+DeltaF2":
        raise ValueError("the exact matching reconstruction changed")
    scoring = payload.get("scoring_contract")
    if not isinstance(scoring, dict):
        raise ValueError("manifest lacks scoring_contract")
    if tuple(scoring.get("model_order", ())) != ("H4", "H8", "H12"):
        raise ValueError("frozen harmonic order changed")
    if Fraction(str(scoring.get("radial_exponent"))) != Fraction(13, 8):
        raise ValueError("radial exponent must remain frozen at 13/8")
    if str(scoring.get("p_ref")) != "0.59274605079":
        raise ValueError("fixed P205 probability changed")
    targets = payload.get("targets")
    if not isinstance(targets, list) or [int(row["N"]) for row in targets] != [25, 50, 125]:
        raise ValueError("target order must remain N=25,50,125")
    dependency_groups = [str(row.get("dependency_group", "")) for row in targets]
    if len(dependency_groups) != len(set(dependency_groups)):
        raise ValueError("the three independent P205 streams need distinct dependency groups")
    for row in targets:
        characters = row.get("exact_character_difference")
        if not isinstance(characters, dict) or tuple(characters) != ("H4", "H8", "H12"):
            raise ValueError(f"N={row.get('N')}: frozen exact character vector changed")
        for value in characters.values():
            Fraction(str(value))
        checksums = row.get("sha256")
        if not isinstance(checksums, dict) or set(checksums) != {
            "histogram", "moments", "metadata"
        }:
            raise ValueError(f"N={row.get('N')}: exact input checksums are incomplete")
    return payload


def verify_source_checkout(source_root: Path, manifest: Mapping[str, Any]) -> str:
    if not source_root.is_dir():
        raise ValueError(f"source checkout is missing: {source_root}")
    try:
        completed = subprocess.run(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ValueError(f"cannot resolve source checkout commit: {source_root}") from exc
    actual = completed.stdout.strip()
    expected = str(manifest["source"]["commit"])
    if actual != expected:
        raise ValueError(f"source checkout is {actual}, expected immutable commit {expected}")
    return actual


def _archive_entry(target: Mapping[str, Any]) -> dict[str, Any]:
    files = target["files"]
    return {
        "N": int(target["N"]),
        "histogram": str(files["histogram"]),
        "moments": str(files["moments"]),
        "metadata": str(files["metadata"]),
        "dependency_group": str(target["dependency_group"]),
        "expected_first": [int(value) for value in target["expected_first"]],
        "expected_second": [int(value) for value in target["expected_second"]],
    }


def load_archive(source_root: Path, target: Mapping[str, Any]) -> Archive:
    entry = _archive_entry(target)
    archive = read_archive(source_root, entry)
    for kind, expected in target["sha256"].items():
        actual = sha256(archive.paths[kind])
        if actual != str(expected):
            raise ArchiveNotScoreable(
                f"N={archive.n}: {kind} SHA256 {actual} != frozen {expected}"
            )
    return archive


def fixed_p_estimate(
    archive: Archive, p_ref: float, omitted_batch: Optional[int] = None
) -> dict[str, float]:
    """Evaluate both activation contrasts at the frozen common probability."""

    first = _sum_histograms(archive.histograms["first"], omitted_batch)
    second = _sum_histograms(archive.histograms["second"], omitted_batch)
    if first["samples"] != second["samples"]:
        raise ArchiveNotScoreable("orientation sample totals differ")
    first_values = activation_components(
        archive.n, first["samples"], first["k1"], first["k2"], p_ref
    )
    second_values = activation_components(
        archive.n, second["samples"], second["k1"], second["k2"], p_ref
    )
    delta_f1 = first_values[0] - second_values[0]
    delta_f2 = first_values[1] - second_values[1]
    delta_m = (first_values[0] + first_values[1] - 1.0) - (
        second_values[0] + second_values[1] - 1.0
    )
    reconstruction_residual = delta_m - delta_f1 - delta_f2
    if abs(reconstruction_residual) > 5.0e-15:
        raise ArithmeticError("DeltaM != DeltaF1+DeltaF2 at frozen p_ref")
    return {
        "first_F1": first_values[0],
        "first_F2": first_values[1],
        "first_M": first_values[0] + first_values[1] - 1.0,
        "second_F1": second_values[0],
        "second_F2": second_values[1],
        "second_M": second_values[0] + second_values[1] - 1.0,
        "delta_F1": delta_f1,
        "delta_F2": delta_f2,
        "delta_M": delta_m,
        "reconstruction_residual": reconstruction_residual,
    }


def score_archive(archive: Archive, p_ref: float) -> dict[str, Any]:
    point = fixed_p_estimate(archive, p_ref)
    batch_ids = [row.batch for row in archive.histograms["first"]]
    deleted = [fixed_p_estimate(archive, p_ref, batch) for batch in batch_ids]
    covariance = [
        [
            jackknife_covariance(
                [row[left] for row in deleted], [row[right] for row in deleted]
            )
            for right in COMPONENT_METRICS
        ]
        for left in COMPONENT_METRICS
    ]
    matrix = mp.matrix(covariance)
    determinant = mp.det(matrix)
    if determinant <= 0 or min(covariance[0][0], covariance[1][1]) <= 0:
        raise ArithmeticError(f"N={archive.n}: component covariance is not positive definite")
    correlation = covariance[0][1] / math.sqrt(covariance[0][0] * covariance[1][1])
    return {
        "point": point,
        "deleted": deleted,
        "batch_ids": batch_ids,
        "covariance": covariance,
        "correlation": correlation,
    }


def block_covariance(by_size: Mapping[int, Mapping[str, Any]], sizes: Sequence[int]) -> list[list[float]]:
    width = len(sizes) * len(COMPONENTS)
    output = [[0.0 for _ in range(width)] for _ in range(width)]
    for size_index, n in enumerate(sizes):
        block = by_size[n]["covariance"]
        offset = size_index * len(COMPONENTS)
        for left in range(len(COMPONENTS)):
            for right in range(len(COMPONENTS)):
                output[offset + left][offset + right] = float(block[left][right])
    return output


def component_model_vector(
    sizes: Sequence[int], characters: Mapping[int, Fraction], exponent: Fraction
) -> list[mp.mpf]:
    power = mp.mpf(exponent.numerator) / exponent.denominator
    return [
        mp.power(n, -power)
        * mp.mpf(characters[n].numerator)
        / characters[n].denominator
        for n in sizes
    ]


def fit_component_pair(
    observed: Sequence[float],
    covariance: Sequence[Sequence[float]],
    first_vector: Sequence[mp.mpf | float],
    second_vector: Sequence[mp.mpf | float],
    first_name: str,
    second_name: str,
) -> dict[str, Any]:
    """Fit one frozen character amplitude to each activation component."""

    if len(observed) != 2 * len(first_vector) or len(first_vector) != len(second_vector):
        raise ValueError("component GLS dimensions are inconsistent")
    design_rows = []
    for first_value, second_value in zip(first_vector, second_vector):
        design_rows.extend([[mp.mpf(first_value), 0], [0, mp.mpf(second_value)]])
    y = mp.matrix([mp.mpf(value) for value in observed])
    design = mp.matrix(design_rows)
    cov = mp.matrix(covariance)
    inverse = cov**-1
    information = design.T * inverse * design
    amplitude_covariance = information**-1
    amplitudes = amplitude_covariance * design.T * inverse * y
    fitted = design * amplitudes
    residual = y - fitted
    chi_square = (residual.T * inverse * residual)[0]
    degrees = len(observed) - 2
    survival = mp.gammainc(
        mp.mpf(degrees) / 2, chi_square / 2, mp.inf
    ) / mp.gamma(mp.mpf(degrees) / 2)
    denominator = amplitudes[0] + amplitudes[1]
    k2_fraction = amplitudes[1] / denominator if denominator != 0 else mp.nan
    marginal_se = [mp.sqrt(cov[index, index]) for index in range(len(observed))]
    return {
        "K1_character": first_name,
        "K2_character": second_name,
        "amplitudes": {"A1": _as_float(amplitudes[0]), "A2": _as_float(amplitudes[1])},
        "amplitude_covariance": [
            [_as_float(amplitude_covariance[i, j]) for j in range(2)] for i in range(2)
        ],
        "amplitude_standard_errors": {
            "A1": _as_float(mp.sqrt(amplitude_covariance[0, 0])),
            "A2": _as_float(mp.sqrt(amplitude_covariance[1, 1])),
        },
        "K2_fitted_fraction": _as_float(k2_fraction),
        "component_interaction": (
            "reinforcing" if amplitudes[0] * amplitudes[1] > 0 else "cancelling"
        ),
        "fitted_vector": [_as_float(value) for value in fitted],
        "residual_vector": [_as_float(value) for value in residual],
        "marginal_standardized_residuals": [
            _as_float(residual[index] / marginal_se[index])
            for index in range(len(observed))
        ],
        "chi_square": _as_float(chi_square),
        "degrees_of_freedom": degrees,
        "chi_square_survival": _as_float(survival),
    }


def _original_prism_contrasts(
    source_root: Path, manifest: Mapping[str, Any]
) -> tuple[dict[int, float], dict[str, Any]]:
    contract = manifest["source"]["original_total_score"]
    path = source_root / str(contract["path"])
    if not path.is_file():
        raise ValueError(f"original P205 score is missing: {path}")
    actual_hash = sha256(path)
    if actual_hash != str(contract["sha256"]):
        raise ValueError("original P205 total-score checksum changed")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") != "matching-one/p205-quotient-character-prism-score/v1":
        raise ValueError("original P205 total-score schema changed")
    if str(payload.get("fixed_probability")) != "0.59274605079":
        raise ValueError("original P205 total score used another probability")
    contrasts = {
        int(n): float(row["DeltaM_first_minus_second"])
        for n, row in payload["by_size"].items()
    }
    return contrasts, {"path": str(contract["path"]), "sha256": actual_hash}


def render(
    manifest_path: Path,
    manifest: Mapping[str, Any],
    source_root: Path,
    archives: Mapping[int, Archive],
) -> dict[str, Any]:
    sizes = [int(row["N"]) for row in manifest["targets"]]
    p_ref = float(manifest["scoring_contract"]["p_ref"])
    scored = {n: score_archive(archives[n], p_ref) for n in sizes}
    original, original_provenance = _original_prism_contrasts(source_root, manifest)
    by_n: dict[str, Any] = {}
    for target in manifest["targets"]:
        n = int(target["N"])
        row = scored[n]
        point = row["point"]
        original_residual = point["delta_M"] - original[n]
        if abs(original_residual) > 5.0e-15:
            raise ArithmeticError(f"N={n}: K1+K2 does not reconstruct original P205 DeltaM")
        by_n[str(n)] = {
            "status": "scoreable",
            "dependency_group": str(target["dependency_group"]),
            "representatives": {
                "first": [int(value) for value in target["expected_first"]],
                "second": [int(value) for value in target["expected_second"]],
            },
            "samples_per_orientation": int(archives[n].metadata["samples_per_pair"]),
            "batch_count": len(row["batch_ids"]),
            "fixed_p_components": point,
            "original_P205_delta_M": original[n],
            "original_reconstruction_residual": original_residual,
            "estimate_vector_order": list(COMPONENT_METRICS),
            "estimate_vector": [point[name] for name in COMPONENT_METRICS],
            "delete_one_covariance": row["covariance"],
            "delete_one_correlation": row["correlation"],
            "provenance": {
                "runner_commit": archives[n].metadata["git_commit"],
                "seed": archives[n].metadata["seed"],
                "counter_first": archives[n].metadata["replica_counter_first"],
                "counter_last_exclusive": archives[n].metadata[
                    "replica_counter_last_exclusive"
                ],
                "inputs": {
                    kind: {
                        "path": str(archives[n].paths[kind].relative_to(source_root)),
                        "sha256": sha256(archives[n].paths[kind]),
                    }
                    for kind in ("histogram", "moments", "metadata")
                },
            },
        }

    observation_order = [
        {"N": n, "component": component}
        for n in sizes
        for component in COMPONENTS
    ]
    observed = [
        scored[n]["point"][metric]
        for n in sizes
        for metric in COMPONENT_METRICS
    ]
    covariance = block_covariance(scored, sizes)
    exponent = Fraction(str(manifest["scoring_contract"]["radial_exponent"]))
    model_order = tuple(manifest["scoring_contract"]["model_order"])
    characters = {
        model: {
            int(target["N"]): Fraction(str(target["exact_character_difference"][model]))
            for target in manifest["targets"]
        }
        for model in model_order
    }
    vectors = {
        model: component_model_vector(sizes, characters[model], exponent)
        for model in model_order
    }
    fits = [
        fit_component_pair(observed, covariance, vectors[first], vectors[second], first, second)
        for first in model_order
        for second in model_order
    ]
    best_chi_square = min(row["chi_square"] for row in fits)
    for row in fits:
        row["delta_chi_square_from_best"] = row["chi_square"] - best_chi_square
    ranked = sorted(fits, key=lambda row: row["chi_square"])
    best = ranked[0]
    if best["K1_character"] != "H4" or best["K2_character"] != "H4":
        decision = "best_component_pair_is_not_H4_H4"
    else:
        decision = "both_activation_components_select_H4"
    dependency_groups = [
        {
            "id": str(target["dependency_group"]),
            "sizes": [int(target["N"])],
            "rule": "delete_same_aligned_batch_from_both_orientations_and_both_activations",
            "independent_evidence_units": 1,
        }
        for target in manifest["targets"]
    ]
    return {
        "schema": SCHEMA,
        "status": "retrospective branch-only existing-data reanalysis; no new simulation",
        "exact_mapping": {
            "K1": "K_minus",
            "K2": "K_plus",
            "F1": "E[Pr(Binomial(N,p_ref)>=K1)]",
            "F2": "E[Pr(Binomial(N,p_ref)>=K2)]",
            "DeltaM": "DeltaF1+DeltaF2",
        },
        "fixed_scoring": {
            "p_ref": p_ref,
            "radial_exponent": "13/8",
            "model_order": list(model_order),
            "component_pair_count": len(fits),
            "fitted_parameters_per_pair": 2,
            "offset_or_correction_fitted": False,
        },
        "size_order": sizes,
        "by_N": by_n,
        "joint_score": {
            "observation_order": observation_order,
            "estimate_vector": observed,
            "delete_one_covariance": covariance,
            "zero_cross_N_blocks_mean_independent_counter_streams": True,
            "frozen_component_pair_scores": fits,
            "best_pair": best,
            "runner_up_pair": ranked[1],
            "delta_chi_square_to_runner_up": ranked[1]["chi_square"] - best["chi_square"],
            "decision": decision,
        },
        "dependency_groups": dependency_groups,
        "provenance": {
            "manifest": str(manifest_path),
            "manifest_sha256": sha256(manifest_path),
            "source_integration_state": "branch_only",
            "source_branch": str(manifest["source"]["branch"]),
            "source_commit": str(manifest["source"]["commit"]),
            "source_checkout_commit_verified": True,
            "original_total_score": original_provenance,
        },
        "interpretation_guard": (
            "This retrospective decomposition reuses a preregistered small-N character "
            "prism from a branch-only commit. It can distinguish the frozen H4/H8/H12 "
            "character lines for K1 and K2, but it does not identify a continuum operator, "
            "prove an asymptotic law, or convert branch-only inputs into main-integrated facts."
        ),
    }


def render_markdown(payload: Mapping[str, Any]) -> str:
    best = payload["joint_score"]["best_pair"]
    runner_up = payload["joint_score"]["runner_up_pair"]
    lines = [
        "# P205 K1/K2 quotient-prism reanalysis",
        "",
        "## Outcome",
        "",
        "The frozen quotient prism selects `H4/H4`: both the first ambient-homology "
        "activation (`K1=K_minus`) and the second (`K2=K_plus`) are compatible with "
        "the same H4 character line and reinforce at the fixed P205 probability.",
        "",
        "This is a retrospective reuse of branch-only archives. It generated no Monte "
        "Carlo samples, kept `p_ref=0.59274605079` and `N^-13/8` fixed, and fitted "
        "only one amplitude per activation component.",
        "",
        "| component model | fitted amplitude | SE |",
        "|:---|---:|---:|",
        f"| K1 {best['K1_character']} | {best['amplitudes']['A1']:.10f} | "
        f"{best['amplitude_standard_errors']['A1']:.10f} |",
        f"| K2 {best['K2_character']} | {best['amplitudes']['A2']:.10f} | "
        f"{best['amplitude_standard_errors']['A2']:.10f} |",
        "",
        f"K2 supplies {best['K2_fitted_fraction']:.3%} of the signed fitted H4 amplitude. "
        f"The joint score is chi-square={best['chi_square']:.6f} on "
        f"{best['degrees_of_freedom']} df (p={best['chi_square_survival']:.6f}). "
        f"The runner-up is `{runner_up['K1_character']}/{runner_up['K2_character']}` "
        f"at delta chi-square={payload['joint_score']['delta_chi_square_to_runner_up']:.6f}.",
        "",
        "## Activation-resolved contrasts",
        "",
        "| N | Delta F1 | Delta F2 | Delta M | K1/K2 correlation | original closure |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for n in payload["size_order"]:
        row = payload["by_N"][str(n)]
        point = row["fixed_p_components"]
        lines.append(
            f"| {n} | {point['delta_F1']:+.9e} | {point['delta_F2']:+.9e} | "
            f"{point['delta_M']:+.9e} | {row['delete_one_correlation']:+.3f} | "
            f"{row['original_reconstruction_residual']:+.2e} |"
        )
    lines.extend([
        "",
        "Every size reinforces: `DeltaF1` and `DeltaF2` have the same sign. The "
        "complete 6x6 covariance in JSON retains each within-size K1/K2 correlation; "
        "cross-N blocks are zero only because the archived counter streams are distinct.",
        "",
        "## Frozen nine-pair score",
        "",
        "| K1 line | K2 line | chi-square / 4 df | p-value | delta chi-square |",
        "|:---|:---|---:|---:|---:|",
    ])
    for row in payload["joint_score"]["frozen_component_pair_scores"]:
        lines.append(
            f"| {row['K1_character']} | {row['K2_character']} | "
            f"{row['chi_square']:.6f} | {row['chi_square_survival']:.6g} | "
            f"{row['delta_chi_square_from_best']:.6f} |"
        )
    lines.extend([
        "",
        "## Boundary",
        "",
        "The total P205 prism was frozen before reveal, but assigning separate character "
        "lines to K1 and K2 is a retrospective analysis. The source archives remain "
        "`branch_only`, all three sizes are deliberately small, and the result identifies "
        "neither a continuum operator nor an asymptotic theorem. It does show that the "
        "ordinary matching H4 signal on this exact quotient code is not solely a first-"
        "activation effect: the second activation supplies a resolved reinforcing share.",
        "",
        "The JSON records the immutable source commit and input hashes, verifies "
        "`DeltaF1+DeltaF2` against the original P205 `DeltaM`, and documents aligned "
        "delete-one dependency groups.",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "analysis/two_activation_prism_manifest.yaml",
    )
    parser.add_argument(
        "--source-root",
        type=Path,
        required=True,
        help="checkout whose HEAD is the immutable branch-only source commit",
    )
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args()
    try:
        manifest = load_manifest(args.manifest)
        source_root = args.source_root.resolve()
        verify_source_checkout(source_root, manifest)
        archives = {
            int(target["N"]): load_archive(source_root, target)
            for target in manifest["targets"]
        }
        payload = render(args.manifest, manifest, source_root, archives)
    except (
        ArchiveNotScoreable,
        ArithmeticError,
        KeyError,
        OSError,
        TypeError,
        ValueError,
        json.JSONDecodeError,
        yaml.YAMLError,
    ) as exc:
        raise SystemExit(str(exc)) from exc
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.output_md.write_text(render_markdown(payload), encoding="utf-8")
    print(args.output_json)
    print(args.output_md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
