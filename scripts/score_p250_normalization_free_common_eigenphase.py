#!/usr/bin/env python3
"""Score a normalization-free common Z5 charged-cubic hand eigenphase."""

from __future__ import annotations

import argparse
import cmath
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from score_z5_charged_multiseparation import (
    denominator_rows,
    jackknife,
    read_batches,
)
from score_z5_charged_threepoint import zero_score
from z5_projective_leg_multiseparation_mc import SCHEMA as RESPONSE_SCHEMA


MANIFEST_SCHEMA = "matching-one/p250-normalization-free-common-eigenphase-retrospective/v1"
SCORE_SCHEMA = "matching-one/p250-normalization-free-common-eigenphase-score/v1"
SEPARATION = 1
CONTRAST_ORDER = ("cross_re", "cross_im", "norm_difference")
OMEGA_ORDER = (
    "Omega113_plus_re", "Omega113_plus_im",
    "Omega113_minus_re", "Omega113_minus_im",
    "Omega122_plus_re", "Omega122_plus_im",
    "Omega122_minus_re", "Omega122_minus_im",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def complex_channels(vector: Sequence[float]) -> tuple[complex, complex, complex, complex]:
    if len(vector) != 8:
        raise ValueError("normalized charged vector must have eight real coordinates")
    return (
        complex(vector[0], vector[1]),
        complex(vector[2], vector[3]),
        complex(vector[4], vector[5]),
        complex(vector[6], vector[7]),
    )


def common_eigenphase_contrast(vector: Sequence[float]) -> list[float]:
    a_plus, a_minus, b_plus, b_minus = complex_channels(vector)
    cross = a_plus * b_minus - a_minus * b_plus
    norm_difference = (
        abs(a_plus) ** 2 + abs(b_plus) ** 2
        - abs(a_minus) ** 2 - abs(b_minus) ** 2
    )
    return [cross.real, cross.imag, norm_difference]


def q_projection(vector: Sequence[float]) -> complex:
    a_plus, a_minus, b_plus, b_minus = complex_channels(vector)
    denominator = abs(a_minus) ** 2 + abs(b_minus) ** 2
    if denominator <= 0.0:
        raise ValueError("minus-hand normalized cubic support is zero")
    return (a_plus * a_minus.conjugate() + b_plus * b_minus.conjugate()) / denominator


def channel_ratios(vector: Sequence[float]) -> tuple[complex, complex]:
    a_plus, a_minus, b_plus, b_minus = complex_channels(vector)
    if not a_minus or not b_minus:
        raise ValueError("channelwise minus-hand support is zero")
    return a_plus / a_minus, b_plus / b_minus


def jackknife_covariance(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    if len(rows) < 2 or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError("jackknife rows must be rectangular and nontrivial")
    center = [sum(row[index] for row in rows) / len(rows) for index in range(len(rows[0]))]
    factor = (len(rows) - 1) / len(rows)
    return [[
        factor * sum(
            (row[i] - center[i]) * (row[j] - center[j]) for row in rows
        )
        for j in range(len(center))
    ] for i in range(len(center))]


def quadratic(covariance: Sequence[Sequence[float]], gradient: Sequence[float]) -> float:
    return sum(
        gradient[i] * covariance[i][j] * gradient[j]
        for i in range(len(gradient)) for j in range(len(gradient))
    )


def complex_summary(point: complex, deleted: Sequence[complex]) -> dict[str, object]:
    covariance = jackknife_covariance([[value.real, value.imag] for value in deleted])
    amplitude = abs(point)
    if amplitude <= 0.0:
        raise ValueError("complex descriptive projection has zero amplitude")
    amplitude_gradient = [point.real / amplitude, point.imag / amplitude]
    phase_gradient = [-point.imag / amplitude**2, point.real / amplitude**2]
    return {
        "complex_re_im": [point.real, point.imag],
        "covariance_re_im": covariance,
        "amplitude": amplitude,
        "amplitude_standard_error": math.sqrt(max(0.0, quadratic(covariance, amplitude_gradient))),
        "phase_radians": cmath.phase(point),
        "phase_standard_error": math.sqrt(max(0.0, quadratic(covariance, phase_gradient))),
    }


def conjugacy_residual(rows: Sequence[Mapping[str, object]]) -> float:
    maximum = 0.0
    for row in rows:
        for hand in ("plus", "minus"):
            for primary, conjugate in (("C113", "C244"), ("C122", "C334")):
                first = complex(
                    float(row[f"d1_{primary}_{hand}_re"]),
                    float(row[f"d1_{primary}_{hand}_im"]),
                )
                second = complex(
                    float(row[f"d1_{conjugate}_{hand}_re"]),
                    float(row[f"d1_{conjugate}_{hand}_im"]),
                )
                maximum = max(maximum, abs(second - first.conjugate()))
    return maximum


def validate_manifest(root: Path, manifest_path: Path) -> dict:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("wrong common-eigenphase manifest schema")
    if manifest.get("status") != "retrospective_protocol_locked_before_common_eigenphase_score":
        raise ValueError("common-eigenphase protocol is not locked")
    for key in ("response", "batches", "prior_support_score", "producer_manifest"):
        path = root / manifest["input"][key]
        if sha256(path) != manifest["input"][key + "_sha256"]:
            raise ValueError(f"frozen digest mismatch for {key}")
    if tuple(manifest["primary_contrast"]["feature_order"]) != CONTRAST_ORDER:
        raise ValueError("primary contrast order changed")
    if int(manifest["eligibility"]["separation"]) != SEPARATION:
        raise ValueError("only the frozen d1 support row is eligible")
    return manifest


def score(root: Path, manifest_path: Path) -> dict[str, object]:
    manifest = validate_manifest(root, manifest_path)
    response = json.loads((root / manifest["input"]["response"]).read_text(encoding="utf-8"))
    if response.get("schema") != RESPONSE_SCHEMA or not response["exact_gate"].get("passed"):
        raise ValueError("wrong or failed projective-leg response")
    observed_run = response["run"]
    expected = manifest["input"]
    run_identity = (
        ("samples", "samples"), ("batches", "batches_count"),
        ("seed", "seed"), ("replica_offset", "counter"),
    )
    for observed_key, expected_key in run_identity[:3]:
        if observed_run[observed_key] != expected[expected_key]:
            raise ValueError(f"run identity mismatch for {observed_key}")
    if [observed_run["replica_offset"], observed_run["replica_last_exclusive"]] != expected["counter"]:
        raise ValueError("counter interval mismatch")

    batches = read_batches(root / manifest["input"]["batches"])
    if len(batches) != expected["batches_count"] or sum(row["samples"] for row in batches) != expected["samples"]:
        raise ValueError("batch coverage mismatch")
    if [row["batch"] for row in batches] != list(range(len(batches))):
        raise ValueError("batch ids are not zero-based contiguous")
    conjugacy = conjugacy_residual(batches)
    if conjugacy > 1e-12:
        raise ValueError("charged conjugacy gate failed")

    denominators = denominator_rows(batches, SEPARATION)
    minimum_pair_z = min(row["abs_z"] for row in denominators.values())
    if minimum_pair_z < float(manifest["eligibility"]["prior_d1_minimum_pair_abs_z"]):
        raise ValueError("recomputed d1 denominator gate weakened")
    local_point, local_covariance, _ = jackknife(batches, SEPARATION, "local_variance")
    support = zero_score(local_point, local_covariance)
    if support["survival_p"] >= float(manifest["eligibility"]["support_alpha"]):
        raise ValueError("d1 cubic support gate no longer passes")

    omega_point, omega_covariance, omega_deleted = jackknife(
        batches, SEPARATION, "separation"
    )
    contrast_point = common_eigenphase_contrast(omega_point)
    contrast_deleted = [common_eigenphase_contrast(row) for row in omega_deleted]
    contrast_covariance = jackknife_covariance(contrast_deleted)
    primary = zero_score(contrast_point, contrast_covariance)
    alpha = float(manifest["primary_contrast"]["decision_alpha"])

    q_point = q_projection(omega_point)
    q_deleted = [q_projection(row) for row in omega_deleted]
    ratios = channel_ratios(omega_point)
    deleted_ratios = [channel_ratios(row) for row in omega_deleted]
    ratio_rows = [
        [pair[0].real, pair[0].imag, pair[1].real, pair[1].imag]
        for pair in deleted_ratios
    ]
    return {
        "schema": SCORE_SCHEMA,
        "status": "retrospective_normalization_free_reveal",
        "manifest": str(manifest_path.relative_to(root)),
        "manifest_sha256": sha256(manifest_path),
        "operator": manifest["eligibility"]["operator"],
        "separation": SEPARATION,
        "input_gates": {
            "samples": expected["samples"],
            "batches": expected["batches_count"],
            "minimum_pair_abs_z": minimum_pair_z,
            "cubic_support_zero_score": support,
            "DFT_conjugacy_max_abs": conjugacy,
        },
        "normalization": manifest["normalization"],
        "omega_order": list(OMEGA_ORDER),
        "omega_point": omega_point,
        "omega_full_shared_batch_covariance_8x8": omega_covariance,
        "primary_common_unit_eigenphase": {
            "null": manifest["primary_contrast"]["null"],
            "contrast_order": list(CONTRAST_ORDER),
            "contrast_point": contrast_point,
            "contrast_full_shared_batch_covariance_3x3": contrast_covariance,
            **primary,
            "decision_alpha": alpha,
            "decision": (
                "common_unit_eigenphase_eliminated"
                if primary["survival_p"] < alpha
                else "common_unit_eigenphase_survives"
            ),
        },
        "descriptive_eigenphase": {
            "common_projection": complex_summary(q_point, q_deleted),
            "channel_order": ["q113_re", "q113_im", "q122_re", "q122_im"],
            "channel_point": [ratios[0].real, ratios[0].imag, ratios[1].real, ratios[1].imag],
            "channel_full_shared_batch_covariance_4x4": jackknife_covariance(ratio_rows),
        },
        "conjugate_channels": {
            "role": "exact controls, not independent observations",
            "relations": ["C244=conj(C113)", "C334=conj(C122)"],
            "maximum_batch_residual": conjugacy,
        },
        "claim_boundary": [
            "This is the charged projective-leg sector, not ordinary global A_top.",
            "Positive two-point normalization cancels field amplitudes while the exact transported Z5 basis retains complex phase.",
            "A surviving relation would not identify a local primary or universal OPE coefficient; rejection removes only this common unit-eigenphase completion.",
            "No H4/H8/H12 model was scored and no new simulation was used.",
        ],
    }


def render(result: Mapping[str, object]) -> str:
    primary = result["primary_common_unit_eigenphase"]
    common = result["descriptive_eigenphase"]["common_projection"]
    lines = [
        "# P250 normalization-free charged common eigenphase", "",
        "The d=1 projective-leg row is used because its pair denominators and eight-real cubic support passed the previously frozen support gate.", "",
        "## Input gates", "",
        f"- minimum two-point denominator: `{result['input_gates']['minimum_pair_abs_z']:.6g} sigma`",
        f"- cubic support: `{result['input_gates']['cubic_support_zero_score']['chi_square']}/{result['input_gates']['cubic_support_zero_score']['degrees_of_freedom']}`, p `{result['input_gates']['cubic_support_zero_score']['survival_p']:.6g}`",
        f"- conjugacy maximum residual: `{result['input_gates']['DFT_conjugacy_max_abs']:.3g}`", "",
        "## Frozen common-unit-eigenphase contrast", "",
        f"- contrast `{primary['contrast_point']}`",
        f"- chi-square `{primary['chi_square']}/{primary['degrees_of_freedom']}`, p `{primary['survival_p']:.6g}`",
        f"- decision: **`{primary['decision']}`**", "",
        "## Descriptive transported phase", "",
        f"- common q projection: `{common['complex_re_im']}`",
        f"- amplitude: `{common['amplitude']:.6g} +/- {common['amplitude_standard_error']:.3g}`",
        f"- phase: `{common['phase_radians']:.6g} +/- {common['phase_standard_error']:.3g}` radians", "",
        "## Boundary", "",
        *[f"- {line}" for line in result["claim_boundary"]], "",
    ]
    return "\n".join(lines)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path,
        default=root / "analysis/p250_normalization_free_common_eigenphase_20260830.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = score(root, args.manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    args.markdown.write_text(render(result), encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
