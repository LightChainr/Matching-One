#!/usr/bin/env python3
"""Audit and freeze the semantics-matched Gaussian x annulus rectangle.

The revealed PR277 and P253 numbers are used only for a post-reveal design
audit.  Incompatible rows are never rescaled or renamed into numerical cells.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Sequence

import numpy as np


LAMBDAS = (0.0, 0.5, 1.0)
RADII = (2, 4, 7, 8)
CHANNELS = ("A_plus", "A_minus")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def radial_basis(lam: float, coordinate: float) -> np.ndarray:
    if lam == 0:
        third = 1.0 if abs(coordinate) < 1e-14 else 0.0
    elif lam == 1:
        # The confluent lambda->1 limit is the rank-3/quadratic adversary.
        third = coordinate * coordinate
    else:
        third = lam ** coordinate
    return np.asarray((1.0, coordinate, third), dtype=float)


def radial_holdout_row(lam: float) -> tuple[np.ndarray, np.ndarray]:
    coordinate = np.log2(np.asarray(RADII, dtype=float) / RADII[0])
    calibration = np.stack([radial_basis(lam, value) for value in coordinate[:3]])
    interpolation = radial_basis(lam, coordinate[3]) @ np.linalg.inv(calibration)
    return np.r_[-interpolation, 1.0], interpolation


def gaussian_separation(score: dict) -> dict:
    entries = {
        0.0: score["fixed_diagnostics_after_primary"]["lambda_0_pure_Jordan"],
        0.5: score["primary_lambda_half"],
        1.0: score["fixed_diagnostics_after_primary"]["lambda_1_persistent_curvature"],
    }
    output = {}
    for view in ("scalar_U", "thermal_jet"):
        rows = {lam: entries[lam][view] for lam in LAMBDAS}
        residual = {lam: np.asarray(rows[lam]["residual"], dtype=float) for lam in LAMBDAS}
        covariance = {lam: np.asarray(rows[lam]["covariance"], dtype=float) for lam in LAMBDAS}
        if not np.allclose(residual[0.5], (residual[0.0] + residual[1.0]) / 2, atol=1e-12):
            raise ValueError(f"{view}: lambda residuals are not affine")
        curvature = residual[1.0] - residual[0.0]
        curvature_covariance = 2 * (
            covariance[1.0] + covariance[0.0] - 2 * covariance[0.5]
        )
        labels = rows[0.0]["residual_order"]
        metrics = []
        for index, label in enumerate(labels):
            variance = curvature_covariance[index, index]
            distance_0_to_1 = curvature[index] ** 2 / variance if variance > 0 else 0.0
            metrics.append({
                "row": label,
                "curvature_lambda_1_minus_0": float(curvature[index]),
                "curvature_standard_error": float(math.sqrt(max(variance, 0))),
                "mahalanobis_squared_lambda_0_vs_1": float(distance_0_to_1),
                "minimum_adjacent_0_half_or_half_1": float(distance_0_to_1 / 4),
            })
        precision = np.linalg.pinv(curvature_covariance, rcond=1e-10, hermitian=True)
        output[view] = {
            "rows": metrics,
            "best_individual_row": max(
                metrics, key=lambda row: row["minimum_adjacent_0_half_or_half_1"]
            ),
            "full_view_mahalanobis_squared_lambda_0_vs_1": float(
                curvature @ precision @ curvature
            ),
            "semantic_warning": (
                "This is a design-sensitivity diagnostic inside PR277 only. "
                "It does not make the row synonymous with a P253 landing-H4 row."
            ),
        }
    return output


def annulus_selection(analysis: dict) -> dict:
    block = analysis["contrast_vector"]
    labels = block["order"]
    point = np.asarray(block["point"], dtype=float)
    covariance = np.asarray(block["covariance"], dtype=float)
    result = {}
    for size in (325, 425):
        indices = {
            channel: [labels.index(f"N{size}_R{radius}_Delta_{channel}") for radius in RADII]
            for channel in CHANNELS
        }
        by_channel = {}
        for channel in CHANNELS:
            local = indices[channel]
            values = point[local]
            local_covariance = covariance[np.ix_(local, local)]
            predictions = {}
            for lam in LAMBDAS:
                row, interpolation = radial_holdout_row(lam)
                predictions[str(lam)] = {
                    "R8_prediction_from_R2_R4_R7": float(interpolation @ values[:3]),
                    "R8_residual": float(row @ values),
                    "R8_residual_standard_error": float(
                        math.sqrt(max(row @ local_covariance @ row, 0))
                    ),
                }
            pairs = []
            for first, second in ((0.0, 0.5), (0.5, 1.0), (0.0, 1.0)):
                row_first, h_first = radial_holdout_row(first)
                row_second, h_second = radial_holdout_row(second)
                delta = float((h_first - h_second) @ values[:3])
                pooled = 0.5 * (
                    row_first @ local_covariance @ row_first
                    + row_second @ local_covariance @ row_second
                )
                pairs.append({
                    "lambda_pair": [first, second],
                    "prediction_difference": delta,
                    "pooled_residual_variance": float(pooled),
                    "mahalanobis_squared": float(delta * delta / pooled),
                })
            by_channel[channel] = {
                "predictions": predictions,
                "pairwise_separation": pairs,
                "minimum_adjacent_mahalanobis_squared": min(
                    row["mahalanobis_squared"] for row in pairs
                    if row["lambda_pair"] in ([0.0, 0.5], [0.5, 1.0])
                ),
            }

        joint_indices = indices["A_plus"] + indices["A_minus"]
        joint_point = point[joint_indices]
        joint_covariance = covariance[np.ix_(joint_indices, joint_indices)]
        joint_pairs = []
        for first, second in ((0.0, 0.5), (0.5, 1.0), (0.0, 1.0)):
            row_first, h_first = radial_holdout_row(first)
            row_second, h_second = radial_holdout_row(second)
            residual_first = np.zeros((2, 8))
            residual_second = np.zeros((2, 8))
            prediction_difference = np.zeros((2, 8))
            residual_first[0, :4] = row_first
            residual_first[1, 4:] = row_first
            residual_second[0, :4] = row_second
            residual_second[1, 4:] = row_second
            prediction_difference[0, :3] = h_first - h_second
            prediction_difference[1, 4:7] = h_first - h_second
            delta = prediction_difference @ joint_point
            pooled = 0.5 * (
                residual_first @ joint_covariance @ residual_first.T
                + residual_second @ joint_covariance @ residual_second.T
            )
            joint_pairs.append({
                "lambda_pair": [first, second],
                "prediction_difference_A_plus_A_minus": delta.tolist(),
                "mahalanobis_squared": float(
                    delta @ np.linalg.pinv(pooled, rcond=1e-12, hermitian=True) @ delta
                ),
            })
        result[f"N{size}"] = {
            "channels": by_channel,
            "joint_pairwise_separation": joint_pairs,
            "joint_minimum_adjacent_mahalanobis_squared": min(
                row["mahalanobis_squared"] for row in joint_pairs
                if row["lambda_pair"] in ([0.0, 0.5], [0.5, 1.0])
            ),
        }
    selected = max(
        result, key=lambda key: result[key]["joint_minimum_adjacent_mahalanobis_squared"]
    )
    return {
        "method": (
            "Post-reveal design metric: calibrate the confluent basis on R2/R4/R7, "
            "predict R8, and use the complete A_plus/A_minus covariance."
        ),
        "lambda_basis": {
            "0": "1,n,indicator(n=0)",
            "1/2": "1,n,(1/2)^n",
            "1": "1,n,n^2 confluent persistent-curvature limit",
        },
        "candidates": result,
        "selected_existing_annulus_context": selected,
        "selected_rows": ["A_plus ordinary/matching-even", "A_minus matching-odd"],
        "warning": (
            "The selected N425 block has weak absolute separation; it selects semantics "
            "and geometry for a future rectangle, not a mechanism from existing data."
        ),
    }


def semantic_crosswalk() -> dict:
    return {
        "PR277_scalar_U": {
            "state": "global rank-histogram thermal response",
            "source": "Bernoulli-p derivative reconstructed from a complete threshold curve",
            "observer": "N^(13/8) P4[S_prime]/Mbar_prime, global orientation/rank projector",
            "normalization": "global thermal slope and canonical width",
            "transfer_coordinate": "Gaussian cover generation k=log2(N/N0)",
            "batch_contract": "100 threshold-histogram batches; source lineages partly aligned",
            "matching_sector": "matching-odd thermal scalar side view; not an ordinary A_plus row",
        },
        "PR277_thermal_jet": {
            "state": "width-normalized Hermite-Krawtchouk jet orders 2 through 6",
            "source": "higher Bernoulli-p derivatives of the global threshold curve",
            "observer": "alternating P4[D] even orders and P4[S] odd orders",
            "normalization": "canonical dimensionless width to each derivative order",
            "transfer_coordinate": "Gaussian cover generation k=log2(N/N0)",
            "batch_contract": "100 threshold-histogram/moment batches with complete jet covariance",
            "matching_sector": "matching-odd thermal jet",
        },
        "P253_A_plus": {
            "state": "local landing-shell pivotal H4 amplitude",
            "source": "single-root occupied/unoccupied toggle at fixed p",
            "observer": "(primal_h4+matching_h4)/(primal_pivotal+matching_pivotal), then orientation contrast",
            "normalization": "conditional pivotal-event mass",
            "transfer_coordinate": "annulus n=log2(R/2) at fixed quotient",
            "batch_contract": "200 common-field batches across radii/designs",
            "matching_sector": "ordinary/matching-even",
        },
        "P253_A_minus": {
            "state": "local landing-shell pivotal H4 amplitude",
            "source": "single-root occupied/unoccupied toggle at fixed p",
            "observer": "(primal_h4-matching_h4)/(primal_pivotal+matching_pivotal), then orientation contrast",
            "normalization": "conditional pivotal-event mass",
            "transfer_coordinate": "annulus n=log2(R/2) at fixed quotient",
            "batch_contract": "200 common-field batches across radii/designs",
            "matching_sector": "matching-odd",
        },
        "exact_pairwise_verdict": {
            "same_source_rows": [],
            "same_observer_rows": [],
            "same_normalization_rows": [],
            "numerically_eligible_cross_context_pairs": [],
            "reason": (
                "A global product-measure derivative is not a fixed-p root-toggle conditional "
                "amplitude. Multiplication by N, a pivotal mass, or a power of R cannot change "
                "the sigma-algebra or sufficient statistic."
            ),
        },
    }


def build_report(gaussian_path: Path, annulus_path: Path) -> dict:
    gaussian = json.loads(gaussian_path.read_text(encoding="utf-8"))
    annulus = json.loads(annulus_path.read_text(encoding="utf-8"))
    if gaussian.get("schema") != "matching-one/norm4-generation4-pilot-score/v1":
        raise ValueError("unexpected PR277 score schema")
    if annulus.get("schema") != "matching-one/norm5-multiradius-pivotal-score/v1":
        raise ValueError("unexpected P253 source score schema")
    return {
        "schema": "matching-one/p267-gaussian-annulus-semantic-crosswalk/v1",
        "status": "post_reveal_design_audit_no_numeric_rectangle",
        "issues": [154, 253, 255],
        "source_commits": {"Gaussian_PR277": "3e855ce", "annulus_P253": "3123b73"},
        "source_files": {
            "Gaussian_PR277_score_sha256": sha256(gaussian_path),
            "annulus_old_score_sha256": sha256(annulus_path),
        },
        "semantic_crosswalk": semantic_crosswalk(),
        "Gaussian_PR277_lambda_sensitivity": gaussian_separation(gaussian),
        "annulus_row_selection": annulus_selection(annulus),
        "rectangle": {
            "row_order": ["A_plus ordinary/matching-even", "A_minus matching-odd"],
            "context_order": ["annulus_radius_doubling", "Gaussian_cover_doubling"],
            "cells": {
                "annulus_radius_doubling:A_plus": "existing P253 N425 R2/R4/R7/R8",
                "annulus_radius_doubling:A_minus": "existing P253 N425 R2/R4/R7/R8",
                "Gaussian_cover_doubling:A_plus": "missing",
                "Gaussian_cover_doubling:A_minus": "missing",
            },
            "numerical_score_now": "forbidden_no_common_source_readout_rectangle",
        },
        "frozen_future_score": {
            "Gaussian_residual": "x3-2*x2+x1-lambda*(x2-2*x1+x0) at fixed R=2",
            "annulus_residual": "R8 minus the R2/R4/R7 confluent-basis prediction",
            "lambda_candidates_in_frozen_order": [0, 0.5, 1],
            "shared_generator": "one common lambda for both contexts and both rows",
            "context_enriched": "one lambda per context, shared by A_plus/A_minus inside that context",
            "comparison": (
                "Report all three diagonal and all nine context-pair GLS scores. "
                "Use Delta=min(diagonal)-min(all) only with a frozen Gaussian parametric "
                "bootstrap and the worst-case p over the three shared nulls."
            ),
            "covariance": (
                "complete A_plus/A_minus covariance inside each context; zero between the "
                "revealed P253 counters and the future independent Gaussian counters"
            ),
            "claim_boundary": (
                "A context gain distinguishes one shared generator from a context-dependent "
                "effective generator. It does not by itself prove path/state memory."
            ),
        },
        "scientific_card": [
            "MECHANISM QUESTION: can one low-dimensional generator transport the same local pivotal rows through cover and radius contexts?",
            "SEMANTIC VERDICT: no PR277 global-rank row is identical to a P253 local landing-H4 row, so no current numerical rectangle exists.",
            "SELECTED ROWS: N425 A_plus ordinary and A_minus matching-odd, chosen jointly with the full annulus covariance.",
            "MISSING CELLS: the same fixed-p root-toggle A_plus/A_minus readout along one norm-4 Gaussian cover chain.",
            "DECISION: compare shared lambda with context-specific lambda only after the missing cells pass the general-period runner preflight.",
        ],
    }


def render_markdown(report: dict) -> str:
    annulus = report["annulus_row_selection"]
    lines = [
        "# Gaussian x annulus semantic crosswalk",
        "",
        "## Verdict",
        "",
        "No revealed PR277 row and P253 row has the same source, observer, normalization, and transfer coordinate. "
        "There is therefore no numerical 2x2 rectangle to score from the old archives.",
        "",
        "| archive row | source/readout | sector | transfer coordinate | cross-context eligible? |",
        "|---|---|---|---|---|",
        "| PR277 scalar U | global threshold-rank thermal derivative | matching-odd scalar side view | cover generation | no |",
        "| PR277 r2-r6 | global Hermite-Krawtchouk derivative jet | matching-odd | cover generation | no |",
        "| P253 A_plus | fixed-p root-toggle landing H4 | ordinary/matching-even | log2(R/2) | target row |",
        "| P253 A_minus | fixed-p root-toggle landing H4 | matching-odd | log2(R/2) | target row |",
        "",
        "Scaling or renaming cannot turn a product-measure derivative into a conditional root-toggle observable.",
        "",
        "## Frozen row choice",
        "",
        f"The existing annulus context selected `{annulus['selected_existing_annulus_context']}` for the correlated "
        "`(A_plus,A_minus)` pair. Its minimum adjacent-lambda Mahalanobis separation is "
        f"`{annulus['candidates'][annulus['selected_existing_annulus_context']]['joint_minimum_adjacent_mahalanobis_squared']:.6g}`. "
        "This is weak; it freezes the row/geometry, not a mechanism conclusion.",
        "",
        "## Missing cells",
        "",
        "| context | A_plus ordinary | A_minus matching-odd |",
        "|---|---|---|",
        "| annulus radius doubling | existing P253 N425 | existing P253 N425 |",
        "| Gaussian cover doubling | missing | missing |",
        "",
        "The missing cells must use the P253 fixed-p root-toggle and landing-shell definitions verbatim on a norm-4 cover chain. "
        "The existing cyclic multiradius runner accepts the primitive N85/N170 parents but rejects the nonprimitive N340/N680 children, "
        "so the frozen acquisition requires only a general-period geometry adapter, not a new observable framework.",
        "",
        "## Frozen score",
        "",
        "The shared model uses one of `lambda={0,1/2,1}` in both contexts and both rows. "
        "The minimal context-enriched adversary allows one lambda per context while keeping it shared across A_plus/A_minus. "
        "All fixed scores and the predeclared bootstrap comparison use complete within-context covariance; the future Gaussian block is independent of P253.",
        "",
        "## Scientific card",
        "",
    ]
    lines.extend(f"- {row}" for row in report["scientific_card"])
    lines.append("")
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gaussian-score", type=Path, required=True)
    parser.add_argument("--annulus-score", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    report = build_report(args.gaussian_score, args.annulus_score)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(report), encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
