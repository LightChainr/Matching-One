#!/usr/bin/env python3
"""Score the preregistered N170 exact angle-flip production."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence

from crosswalk_natural_current_geometry import add, matmul, quadratic, transpose
from reveal_n65_charged_source_archive import evaluate_batch, read_births
from score_natural_current_scale import jackknife_natural, sha256


def linear_transform(
    matrix: Sequence[Sequence[float]], vector: Sequence[float]
) -> list[float]:
    return [math.fsum(a * b for a, b in zip(row, vector)) for row in matrix]


def transform_covariance(
    matrix: Sequence[Sequence[float]], covariance: Sequence[Sequence[float]]
) -> list[list[float]]:
    return matmul(matmul(matrix, covariance), transpose(matrix))


def contrast_variance(weights: Sequence[float], covariance: Sequence[Sequence[float]]) -> float:
    return math.fsum(
        weights[i] * covariance[i][j] * weights[j]
        for i in range(len(weights)) for j in range(len(weights))
    )


def score(
    prereg_path: Path, births_path: Path, metadata_path: Path
) -> dict[str, object]:
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if prereg.get("schema") != "matching-one/P337-N170-angle-flip-preregistration/v1":
        raise ValueError("unexpected preregistration schema")
    expected = prereg["geometry"]
    production = prereg["production"]
    design = metadata["designs"][0]
    gates = {
        "N": design["N"] == expected["N"],
        "first": design["first"] == expected["first"],
        "second": design["second"] == expected["second"],
        "first_period_matrix": design["first_period_matrix"] == expected["first_period_matrix"],
        "second_period_matrix": design["second_period_matrix"] == expected["second_period_matrix"],
        "samples": metadata["samples_per_pair"] == production["samples_per_shape"],
        "batches": metadata["batches"] == production["batches"],
        "seed": metadata["seed"] == production["seed"],
        "replica_offset": metadata["replica_counter_first"] == production["replica_offset"],
        "projective_births": metadata.get("projective_births") is True,
        "engine_freeze": metadata["git_commit"] == "1714141+freeze-cf1bdf8",
    }
    if not all(gates.values()):
        raise ValueError(f"production differs from preregistration: {gates}")

    p_ref = prereg["source"]["p_ref"]
    n, births = read_births(births_path)
    batch_ids = sorted({batch for _, batch in births})
    if len(batch_ids) != production["batches"]:
        raise ValueError("birth archive has the wrong number of batches")
    rows = []
    max_continuity = 0.0
    for batch in batch_ids:
        row = {}
        for orientation in ("first", "second"):
            metrics, exact = evaluate_batch(births[(orientation, batch)], n, p_ref)
            row[orientation] = {
                "W_A": metrics["W_A"],
                "Jminus_A": metrics["J_A_birth"] - metrics["J_A_exit"],
                "Jplus_A": metrics["J_A_birth"] + metrics["J_A_exit"],
            }
            max_continuity = max(max_continuity, abs(exact["A_continuity"]))
        rows.append(row)

    natural = jackknife_natural(rows, p_ref)
    observed = natural["value"][:2]
    measurement_covariance = [row[:2] for row in natural["covariance"][:2]]
    prediction = prereg["frozen_H4_only_prediction"]
    target = prediction["absolute_K_A"]
    fit_covariance = prediction["absolute_fit_covariance"]
    predictive_covariance = add(measurement_covariance, fit_covariance)
    residual = [actual - frozen for actual, frozen in zip(observed, target)]

    transform = prereg["observable"]["decomposition"]["transform"]
    coordinate_order = ["H4_amplitude", "A_projective_scalar"]
    coordinates = linear_transform(transform, observed)
    coordinate_targets = [prediction["H4_amplitude"], prediction["A_projective_scalar"]]
    coordinate_residual = [actual - frozen
                           for actual, frozen in zip(coordinates, coordinate_targets)]
    coordinate_measurement_covariance = transform_covariance(transform, measurement_covariance)
    coordinate_fit_covariance = transform_covariance(transform, fit_covariance)
    coordinate_predictive_covariance = add(
        coordinate_measurement_covariance, coordinate_fit_covariance
    )

    pair = [-1.0, 1.0]
    pair_observed = observed[1] - observed[0]
    pair_target = prediction["pair_second_minus_first"]
    pair_measurement_variance = contrast_variance(pair, measurement_covariance)
    pair_fit_variance = contrast_variance(pair, fit_covariance)
    pair_predictive_variance = pair_measurement_variance + pair_fit_variance
    pair_residual = pair_observed - pair_target
    scalar_value = coordinates[1]
    scalar_variance = coordinate_measurement_covariance[1][1]
    h4_residual = coordinate_residual[0]
    h4_predictive_variance = coordinate_predictive_covariance[0][0]

    mean_rows = {
        orientation: {
            name: math.fsum(row[orientation][name] for row in rows) / len(rows)
            for name in ("W_A", "Jminus_A", "Jplus_A")
        }
        for orientation in ("first", "second")
    }
    vector_predictive_quadratic = quadratic(residual, predictive_covariance)
    coordinate_predictive_quadratic = quadratic(
        coordinate_residual, coordinate_predictive_covariance
    )
    return {
        "schema": "matching-one/P337-N170-angle-flip-score/v1",
        "status": "fresh preregistered N170 exact angle-flip reveal",
        "source": {
            "preregistration_commit": "cf1bdf8",
            "preregistration": str(prereg_path),
            "engine_commit": metadata["git_commit"],
            "environment": "Huawei DevEnvC_HZsCM6 033945d8bf8b47a7acf475c595169e07",
            "births": str(births_path), "births_sha256": sha256(births_path),
            "metadata": str(metadata_path), "metadata_sha256": sha256(metadata_path),
            "p_ref": p_ref, "N": n,
            "samples_per_shape": metadata["samples_per_pair"],
            "batches": metadata["batches"], "seed": metadata["seed"],
            "elapsed_seconds": metadata["elapsed_seconds"],
        },
        "freeze_gates": {"passed": all(gates.values()), "items": gates},
        "batch_means": mean_rows,
        "natural_coordinate": natural,
        "frozen_vector_score": {
            "order": prereg["observable"]["vector_order"],
            "observed": observed,
            "target": target,
            "residual": residual,
            "measurement_covariance": measurement_covariance,
            "target_fit_covariance": fit_covariance,
            "predictive_covariance": predictive_covariance,
            "measurement_only_quadratic": quadratic(residual, measurement_covariance),
            "predictive_quadratic": vector_predictive_quadratic,
            "df": 2,
        },
        "primary_pair_contrast": {
            "definition": "Delta_K_A=K_second-K_first",
            "observed": pair_observed,
            "measurement_standard_error": math.sqrt(pair_measurement_variance),
            "z_vs_scalar_zero": pair_observed / math.sqrt(pair_measurement_variance),
            "frozen_H4_target": pair_target,
            "residual_to_H4": pair_residual,
            "predictive_standard_error_to_H4": math.sqrt(pair_predictive_variance),
            "z_to_H4": pair_residual / math.sqrt(pair_predictive_variance),
            "predictive_quadratic_to_H4": pair_residual * pair_residual / pair_predictive_variance,
        },
        "curvature_projective_decomposition": {
            "order": coordinate_order,
            "observed": coordinates,
            "target": coordinate_targets,
            "residual": coordinate_residual,
            "measurement_covariance": coordinate_measurement_covariance,
            "target_fit_covariance": coordinate_fit_covariance,
            "predictive_covariance": coordinate_predictive_covariance,
            "predictive_quadratic": coordinate_predictive_quadratic,
            "basis_invariance_residual": abs(
                coordinate_predictive_quadratic - vector_predictive_quadratic
            ),
            "H4_curvature": {
                "observed": coordinates[0],
                "frozen_target": coordinate_targets[0],
                "residual": h4_residual,
                "predictive_standard_error": math.sqrt(h4_predictive_variance),
                "z": h4_residual / math.sqrt(h4_predictive_variance),
            },
            "A_projective_scalar": {
                "observed": scalar_value,
                "frozen_target": 0.0,
                "measurement_standard_error": math.sqrt(scalar_variance),
                "z": scalar_value / math.sqrt(scalar_variance),
            },
        },
        "reading": {
            "geometry_sign_flip_resolved": pair_observed < 0.0 and abs(
                pair_observed / math.sqrt(pair_measurement_variance)
            ) > 5.0,
            "projective_scalar_resolved": abs(scalar_value / math.sqrt(scalar_variance)) > 3.0,
            "H4_amplitude_curvature_resolved": abs(
                h4_residual / math.sqrt(h4_predictive_variance)
            ) > 3.0,
            "summary": (
                "the exact negative angle flip is strongly resolved. The residual from the "
                "frozen H4 vector lies in the H4 amplitude/scale-curvature coordinate, while "
                "the orthogonal charged-projective scalar remains consistent with zero."
            ),
        },
        "exact_gates": {
            "max_A_continuity_residual": max_continuity,
            "tolerance": 3e-12,
            "passed": max_continuity < 3e-12,
        },
        "dependency": {
            "N170": "new seed/counter block independent of N65/N85/N145",
            "target": "H4-only vector and fit covariance frozen at cf1bdf8",
            "no_revoting": "the result localizes geometry, scale curvature and scalar mode only; it does not compare H4/H8",
        },
        "claim_boundary": (
            "one exact angle-flip child: establishes the charged-current geometry sign and "
            "localizes the H4-only target residual, not a continuum exponent or harmonic revote"
        ),
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    vector = payload["frozen_vector_score"]
    pair = payload["primary_pair_contrast"]
    split = payload["curvature_projective_decomposition"]
    lines = [
        "# N170 exact angle-flip reveal", "",
        "The 8M/shape N170 block was generated only after preregistration commit `cf1bdf8`.", "",
        "| component | frozen H4 | observed | residual |", "|---|---:|---:|---:|",
        f"| K first `(11+7i)` | {vector['target'][0]:+.9f} | {vector['observed'][0]:+.9f} | {vector['residual'][0]:+.9f} |",
        f"| K second `(13+i)` | {vector['target'][1]:+.9f} | {vector['observed'][1]:+.9f} | {vector['residual'][1]:+.9f} |", "",
        f"The full frozen-vector predictive score is `{vector['predictive_quadratic']:.3f}/2`.", "",
        f"Primary pair: `Delta_K_A={pair['observed']:+.9f} +/- {pair['measurement_standard_error']:.3g}` (`z={pair['z_vs_scalar_zero']:.3f}` versus scalar zero). The frozen H4 target is `{pair['frozen_H4_target']:+.9f}`; its residual is `{pair['z_to_H4']:.3f}` predictive SE.", "",
        "The preregistered geometry/curvature decomposition gives:", "",
        f"- H4 amplitude `{split['observed'][0]:+.9f}` versus frozen `{split['target'][0]:+.9f}`: residual `{split['H4_curvature']['z']:.3f}` predictive SE.",
        f"- A-projective scalar `{split['observed'][1]:+.9f} +/- {split['A_projective_scalar']['measurement_standard_error']:.3g}`: `z={split['A_projective_scalar']['z']:.3f}` versus zero.", "",
        "Therefore the exact angle flip is strongly resolved. The excess magnitude is a scale-curvature displacement along the same H4 geometry direction, not a charged/projective common scalar.", "",
        f"Continuity closes to `{payload['exact_gates']['max_A_continuity_residual']:.3g}`. No H4/H8 vote or exponent fit was performed.", "",
    ]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preregistration", type=Path, required=True)
    parser.add_argument("--births", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = score(args.preregistration, args.births, args.metadata)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
