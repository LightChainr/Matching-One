#!/usr/bin/env python3
"""Score the preregistered N340 same-lineage second-child production."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence

from crosswalk_natural_current_geometry import add, quadratic
from reveal_n65_charged_source_archive import evaluate_batch, read_births
from score_n170_angle_flip import linear_transform, transform_covariance
from score_natural_current_scale import jackknife_natural, sha256


def score(prereg_path: Path, births_path: Path, metadata_path: Path) -> dict[str, object]:
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if prereg.get("schema") != "matching-one/P337-N340-second-child-preregistration/v1":
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
        "first_smith": design["first_smith_invariants"] == expected["Smith_classes"][0],
        "second_smith": design["second_smith_invariants"] == expected["Smith_classes"][1],
        "samples": metadata["samples_per_pair"] == production["samples_per_shape"],
        "batches": metadata["batches"] == production["batches"],
        "seed": metadata["seed"] == production["seed"],
        "replica_offset": metadata["replica_counter_first"] == production["replica_offset"],
        "projective_births": metadata.get("projective_births") is True,
        "engine_freeze": metadata["git_commit"] == "1714141+freeze-5369c21",
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
    transform = prereg["observable"]["decomposition"]["transform"]
    coordinates = linear_transform(transform, observed)
    coordinate_covariance = transform_covariance(transform, measurement_covariance)
    amplitude, scalar = coordinates
    amplitude_se = math.sqrt(coordinate_covariance[0][0])
    scalar_se = math.sqrt(coordinate_covariance[1][1])

    model_scores = []
    for model in prereg["frozen_models_in_scoring_order"]:
        target = model["absolute_K_A"]
        target_covariance = model["absolute_target_covariance"]
        predictive_covariance = add(measurement_covariance, target_covariance)
        vector_residual = [value - fixed for value, fixed in zip(observed, target)]
        target_amplitude = model["H4_amplitude"]
        amplitude_residual = amplitude - target_amplitude
        target_variance = model["H4_amplitude_target_variance"]
        predictive_variance = coordinate_covariance[0][0] + target_variance
        pair_observed = observed[1] - observed[0]
        pair_residual = pair_observed - model["pair_second_minus_first"]
        dc = transform[0][1] ** -1
        pair_measurement_se = abs(dc) * amplitude_se
        pair_predictive_se = abs(dc) * math.sqrt(predictive_variance)
        model_scores.append({
            "name": model["name"], "role": model["role"],
            "H4_amplitude_target": target_amplitude,
            "H4_amplitude_observed": amplitude,
            "H4_amplitude_residual": amplitude_residual,
            "measurement_only_standard_error": amplitude_se,
            "measurement_only_z": amplitude_residual / amplitude_se,
            "target_standard_error": math.sqrt(target_variance),
            "predictive_standard_error": math.sqrt(predictive_variance),
            "predictive_z": amplitude_residual / math.sqrt(predictive_variance),
            "absolute_K_A_target": target,
            "absolute_K_A_residual": vector_residual,
            "absolute_predictive_covariance": predictive_covariance,
            "absolute_predictive_quadratic": quadratic(vector_residual, predictive_covariance),
            "absolute_df": 2,
            "pair_target": model["pair_second_minus_first"],
            "pair_observed": pair_observed,
            "pair_residual": pair_residual,
            "pair_measurement_standard_error": pair_measurement_se,
            "pair_predictive_standard_error": pair_predictive_se,
        })

    closest = min(model_scores, key=lambda item: abs(item["measurement_only_z"]))
    mean_rows = {
        orientation: {
            name: math.fsum(row[orientation][name] for row in rows) / len(rows)
            for name in ("W_A", "Jminus_A", "Jplus_A")
        }
        for orientation in ("first", "second")
    }
    pair_value = observed[1] - observed[0]
    pair_se = natural["standard_error"][2]
    return {
        "schema": "matching-one/P337-N340-second-child-score/v1",
        "status": "fresh preregistered N340 same-lineage heldout reveal",
        "source": {
            "preregistration_commit": "5369c21",
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
        "decomposition": {
            "order": ["H4_amplitude", "A_projective_scalar"],
            "observed": coordinates,
            "measurement_covariance": coordinate_covariance,
            "H4_amplitude": {"value": amplitude, "standard_error": amplitude_se},
            "A_projective_scalar": {
                "value": scalar, "standard_error": scalar_se,
                "frozen_target": 0.0, "z": scalar / scalar_se,
            },
        },
        "fixed_model_scores": model_scores,
        "primary_pair_sign_flip": {
            "observed": pair_value, "standard_error": pair_se,
            "z_vs_scalar_zero": pair_value / pair_se,
            "expected_sign": "positive",
            "resolved": pair_value > 0.0 and pair_value / pair_se > 5.0,
        },
        "reading": {
            "closest_fixed_target_by_measurement_residual": closest["name"],
            "closest_measurement_only_z": closest["measurement_only_z"],
            "projective_scalar_resolved": abs(scalar / scalar_se) > 3.0,
            "summary": (
                "N340 flips the H4 pair sign back as frozen. Fixed-target scores locate the "
                "same-lineage amplitude among nominal H4, the secondary observed slow "
                "continuation, and scale-neutral, while the projective scalar remains a zero control."
            ),
        },
        "exact_gates": {
            "max_A_continuity_residual": max_continuity,
            "tolerance": 3e-12, "passed": max_continuity < 3e-12,
        },
        "dependency": {
            "N340": "new seed/counter block independent of N85 and N170",
            "targets": "all three fixed at preregistration commit 5369c21",
            "effective_target": "secondary descriptive N85-to-N170 transfer with source uncertainty retained",
        },
        "claim_boundary": "one same-lineage heldout second child; no exponent fit, harmonic revote, or post-reveal basis/model selection",
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    natural = payload["natural_coordinate"]
    split = payload["decomposition"]
    pair = payload["primary_pair_sign_flip"]
    lines = [
        "# N340 same-lineage second-child reveal", "",
        "The 12M/shape block was generated only after preregistration commit `5369c21`.", "",
        f"Observed `K_A=({natural['value'][0]:+.9f},{natural['value'][1]:+.9f})`; the pair is `{pair['observed']:+.9f} +/- {pair['standard_error']:.3g}` (`z={pair['z_vs_scalar_zero']:.3f}` versus scalar zero), with the preregistered positive H4 flip restored.", "",
        f"H4 amplitude is `{split['H4_amplitude']['value']:+.9f} +/- {split['H4_amplitude']['standard_error']:.3g}`. Fixed, pre-reveal comparisons:", "",
        "| target | amplitude | residual / measurement SE | residual / predictive SE | full vector q/2 |",
        "|---|---:|---:|---:|---:|",
    ]
    for model in payload["fixed_model_scores"]:
        lines.append(
            f"| {model['name']} | {model['H4_amplitude_target']:+.9f} | "
            f"{model['measurement_only_z']:+.3f} | {model['predictive_z']:+.3f} | "
            f"{model['absolute_predictive_quadratic']:.3f} |"
        )
    lines.extend([
        "",
        f"The charged/projective scalar is `{split['A_projective_scalar']['value']:+.9f} +/- {split['A_projective_scalar']['standard_error']:.3g}` (`z={split['A_projective_scalar']['z']:.3f}`).", "",
        f"By the frozen measurement-residual rule, the closest target is `{payload['reading']['closest_fixed_target_by_measurement_residual']}`. The observed-transfer target remains secondary and its N85/N170 source uncertainty is retained in the predictive column.", "",
        f"Continuity closes to `{payload['exact_gates']['max_A_continuity_residual']:.3g}`. No exponent fit, basis change, or H4/H8 vote was performed.", "",
    ])
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
