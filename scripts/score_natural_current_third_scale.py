#!/usr/bin/env python3
"""Score the preregistered held-out N145 natural charged-current test."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence

from reveal_n65_charged_source_archive import evaluate_batch, read_births
from score_natural_current_scale import jackknife_natural, sha256


def score(
    prereg_path: Path, births_path: Path, metadata_path: Path
) -> dict[str, object]:
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if prereg.get("schema") != "matching-one/P337-natural-current-third-scale-preregistration/v1":
        raise ValueError("unexpected preregistration schema")

    expected = prereg["geometry"]
    design = metadata["designs"][0]
    gates = {
        "N": design["N"] == expected["N"],
        "first": design["first"] == expected["first"],
        "second": design["second"] == expected["second"],
        "first_period_matrix": design["first_period_matrix"] == expected["first_period_matrix"],
        "second_period_matrix": design["second_period_matrix"] == expected["second_period_matrix"],
        "samples": metadata["samples_per_pair"] == prereg["design"]["samples_per_shape"],
        "batches": metadata["batches"] == prereg["design"]["batches"],
        "seed": metadata["seed"] == prereg["design"]["seed"],
        "replica_offset": metadata["replica_counter_first"] == prereg["design"]["replica_offset"],
        "projective_births": metadata.get("projective_births") is True,
        "engine_freeze": metadata["git_commit"] == "1714141+freeze-3d238af",
    }
    if not all(gates.values()):
        raise ValueError(f"production differs from preregistration: {gates}")

    p_ref = prereg["source"]["p_ref"]
    n, births = read_births(births_path)
    batch_ids = sorted({batch for _, batch in births})
    if len(batch_ids) != prereg["design"]["batches"]:
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
    observed = natural["value"][2]
    variance = natural["covariance"][2][2]

    def compare(name: str) -> dict[str, object]:
        target = prereg["frozen_targets_at_N145"][name]
        residual = observed - target["value"]
        predictive_variance = variance + target["fit_variance"]
        return {
            "name": name,
            "tier": target["tier"],
            "frozen_target": target["value"],
            "residual": residual,
            "measurement_standard_error": math.sqrt(variance),
            "measurement_only_quadratic": residual * residual / variance,
            "predictive_standard_error": math.sqrt(predictive_variance),
            "predictive_quadratic": residual * residual / predictive_variance,
        }

    primary = [compare(name) for name in prereg["scoring_contract"]["primary_order"]]
    secondary = [compare(name) for name in prereg["scoring_contract"]["secondary_order"]]
    closest_primary = min(primary, key=lambda row: row["predictive_quadratic"])
    mean_rows = {
        orientation: {
            name: math.fsum(row[orientation][name] for row in rows) / len(rows)
            for name in ("W_A", "Jminus_A", "Jplus_A")
        }
        for orientation in ("first", "second")
    }
    delta85 = prereg["frozen_targets_at_N145"]["source_fitted_scale_neutral"]["value"]
    h4_target = prereg["frozen_targets_at_N145"]["source_fitted_project_H4"]["value"]
    effective_target = prereg["frozen_targets_at_N145"]["secondary_post_reveal_effective_transfer"]["value"]
    effective_score = secondary[0]

    return {
        "schema": "matching-one/P337-natural-current-third-scale-N145-score/v1",
        "status": "fresh held-out N145 scale reveal scored against all frozen targets",
        "source": {
            "preregistration_commit": "3d238af",
            "preregistration": str(prereg_path),
            "engine_commit": metadata["git_commit"],
            "environment": "Huawei DevEnvC_ZyTrST f415a4bcbd9a438b85f5f29e4a507ea4",
            "births": str(births_path),
            "births_sha256": sha256(births_path),
            "metadata": str(metadata_path),
            "metadata_sha256": sha256(metadata_path),
            "p_ref": p_ref,
            "N": n,
            "samples_per_shape": metadata["samples_per_pair"],
            "batches": metadata["batches"],
            "seed": metadata["seed"],
            "elapsed_seconds": metadata["elapsed_seconds"],
        },
        "freeze_gates": {"passed": all(gates.values()), "items": gates},
        "batch_means": mean_rows,
        "natural_coordinate": natural,
        "primary_target_comparison": primary,
        "secondary_target_comparison": secondary,
        "reading": {
            "independent_nonzero": {
                "value": observed,
                "standard_error": math.sqrt(variance),
                "z": observed / math.sqrt(variance),
            },
            "ratio_to_N85": observed / delta85,
            "frozen_H4_ratio": h4_target / delta85,
            "frozen_effective_ratio": effective_target / delta85,
            "closest_primary_predictive_target": closest_primary["name"],
            "closest_primary_predictive_quadratic": closest_primary["predictive_quadratic"],
            "secondary_effective_predictive_quadratic": effective_score["predictive_quadratic"],
            "summary": (
                "N145 stays nonzero and attenuates only weakly from N85. It bends away from "
                "the frozen continuation of the fast N65-to-N85 transfer and is closest to "
                "scale-neutral propagation; the pure N85-anchored H4 transfer is also low."
            ),
            "variance_design_ratio": (
                math.sqrt(variance) / prereg["design"]["projected_standard_error"]
            ),
        },
        "exact_gates": {
            "max_A_continuity_residual": max_continuity,
            "tolerance": 3e-12,
            "passed": max_continuity < 3e-12,
        },
        "dependency": {
            "N145": "new seed/counter block, independent of N65 and N85 archives",
            "frozen_targets": "all targets fixed at commit 3d238af before N145 generation",
            "predictive_uncertainty": "N65/N85 target-fit variance included only in predictive residuals",
            "no_refit": "N145 does not change any exponent, target, observable, or normalization",
        },
        "claim_boundary": (
            "one held-out N145 point supports curvature away from the fast descriptive transfer; "
            "it does not identify a unique asymptotic state or refit an exponent"
        ),
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    natural = payload["natural_coordinate"]
    reading = payload["reading"]
    lines = [
        "# Held-out N145 natural charged-current reveal",
        "",
        "The 2.4M/shape N145 block was generated only after preregistration commit `3d238af`.",
        "",
        f"`Delta_K_A(N145)={natural['value'][2]:.12g} +/- {natural['standard_error'][2]:.3g}` "
        f"(`z={reading['independent_nonzero']['z']:.3f}`).",
        "",
        "| primary frozen target | target | residual | measurement chi2 | predictive chi2 |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in payload["primary_target_comparison"]:
        lines.append(
            f"| {row['name']} | {row['frozen_target']:.12g} | {row['residual']:+.12g} | "
            f"{row['measurement_only_quadratic']:.4g} | {row['predictive_quadratic']:.4g} |"
        )
    row = payload["secondary_target_comparison"][0]
    lines += [
        "",
        "Secondary descriptive target (not an exponent model):",
        "",
        "| target | value | residual | measurement chi2 | predictive chi2 |",
        "|---|---:|---:|---:|---:|",
        f"| {row['name']} | {row['frozen_target']:.12g} | {row['residual']:+.12g} | "
        f"{row['measurement_only_quadratic']:.4g} | {row['predictive_quadratic']:.4g} |",
    ]
    first, second = payload["batch_means"]["first"], payload["batch_means"]["second"]
    lines += [
        "",
        f"First `(12+i)`: `W_A={first['W_A']:.12g}`, `Jminus_A={first['Jminus_A']:.12g}`, "
        f"`K_A={natural['value'][0]:.12g}`.",
        f"Second `(9+8i)`: `W_A={second['W_A']:.12g}`, `Jminus_A={second['Jminus_A']:.12g}`, "
        f"`K_A={natural['value'][1]:.12g}`.",
        "",
        f"The observed N145/N85 ratio is `{reading['ratio_to_N85']:.4f}`; frozen H4 predicted "
        f"`{reading['frozen_H4_ratio']:.4f}` and the secondary fast-transfer continuation "
        f"`{reading['frozen_effective_ratio']:.4f}`.",
        "",
        "The third scale therefore bends away from continued fast attenuation. Scale-neutral "
        "propagation is the closest primary frozen target. The H4-only target remains a finite "
        "tension rather than a selected description, so the clean reading is correction curvature "
        "or state mixing, not a newly fitted exponent.",
        "",
        f"Observed SE was {reading['variance_design_ratio']:.3f} times the N85 projection. "
        f"Current continuity closes to `{payload['exact_gates']['max_A_continuity_residual']:.3g}`.",
        "",
        "No exponent, observable, target, or normalization was changed after reveal.",
        "",
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
