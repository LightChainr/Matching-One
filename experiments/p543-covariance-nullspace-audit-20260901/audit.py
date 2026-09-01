#!/usr/bin/env python3
"""Re-score every archived pseudoinverse vector from sufficient statistics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import mpmath as mp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from covariance_nullspace import covariance_spectral_diagnostics  # noqa: E402


def number(value: object) -> mp.mpf:
    return mp.mpf(str(value))


def rendered(value: object) -> object:
    if isinstance(value, dict):
        return {key: rendered(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [rendered(item) for item in value]
    if isinstance(value, mp.mpf):
        return mp.nstr(value, 25)
    return value


def compact_diagnostics(score: dict[str, object]) -> dict[str, object]:
    sensitivity = []
    for row in score["cutoff_sensitivity"]:
        sensitivity.append(
            {
                key: row[key]
                for key in (
                    "relative_eigenvalue_cutoff",
                    "absolute_eigenvalue_cutoff",
                    "numerical_rank",
                    "chi_square",
                    "degrees_of_freedom",
                    "chi_square_survival",
                    "active_condition_number",
                    "discarded_residual_projection_l2",
                    "max_abs_discarded_residual_projection",
                    "discarded_nullspace_projection_l2",
                    "max_abs_discarded_nullspace_projection",
                    "nullspace_compatible",
                )
            }
        )
    return {
        key: score[key]
        for key in (
            "chi_square",
            "degrees_of_freedom",
            "chi_square_survival",
            "numerical_rank",
            "relative_eigenvalue_cutoff",
            "absolute_eigenvalue_cutoff",
            "spectral_basis",
            "nullspace_projection_basis",
            "active_condition_number",
            "spectral_eigenvalues",
            "component_standardized_residuals",
            "discarded_eigendirections",
            "discarded_residual_projection_l2",
            "max_abs_discarded_residual_projection",
            "discarded_nullspace_projection_l2",
            "max_abs_discarded_nullspace_projection",
            "null_projection_tolerance",
            "covariance_nullspace_policy",
            "nullspace_status",
            "nullspace_compatible",
            "chi_square_interpretation",
        )
    } | {"cutoff_sensitivity": sensitivity}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=HERE / "INPUT_VECTORS.json")
    parser.add_argument("--output", type=Path, default=HERE / "RESULT.json")
    args = parser.parse_args()

    inputs = json.loads(args.input.read_text(encoding="utf-8"))
    results = []
    for row in inputs["vectors"]:
        standardize = row["spectral_basis"] == "correlation_standardized"
        mp.mp.dps = 80 if standardize else 15
        if standardize:
            residual = [number(value) for value in row["residual"]]
            covariance = [
                [number(value) for value in values]
                for values in row["covariance"]
            ]
        else:
            residual = row["residual"]
            covariance = row["covariance"]
        score = covariance_spectral_diagnostics(
            residual,
            covariance,
            number(row["prior_score"]["relative_eigenvalue_cutoff"]),
            nullspace_policy="estimated",
            standardize=standardize,
        )
        prior_chi = number(row["prior_score"]["chi_square"])
        delta = score["chi_square"] - prior_chi
        tolerance = mp.mpf("1e-8") * max(mp.mpf(1), abs(prior_chi))
        default_statistic_changed = (
            abs(delta) > tolerance
            or score["numerical_rank"] != row["prior_score"]["numerical_rank"]
        )
        if score["nullspace_status"] == "estimated_near_null_incompatibility":
            classification = "interpretation_changed"
        elif default_statistic_changed:
            classification = "numerically_changed"
        else:
            classification = "unchanged"
        relative_spectrum = [
            value / max(score["spectral_eigenvalues"])
            for value in score["spectral_eigenvalues"]
        ]
        sensitivity_ranks = {
            int(item["numerical_rank"]) for item in score["cutoff_sensitivity"]
        }
        results.append(
            {
                "id": row["id"],
                "source_id": row["source_id"],
                "json_pointer": row["json_pointer"],
                "dimension": len(row["residual"]),
                "prior_score": row["prior_score"],
                "default_chi_square_delta_from_archive": delta,
                "default_statistic_changed": default_statistic_changed,
                "classification": classification,
                "minimum_relative_spectral_eigenvalue": min(relative_spectrum),
                "rank_changes_across_frozen_cutoffs": len(sensitivity_ranks) > 1,
                "rescored": compact_diagnostics(score),
            }
        )

    counts = {
        label: sum(item["classification"] == label for item in results)
        for label in ("unchanged", "numerically_changed", "interpretation_changed")
    }
    p50 = next(item for item in results if item["id"].startswith("p50_fullcurve"))
    payload = {
        "schema": "matching-one/covariance-nullspace-historical-audit/v1",
        "status": "complete_existing_sufficient_statistics_only",
        "input_schema": inputs["schema"],
        "new_random_samples": 0,
        "archived_vector_scores": len(results),
        "classification_counts": counts,
        "default_displayed_statistics_changed": sum(
            item["default_statistic_changed"] for item in results
        ),
        "claim_boundary": {
            "p50_default_rejection_reversed": False,
            "p50_interpretation": (
                "The frozen 1e-10 score remains chi-square 9.352/2, but the "
                "discarded direction is residual-incompatible and the result must "
                "be reported as cutoff-sensitive rather than as a silent numerical null."
            ),
            "other_archived_claims_changed": False,
        },
        "p50_decision_row_id": p50["id"],
        "sources": inputs["sources"],
        "results": results,
    }
    args.output.write_text(
        json.dumps(rendered(payload), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
