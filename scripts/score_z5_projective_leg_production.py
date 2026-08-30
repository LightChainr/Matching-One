#!/usr/bin/env python3
"""Fresh P250 production scorer with a hard support-before-phase lock."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from score_z5_charged_multiseparation import (
    closure,
    closure_score,
    denominator_rows,
    jackknife,
    read_batches,
)
from score_z5_charged_threepoint import zero_score
from score_z5_projective_leg import covariance_nondegenerate
from z5_projective_leg_multiseparation_mc import SCHEMA


TARGET_SEPARATIONS = (1, 2)


def covariance_of_jackknife(values: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(values)
    center = [sum(row[index] for row in values) / count for index in range(len(values[0]))]
    factor = (count - 1) / count
    return [
        [
            factor * sum(
                (row[first] - center[first]) * (row[second] - center[second])
                for row in values
            )
            for second in range(len(center))
        ]
        for first in range(len(center))
    ]


def support_then_phase(
    payload: dict,
    batches: Sequence[dict],
    *,
    minimum_pair_z: float,
    support_alpha: float,
) -> dict:
    support_rows = {}
    support_pass = True
    for separation in TARGET_SEPARATIONS:
        denominators = denominator_rows(batches, separation)
        pair_ready = all(value["abs_z"] >= minimum_pair_z for value in denominators.values())
        local_point, local_covariance, _ = jackknife(batches, separation, "local_variance")
        nondegenerate, reason = covariance_nondegenerate(local_covariance)
        support = zero_score(local_point, local_covariance) if nondegenerate else None
        detected = bool(support and support["survival_p"] < support_alpha)
        passed = pair_ready and nondegenerate and detected
        support_pass = support_pass and passed
        support_rows[str(separation)] = {
            "two_point_denominators": denominators,
            "minimum_two_point_abs_z": min(value["abs_z"] for value in denominators.values()),
            "denominator_gate_abs_z": minimum_pair_z,
            "all_four_denominators_pass": pair_ready,
            "local_variance_normalized_cubic_point": local_point,
            "local_variance_normalized_cubic_covariance": local_covariance,
            "cubic_covariance_nondegenerate": nondegenerate,
            "cubic_covariance_reason": reason,
            "support_zero_score": support,
            "support_alpha": support_alpha,
            "support_detected": detected,
            "support_stage_pass": passed,
        }
    if not support_pass:
        phase = {
            "status": "locked_support_gate_failed",
            "computed": False,
            "reason": "phase code path is not entered unless denominator, covariance and support gates pass at both d1 and d2",
        }
    else:
        points = {}
        delete_rows = {}
        per_separation = {}
        for separation in TARGET_SEPARATIONS:
            point, _covariance, deleted = jackknife(batches, separation, "separation")
            points[separation] = point
            delete_rows[separation] = deleted
            per_separation[str(separation)] = closure_score(point, deleted)
        joint_point = []
        for separation in TARGET_SEPARATIONS:
            value = closure(points[separation])
            joint_point.extend((value.real, value.imag))
        joint_delete = []
        for omitted in range(len(batches)):
            row = []
            for separation in TARGET_SEPARATIONS:
                value = closure(delete_rows[separation][omitted])
                row.extend((value.real, value.imag))
            joint_delete.append(row)
        joint_covariance = covariance_of_jackknife(joint_delete)
        joint_score = zero_score(joint_point, joint_covariance)
        phase = {
            "status": "unlocked_after_support_pass",
            "computed": True,
            "relation": "Omega113_plus*Omega122_minus-Omega113_minus*Omega122_plus",
            "joint_order": ["d1_re", "d1_im", "d2_re", "d2_im"],
            "joint_point": joint_point,
            "joint_covariance": joint_covariance,
            "joint_zero_score": joint_score,
            "per_separation": per_separation,
        }
    return {
        "support_stage": support_rows,
        "support_gate_passed": support_pass,
        "phase_closure": phase,
    }


def validate(payload: dict, manifest: Mapping[str, object]) -> None:
    if payload.get("schema") != SCHEMA or not payload["exact_gate"]["passed"]:
        raise ValueError("wrong or failed projective-leg response")
    expected = manifest["run"]
    observed = payload["run"]
    for key in (
        "samples", "batches", "workers", "p", "seed", "replica_offset",
        "replica_last_exclusive",
    ):
        if observed.get(key) != expected.get(key):
            raise ValueError(f"run differs from manifest for {key}")
    if manifest.get("phase_policy") != "locked_until_both_support_stages_pass":
        raise ValueError("manifest does not freeze the phase lock")


def score(payload: dict, batches: Sequence[dict], manifest: Mapping[str, object]) -> dict:
    validate(payload, manifest)
    result = support_then_phase(
        payload,
        batches,
        minimum_pair_z=float(manifest["gates"]["minimum_pair_abs_z"]),
        support_alpha=float(manifest["gates"]["support_alpha"]),
    )
    return {
        "schema": "matching-one/z5-projective-leg-fresh-production-score/v1",
        "status": "fresh_production_reveal",
        "run": payload["run"],
        "score_order": [
            "d1_d2_denominator_gate",
            "d1_d2_8real_cubic_support_gate",
            "conditional_joint_phase_closure",
        ],
        **result,
        "decision": (
            "support_confirmed_phase_revealed" if result["support_gate_passed"]
            else "support_not_confirmed_phase_locked"
        ),
        "claim_boundary": [
            "The support gate and phase lock were frozen before the fresh stream.",
            "Phase closure is downstream of nonzero cubic support and is not a support vote.",
            "The projective-leg row is topological; this score does not identify a local primary or universal OPE coefficient.",
        ],
    }


def render(result: Mapping[str, object]) -> str:
    lines = [
        "# P250 fresh projective-leg production", "",
        "| d | minimum pair z | pair gate | cubic support chi2/8 | p | support pass |",
        "|---:|---:|---|---:|---:|---|",
    ]
    for separation, row in result["support_stage"].items():
        support = row["support_zero_score"]
        lines.append(
            f"| {separation} | {row['minimum_two_point_abs_z']:.3f} | "
            f"{row['all_four_denominators_pass']} | {support['chi_square']:.6g} | "
            f"{support['survival_p']:.6g} | {row['support_stage_pass']} |"
        )
    phase = result["phase_closure"]
    lines.extend(["", f"Support decision: `{result['decision']}`.", ""])
    if phase["computed"]:
        joint = phase["joint_zero_score"]
        lines.extend([
            "Phase closure was unlocked only after both support stages passed.", "",
            f"Joint d1/d2 closure: `{joint['chi_square']}/{joint['degrees_of_freedom']}`, p `{joint['survival_p']}`.", "",
        ])
    else:
        lines.extend([f"Phase closure: `{phase['status']}`; it was not computed.", ""])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("response", type=Path)
    parser.add_argument("batches", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = score(
        json.loads(args.response.read_text()),
        read_batches(args.batches),
        json.loads(args.manifest.read_text()),
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
