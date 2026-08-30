#!/usr/bin/env python3
"""Score the frozen P118 metric-free ratios on the held-out N145 -> N290 block."""
from __future__ import annotations

import argparse
import json
import math
from datetime import datetime
from pathlib import Path
from typing import Mapping, Sequence

from analyze_p48_retrospective import covariance_of_mean
from score_p50_fullcurve_n290 import (
    CHILD_N,
    LINEAGE_SIGN,
    PARENT_N,
    generalized_covariance_score,
    grouped,
    load_metadata,
    pseudovalue_vectors,
    read_one_size,
    rng_group,
    sha256,
    size_statistics,
)


FEATURE_ORDER = ("R_I", "R_T")
SCHEMA = "matching-one/p118-metric-free-ratio-heldout-score/v1"
MANIFEST_SCHEMA = "matching-one/p118-metric-free-ratio-heldout-manifest/v1"


def ratio_coordinates(stat: Mapping[str, float]) -> dict[str, float]:
    """Form the frozen nonlinear coordinates from one full or delete-one estimate."""
    slope = float(stat["mean_slope"])
    p4_s = float(stat["P4_S"])
    p4_d = float(stat["P4_D"])
    denominators = {
        "R_I": p4_s * slope,
        "R_T": p4_d * slope,
    }
    if any(not math.isfinite(value) or value == 0.0 for value in denominators.values()):
        raise ValueError(f"invalid metric-free denominator: {denominators}")
    coordinates = {
        "R_I": float(stat["P4_D_prime"]) / denominators["R_I"],
        "R_T": float(stat["P4_S_prime"]) / denominators["R_T"],
    }
    if any(not math.isfinite(value) for value in coordinates.values()):
        raise ValueError(f"nonfinite metric-free coordinate: {coordinates}")
    return coordinates


def ratio_vector(stat: Mapping[str, float]) -> list[float]:
    coordinates = ratio_coordinates(stat)
    return [coordinates[name] for name in FEATURE_ORDER]


def estimate_size(by_orientation, *, lineage_sign: float) -> dict[str, object]:
    """Jackknife one size, recomputing both nonlinear ratios in every delete-one."""
    point_stat = size_statistics(by_orientation, lineage_sign=lineage_sign)
    point_vector = ratio_vector(point_stat)
    batch_ids = [row.batch for row in by_orientation["first"]]
    deleted = [
        ratio_vector(
            size_statistics(
                by_orientation,
                lineage_sign=lineage_sign,
                omitted=batch,
            )
        )
        for batch in batch_ids
    ]
    pseudo = pseudovalue_vectors(point_vector, deleted)
    covariance = covariance_of_mean(pseudo)
    return {
        "batches": len(batch_ids),
        "coordinates": dict(zip(FEATURE_ORDER, point_vector)),
        "covariance": covariance,
        "source_factors": {
            name: float(point_stat[name])
            for name in (
                "p0",
                "mean_slope",
                "P4_S",
                "P4_D",
                "P4_S_prime",
                "P4_D_prime",
            )
        },
    }


def add_covariances(
    first: Sequence[Sequence[float]], second: Sequence[Sequence[float]]
) -> list[list[float]]:
    if len(first) != len(second) or any(
        len(first[i]) != len(second[i]) for i in range(len(first))
    ):
        raise ValueError("covariance dimensions differ")
    return [
        [float(first[i][j]) + float(second[i][j]) for j in range(len(first[i]))]
        for i in range(len(first))
    ]


def marginal_score(residual: float, variance: float, alpha: float) -> dict[str, object]:
    if not math.isfinite(variance) or variance <= 0.0:
        raise ValueError("marginal residual variance must be positive")
    standard_error = math.sqrt(variance)
    signed_z = residual / standard_error
    p_value = math.erfc(abs(signed_z) / math.sqrt(2.0))
    eliminated = p_value < alpha
    return {
        "residual": residual,
        "variance": variance,
        "standard_error": standard_error,
        "signed_z": signed_z,
        "chi_square": signed_z * signed_z,
        "degrees_of_freedom": 1,
        "p_value": p_value,
        "decision": (
            "eliminated at the frozen alpha for this N145-to-N290 block"
            if eliminated
            else "survives this block; not established as a universal ratio"
        ),
    }


def validate_manifest(manifest_path: Path, root: Path) -> Mapping[str, object]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("wrong P118 manifest schema")
    if manifest.get("status") != "protocol_locked_before_p118_score_execution":
        raise ValueError("P118 protocol is not locked for scoring")
    chronology = manifest["chronology"]
    frozen = datetime.fromisoformat(chronology["definition_freeze_committed_at"])
    revealed = datetime.fromisoformat(chronology["first_target_result_committed_at"])
    if not frozen < revealed:
        raise ValueError("definition freeze does not precede first target result")
    freeze_path = root / chronology["definition_freeze_path"]
    if sha256(freeze_path) != chronology["definition_freeze_sha256"]:
        raise ValueError("definition-freeze digest mismatch")
    expected = {
        "R_I": "P4_D_prime / (P4_S * Mbar_prime)",
        "R_T": "P4_S_prime / (P4_D * Mbar_prime)",
    }
    for name, formula in expected.items():
        if manifest["observables"].get(name) != formula:
            raise ValueError(f"unexpected frozen formula for {name}")
    if tuple(manifest["observables"].get("feature_order", [])) != FEATURE_ORDER:
        raise ValueError("frozen feature order changed")
    for key in (
        "parent_histogram",
        "parent_metadata",
        "child_histogram",
        "child_metadata",
    ):
        item = manifest["inputs"][key]
        path = root / item["path"]
        if sha256(path) != item["sha256"]:
            raise ValueError(f"input digest mismatch for {key}")
    return manifest


def render(root: Path, manifest_path: Path) -> dict[str, object]:
    manifest = validate_manifest(manifest_path, root)
    inputs = manifest["inputs"]
    paths = {
        key: root / inputs[key]["path"]
        for key in (
            "parent_histogram",
            "parent_metadata",
            "child_histogram",
            "child_metadata",
        )
    }
    parent_metadata = load_metadata(paths["parent_metadata"])
    child_metadata = load_metadata(paths["child_metadata"])
    if rng_group(parent_metadata) == rng_group(child_metadata):
        raise ValueError("P118 requires independent parent and child replica domains")
    parent_groups = grouped(
        read_one_size(paths["parent_histogram"], PARENT_N), PARENT_N
    )
    child_groups = grouped(
        read_one_size(paths["child_histogram"], CHILD_N), CHILD_N
    )
    parent = estimate_size(parent_groups, lineage_sign=LINEAGE_SIGN[PARENT_N])
    child = estimate_size(child_groups, lineage_sign=LINEAGE_SIGN[CHILD_N])
    parent_vector = [parent["coordinates"][name] for name in FEATURE_ORDER]
    child_vector = [child["coordinates"][name] for name in FEATURE_ORDER]
    residual = [child_vector[i] - parent_vector[i] for i in range(len(FEATURE_ORDER))]
    covariance = add_covariances(parent["covariance"], child["covariance"])
    alpha = float(manifest["estimation"]["decision_alpha"])
    joint = generalized_covariance_score(
        residual,
        covariance,
        float(manifest["estimation"]["relative_eigenvalue_cutoff"]),
    )
    joint_p = float(joint["chi_square_survival"])
    joint["p_value"] = joint_p
    joint["decision"] = (
        "eliminated at the frozen alpha for this N145-to-N290 block"
        if joint_p < alpha
        else "survives this block; not established as a universal ratio"
    )
    marginals = {
        name: marginal_score(residual[i], covariance[i][i], alpha)
        for i, name in enumerate(FEATURE_ORDER)
    }
    child_over_parent = {
        name: child_vector[i] / parent_vector[i]
        for i, name in enumerate(FEATURE_ORDER)
    }
    return {
        "schema": SCHEMA,
        "status": "frozen held-out score completed",
        "issue": 118,
        "feature_order": list(FEATURE_ORDER),
        "decision_alpha": alpha,
        "chronology": manifest["chronology"],
        "protocol": {
            "manifest": str(manifest_path.relative_to(root)),
            "manifest_sha256": sha256(manifest_path),
            "nonlinear_delete_one_rule": manifest["estimation"]["rule"],
            "covariance_rule": manifest["estimation"]["parent_child_covariance"],
            "scoring_order": manifest["scoring_order"],
        },
        "provenance": {
            key: {
                "path": str(paths[key].relative_to(root)),
                "sha256": sha256(paths[key]),
            }
            for key in paths
        },
        "replica_domains": {
            "parent": list(rng_group(parent_metadata)),
            "child": list(rng_group(child_metadata)),
            "independent": True,
        },
        "observations": {
            f"N{PARENT_N}": parent,
            f"N{CHILD_N}": child,
        },
        "constant_response_null": {
            "residual_definition": "R(N290)-R(N145)",
            "residual": dict(zip(FEATURE_ORDER, residual)),
            "residual_covariance": covariance,
            "joint": joint,
            "marginals": marginals,
        },
        "descriptive_child_over_parent": child_over_parent,
        "interpretation_boundary": (
            "This held-out block reuses P50 raw and is not independent primary evidence. "
            "Survival does not establish cross-microscopic universality; rejection only "
            "rejects a size-constant response over N145 to N290 for this lineage."
        ),
    }


def report_markdown(payload: Mapping[str, object]) -> str:
    obs = payload["observations"]
    null = payload["constant_response_null"]
    joint = null["joint"]
    lines = [
        "# P118 held-out metric-free ratio score",
        "",
        "The definition freeze preceded the first target result by **1 hour 10 minutes 21 seconds**. "
        "Both nonlinear ratios were recomputed inside every size-local delete-one replicate.",
        "",
        "| coordinate | N145 | N290 | N290-N145 | z | marginal p | decision |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for name in FEATURE_ORDER:
        marginal = null["marginals"][name]
        lines.append(
            f"| {name} | {obs['N145']['coordinates'][name]:.9g} | "
            f"{obs['N290']['coordinates'][name]:.9g} | "
            f"{null['residual'][name]:.9g} | {marginal['signed_z']:.4f} | "
            f"{marginal['p_value']:.6g} | {marginal['decision']} |"
        )
    lines.extend(
        [
            "",
            "## Frozen joint decision",
            "",
            f"The joint 2-vector score is chi-square={joint['chi_square']:.6g} "
            f"on {joint['degrees_of_freedom']} dof, p={joint['p_value']:.6g}: "
            f"**{joint['decision']}**.",
            "",
            "The parent and child streams are independent, so the residual covariance is "
            "the sum of their size-local nonlinear-jackknife covariance matrices.",
            "",
            "## Boundary",
            "",
            payload["interpretation_boundary"],
            "",
        ]
    )
    return "\n".join(lines)


def scientific_card(payload: Mapping[str, object]) -> str:
    null = payload["constant_response_null"]
    joint = null["joint"]
    marginals = null["marginals"]
    return "\n".join(
        [
            "# P118 scientific card",
            "",
            "- **Question:** Are the frozen metric-cancelling coordinates R_I and R_T constant from N=145 to N=290?",
            "- **Design:** Pre-target formula freeze; independent size streams; each nonlinear ratio rebuilt in every delete-one; joint 2-vector scored first.",
            f"- **Result:** joint chi-square={joint['chi_square']:.6g}/{joint['degrees_of_freedom']} dof, p={joint['p_value']:.6g}; R_I p={marginals['R_I']['p_value']:.6g}; R_T p={marginals['R_T']['p_value']:.6g}.",
            f"- **Meaning:** {joint['decision']}.",
            "- **Boundary:** This is a held-out reuse of P50 raw, not independent primary evidence and not a cross-model universality claim.",
            "",
        ]
    )


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "analysis/p118_metric_free_ratio_n145_n290_manifest.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--card", type=Path, required=True)
    args = parser.parse_args()
    payload = render(root, args.manifest)
    for path in (args.output, args.markdown, args.card):
        path.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.markdown.write_text(report_markdown(payload), encoding="utf-8")
    args.card.write_text(scientific_card(payload), encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
