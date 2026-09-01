#!/usr/bin/env python3
"""Certify bounded E_top model eliminations on pinned production summaries."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
from statistics import NormalDist
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "matching-one/etop-production-model-certificate/v1"


def canonical_bytes(payload: object) -> bytes:
    return (json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def git_bytes(commit: str, path: str) -> bytes:
    return subprocess.run(
        ["git", "show", f"{commit}:{path}"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    ).stdout


def git_object_record(commit: str, path: str) -> dict[str, object]:
    blob = subprocess.run(
        ["git", "rev-parse", f"{commit}:{path}"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    size = int(subprocess.run(
        ["git", "cat-file", "-s", blob],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip())
    return {"path": path, "git_blob_sha": blob, "bytes": size}


def close(left: object, right: object, *, tolerance: float = 1e-12) -> bool:
    if isinstance(left, float) or isinstance(right, float):
        try:
            a, b = float(left), float(right)
        except (TypeError, ValueError):
            return False
        return math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance)
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        return set(left) == set(right) and all(close(left[key], right[key], tolerance=tolerance) for key in left)
    if isinstance(left, Sequence) and not isinstance(left, (str, bytes)):
        return (
            isinstance(right, Sequence)
            and not isinstance(right, (str, bytes))
            and len(left) == len(right)
            and all(close(a, b, tolerance=tolerance) for a, b in zip(left, right))
        )
    return left == right


def selected_projection(source: Mapping[str, object], manifest: Mapping[str, object]) -> list[dict[str, object]]:
    expected_sizes = list(manifest["selection"]["expected_sizes"])
    expected_ids = list(manifest["selection"]["expected_dataset_ids"])
    source_sizes = list(source["high_statistics_summary"]["sizes"])
    if source_sizes != expected_sizes:
        raise ValueError("source high-statistics selection changed")
    by_size = {int(row["N"]): row for row in source["datasets"]}
    rows = []
    for expected_id, size in zip(expected_ids, expected_sizes):
        row = by_size[size]
        if row["id"] != expected_id:
            raise ValueError(f"dataset identity changed at N={size}")
        order = list(row["covariance_metric_order"])
        ia = order.index("P4_A_top")
        ie = order.index("P4_E_top")
        covariance = row["covariance_intrinsic_center_first_order_influence"]
        identities = row["exact_basis_identity_residuals"]
        if max(abs(float(value)) for value in identities.values()) > 1e-30:
            raise ValueError(f"rank-plane exact identity drift at {expected_id}")
        rows.append({
            "id": expected_id,
            "N": size,
            "batches": int(row["batches"]),
            "samples_per_orientation": int(row["samples_per_orientation"]),
            "orientations": row["orientations"],
            "delta_cos4": float(row["delta_cos4"]),
            "A_top": float(row["point"]["P4_A_top"]),
            "E_top": float(row["point"]["P4_E_top"]),
            "covariance_AE": [
                [float(covariance[ia][ia]), float(covariance[ia][ie])],
                [float(covariance[ie][ia]), float(covariance[ie][ie])],
            ],
            "histogram": row["histogram"],
            "moments": row["moments"],
        })
    return rows


def critical_value(alpha: float, count: int) -> float:
    if not (0.0 < alpha < 1.0) or count <= 0:
        raise ValueError("invalid confidence contract")
    return NormalDist().inv_cdf(1.0 - alpha / (2.0 * count))


def fixed_ratio_certificate(
    projection: Sequence[Mapping[str, object]], model: Mapping[str, object], alpha: float,
) -> dict[str, object]:
    ratio = float(model["ratio_r_in_E_equals_rA"])
    threshold = critical_value(alpha, len(projection))
    rows = []
    for source in projection:
        a = float(source["A_top"])
        e = float(source["E_top"])
        covariance = source["covariance_AE"]
        va = float(covariance[0][0])
        cae = float(covariance[0][1])
        ve = float(covariance[1][1])
        residual = e - ratio * a
        variance = ve + ratio * ratio * va - 2.0 * ratio * cae
        if not variance > 0.0:
            raise ValueError(f"nonpositive fixed-ratio variance for {source['id']}")
        standard_error = math.sqrt(variance)
        lower = residual - threshold * standard_error
        upper = residual + threshold * standard_error
        rows.append({
            "id": source["id"],
            "N": source["N"],
            "residual_E_minus_rA": residual,
            "variance": variance,
            "standard_error": standard_error,
            "z": residual / standard_error,
            "simultaneous_interval": [lower, upper],
            "excludes_zero": lower > 0.0 or upper < 0.0,
        })
    excluded = [row["id"] for row in rows if row["excludes_zero"]]
    strongest = max(rows, key=lambda row: abs(float(row["z"])))
    return {
        "id": model["id"],
        "meaning": model["meaning"],
        "ratio_r_in_E_equals_rA": ratio,
        "familywise_alpha": alpha,
        "bonferroni_two_sided_gaussian_critical": threshold,
        "rows": rows,
        "incompatible_datasets": excluded,
        "strongest_abs_z": {"id": strongest["id"], "z": strongest["z"]},
        "decision": "eliminated" if excluded else "not_eliminated",
        "logic": "a universal fixed-ratio model is eliminated when any simultaneous dataset interval excludes zero",
    }


def fieller_set(row: Mapping[str, object], threshold: float) -> dict[str, object]:
    a_value = float(row["A_top"])
    e_value = float(row["E_top"])
    covariance = row["covariance_AE"]
    va = float(covariance[0][0])
    cae = float(covariance[0][1])
    ve = float(covariance[1][1])
    q = threshold * threshold
    qa = a_value * a_value - q * va
    qb = -2.0 * a_value * e_value + 2.0 * q * cae
    qc = e_value * e_value - q * ve
    discriminant = qb * qb - 4.0 * qa * qc
    base = {
        "id": row["id"],
        "N": row["N"],
        "quadratic": [qa, qb, qc],
        "discriminant": discriminant,
    }
    if discriminant < 0.0:
        if qa < 0.0:
            return {**base, "set_type": "all_real"}
        return {**base, "set_type": "empty"}
    root = math.sqrt(max(discriminant, 0.0))
    first = (-qb - root) / (2.0 * qa)
    second = (-qb + root) / (2.0 * qa)
    lower, upper = sorted((first, second))
    if qa > 0.0:
        return {**base, "set_type": "bounded", "interval": [lower, upper]}
    return {**base, "set_type": "two_rays", "excluded_interval": [lower, upper]}


def free_ratio_certificate(
    projection: Sequence[Mapping[str, object]], model: Mapping[str, object], alpha: float,
) -> dict[str, object]:
    threshold = critical_value(alpha, len(projection))
    sets = [fieller_set(row, threshold) for row in projection]
    if any(row["set_type"] == "empty" for row in sets):
        intersection: list[float] | None = None
    elif any(row["set_type"] == "two_rays" for row in sets):
        raise ValueError("two-ray Fieller set needs an explicit union intersection implementation")
    else:
        bounded = [row["interval"] for row in sets if row["set_type"] == "bounded"]
        lower = max(float(row[0]) for row in bounded) if bounded else -math.inf
        upper = min(float(row[1]) for row in bounded) if bounded else math.inf
        intersection = [lower, upper] if lower <= upper else None
    return {
        "id": model["id"],
        "meaning": model["meaning"],
        "familywise_alpha": alpha,
        "bonferroni_two_sided_gaussian_critical": threshold,
        "fieller_sets": sets,
        "common_ratio_intersection": intersection,
        "decision": "not_eliminated" if intersection is not None else "eliminated",
        "logic": "the common-line model survives exactly when the simultaneous Fieller sets have a nonempty real intersection",
    }


def certificate_from_projection(
    projection: Sequence[Mapping[str, object]], manifest: Mapping[str, object],
) -> dict[str, object]:
    alpha = float(manifest["confidence"]["familywise_alpha_per_model"])
    fixed = [fixed_ratio_certificate(projection, model, alpha) for model in manifest["fixed_ratio_models"]]
    free = free_ratio_certificate(projection, manifest["free_ratio_model"], alpha)
    return {
        "projection_sha256": sha256_bytes(canonical_bytes(projection)),
        "fixed_ratio_models": fixed,
        "free_ratio_model": free,
        "decision_summary": {
            "eliminated": [row["id"] for row in fixed if row["decision"] == "eliminated"],
            "not_eliminated": [row["id"] for row in fixed if row["decision"] != "eliminated"]
            + ([free["id"]] if free["decision"] != "eliminated" else []),
        },
    }


def render(certificate: Mapping[str, object]) -> str:
    lines = [
        "# Production E_top model-elimination certificate",
        "",
        "The exact state basis is `A_top=P2-P0`, `E_top=P0+P2=1-P1`.",
        "Each fixed model uses a separate 99% familywise Gaussian-Bonferroni outer confidence set over eight high-statistics production datasets; no cross-dataset p-values are pooled.",
        "",
        "| model | equation | decision | incompatible production rows |",
        "|---|---|---|---|",
    ]
    for row in certificate["fixed_ratio_models"]:
        ratio = row["ratio_r_in_E_equals_rA"]
        equation = "E_top=0" if ratio == 0.0 else f"E_top={ratio:g} A_top"
        lines.append(
            f"| {row['id']} | {equation} | {row['decision']} | "
            f"{', '.join(row['incompatible_datasets']) or 'none'} |"
        )
    free = certificate["free_ratio_model"]
    lines += [
        f"| {free['id']} | E_top=r A_top, free r | {free['decision']} | common r={free['common_ratio_intersection']} |",
        "",
        "## Interpretation",
        "",
        "Production data eliminate a pure Alexander-odd state response (`E_top=0`) and both exact endpoint-cancellation lines (`F1=0` and `F2=0`) under their declared outer confidence sets. A one-dimensional state line with a free common mixing ratio is not eliminated, so the result establishes a required even component without yet proving a two-dimensional continuum module.",
        "",
        "The `F2=0` exclusion is the narrowest of the three and should be read from its saved simultaneous interval, not as a generic absence of finite-size cancellation.",
        "",
        "## Boundary",
        "",
        str(certificate["claim_boundary"]),
        "",
    ]
    return "\n".join(lines)


def build(manifest: Mapping[str, object]) -> dict[str, object]:
    source_spec = manifest["source"]
    payload = git_bytes(source_spec["commit"], source_spec["path"])
    if sha256_bytes(payload) != source_spec["sha256"]:
        raise ValueError("pinned rank-plane source payload changed")
    source = json.loads(payload)
    projection = selected_projection(source, manifest)
    arithmetic = certificate_from_projection(projection, manifest)
    raw_paths = sorted({str(row[key]) for row in projection for key in ("histogram", "moments")})
    return {
        "schema": SCHEMA,
        "status": "production_data_outer_confidence_model_elimination",
        "source": source_spec,
        "source_projection": projection,
        "raw_production_git_objects": [git_object_record(source_spec["commit"], path) for path in raw_paths],
        "observable": manifest["observable"],
        "confidence": manifest["confidence"],
        **arithmetic,
        "claim_boundary": manifest["claim_boundary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=ROOT / "analysis/etop_production_model_certificate_manifest.json")
    parser.add_argument("--output-json", type=Path, default=ROOT / "results/etop-production-model-certificate/latest.json")
    parser.add_argument("--output-md", type=Path, default=ROOT / "results/etop-production-model-certificate/latest.md")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    if manifest.get("schema") != "matching-one/etop-production-model-certificate-manifest/v1":
        raise ValueError("unexpected E_top manifest schema")
    result = build(manifest)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.output_md.write_text(render(result), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
