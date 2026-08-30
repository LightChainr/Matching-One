#!/usr/bin/env python3
"""Exact-rational confidence certificate over the production E_top rank plane."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "analysis" / "etop_production_rational_confidence_manifest.json"
DEFAULT_OUTPUT = ROOT / "results" / "model-certificates" / "production" / "etop-rational-confidence" / "latest.json"
SCHEMA = "matching-one/etop-production-rational-confidence-certificate/v1"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def git_bytes(commit: str, path: str) -> bytes:
    return subprocess.run(
        ["git", "show", f"{commit}:{path}"], cwd=ROOT, check=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    ).stdout


def git_blob(commit: str, path: str) -> dict[str, Any]:
    blob = subprocess.run(
        ["git", "rev-parse", f"{commit}:{path}"], cwd=ROOT, check=True,
        capture_output=True, text=True,
    ).stdout.strip()
    size = int(subprocess.run(
        ["git", "cat-file", "-s", blob], cwd=ROOT, check=True,
        capture_output=True, text=True,
    ).stdout.strip())
    return {"path": path, "git_blob_sha": blob, "bytes": size}


def fraction(value: object) -> Fraction:
    return Fraction(str(value))


def ftext(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def decimal(value: Fraction, digits: int = 17) -> float:
    return float(format(float(value), f".{digits}g"))


def selected_rows(source: Mapping[str, Any], manifest: Mapping[str, Any]) -> list[dict[str, Any]]:
    sizes = list(manifest["selection"]["expected_sizes"])
    ids = list(manifest["selection"]["expected_ids"])
    require(list(source["high_statistics_summary"]["sizes"]) == sizes,
            "source high-statistics selection changed")
    by_size = {int(row["N"]): row for row in source["datasets"]}
    rows = []
    for expected_id, size in zip(ids, sizes):
        row = by_size[size]
        require(row["id"] == expected_id, f"dataset identity changed at N={size}")
        order = list(row["covariance_metric_order"])
        ia, ie = order.index("P4_A_top"), order.index("P4_E_top")
        covariance = row["covariance_intrinsic_center_first_order_influence"]
        identity_residuals = row["exact_basis_identity_residuals"]
        require(max(abs(float(value)) for value in identity_residuals.values()) <= 1e-30,
                f"exact rank-plane identity failed at {expected_id}")
        rows.append({
            "id": expected_id,
            "N": size,
            "batches": int(row["batches"]),
            "samples_per_orientation": int(row["samples_per_orientation"]),
            "orientations": row["orientations"],
            "A": fraction(row["point"]["P4_A_top"]),
            "E": fraction(row["point"]["P4_E_top"]),
            "var_A": fraction(covariance[ia][ia]),
            "cov_AE": fraction(covariance[ia][ie]),
            "var_E": fraction(covariance[ie][ie]),
            "histogram": str(row["histogram"]),
            "moments": str(row["moments"]),
        })
    return rows


def confidence_row(row: Mapping[str, Any], ratio: Fraction, critical: Fraction) -> dict[str, Any]:
    residual = row["E"] - ratio * row["A"]
    variance = row["var_E"] + ratio * ratio * row["var_A"] - 2 * ratio * row["cov_AE"]
    require(variance > 0, f"nonpositive variance at {row['id']}")
    margin = residual * residual - critical * critical * variance
    return {
        "id": row["id"],
        "N": row["N"],
        "residual": ftext(residual),
        "variance": ftext(variance),
        "squared_feasibility_margin": ftext(margin),
        "margin_decimal": decimal(margin),
        "inside_outer_band": margin <= 0,
        "exact_test": "(E-r*A)^2 <= K^2*(var_E+r^2*var_A-2r*cov_AE)",
    }


def fixed_model(rows: Sequence[Mapping[str, Any]], model: Mapping[str, Any],
                critical: Fraction) -> dict[str, Any]:
    ratio = Fraction(model["r"])
    scored = [confidence_row(row, ratio, critical) for row in rows]
    incompatible = [row["id"] for row in scored if not row["inside_outer_band"]]
    return {
        "id": model["id"],
        "relation": model["relation"],
        "ratio": ftext(ratio),
        "rows": scored,
        "incompatible_rows": incompatible,
        "decision": "eliminated" if incompatible else "not_eliminated",
        "certificate_logic": "one incompatible row excludes a universal fixed-ratio model",
    }


def serialize_projection(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    return [{
        "id": row["id"], "N": row["N"], "batches": row["batches"],
        "samples_per_orientation": row["samples_per_orientation"],
        "orientations": row["orientations"],
        "A_top": ftext(row["A"]), "E_top": ftext(row["E"]),
        "covariance_AE": [
            [ftext(row["var_A"]), ftext(row["cov_AE"])],
            [ftext(row["cov_AE"]), ftext(row["var_E"])],
        ],
        "histogram": row["histogram"], "moments": row["moments"],
    } for row in rows]


def build(manifest_path: Path = DEFAULT_MANIFEST) -> dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    require(manifest.get("schema") == "matching-one/etop-production-rational-confidence-manifest/v1",
            "manifest schema mismatch")
    source_spec = manifest["source"]
    payload = git_bytes(source_spec["commit"], source_spec["path"])
    require(sha256_bytes(payload) == source_spec["sha256"], "pinned source hash mismatch")
    source = json.loads(payload)
    rows = selected_rows(source, manifest)
    critical = Fraction(manifest["confidence"]["critical_value"])
    fixed = [fixed_model(rows, model, critical) for model in manifest["fixed_ratio_models"]]
    free_ratio = Fraction(manifest["free_ratio_model"]["feasible_witness_r"])
    free_rows = [confidence_row(row, free_ratio, critical) for row in rows]
    require(all(row["inside_outer_band"] for row in free_rows),
            "frozen common-line witness left the confidence set")
    projection = serialize_projection(rows)
    raw_paths = sorted({row[key] for row in rows for key in ("histogram", "moments")})
    nominal_critical = statistics.NormalDist().inv_cdf(1 - 0.01 / 16)
    require(float(critical) > nominal_critical, "rational critical is not conservative")
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "robust_statistical",
        "status": "verified_rational_outer_confidence_model_ordering",
        "source": {**source_spec, "git_blob": git_blob(source_spec["commit"], source_spec["path"])["git_blob_sha"]},
        "source_projection": projection,
        "projection_sha256": sha256_bytes(canonical_bytes(projection)),
        "raw_production_git_objects": [git_blob(source_spec["commit"], path) for path in raw_paths],
        "confidence_set": {
            **manifest["confidence"],
            "rational_critical_decimal": float(critical),
            "nominal_critical_check": "runtime NormalDist z_(1-0.01/16) is strictly below 13/4",
            "exact_arithmetic_boundary": (
                "all feasibility margins are exact rational consequences of the pinned decimal "
                "estimates/covariances and K=13/4; Gaussian calibration itself is an assumed outer-set contract"
            ),
        },
        "fixed_ratio_models": fixed,
        "free_ratio_model": {
            "id": manifest["free_ratio_model"]["id"],
            "relation": manifest["free_ratio_model"]["relation"],
            "feasible_witness_r": ftext(free_ratio),
            "rows": free_rows,
            "decision": "not_eliminated",
            "certificate_logic": "one common rational r lies inside all eight simultaneous line-distance bands",
        },
        "decision_summary": {
            "eliminated": [model["id"] for model in fixed if model["decision"] == "eliminated"],
            "not_eliminated": [manifest["free_ratio_model"]["id"]],
            "model_ordering": "free_common_line strictly outranks all three frozen fixed-ratio lines",
        },
        "claim_boundary": {
            "included": "exact rational verification of model feasibility margins inside the declared production Gaussian outer confidence contract",
            "excluded": "finite-sample non-Gaussian coverage, cross-archive independence pooling, exact physical dimension, continuum field identity, exponent or SOS certificate",
            "parent_issue": "remain open",
        },
    }


def validate(value: Mapping[str, Any], manifest_path: Path = DEFAULT_MANIFEST) -> dict[str, Any]:
    expected = build(manifest_path)
    require(value == expected, "production confidence certificate does not exactly reproduce")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_rational_production_confidence_certificate",
        "production_rows": len(value["source_projection"]),
        "eliminated": value["decision_summary"]["eliminated"],
        "surviving_witness_r": value["free_ratio_model"]["feasible_witness_r"],
        "source_sha256": value["source"]["sha256"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        value = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate(value, args.manifest), indent=2, sort_keys=True))
        return 0
    result = build(args.manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
