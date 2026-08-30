#!/usr/bin/env python3
"""Fail-closed validator for the algebraic-model certificate envelope."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "analysis/algebraic_model_certificate_manifest.yaml"
SCHEMA = "matching-one/algebraic-model-certificate-envelope/v1"
CLAIM_LEVELS = {"exact", "robust_statistical", "numerical_diagnostic", "inconclusive"}
TOP_LEVEL_FIELDS = {
    "schema",
    "issue",
    "certificate_id",
    "claim_level",
    "status",
    "inputs",
    "model_class",
    "gauge",
    "confidence_set",
    "certificate",
    "verification",
    "claim_boundary",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _safe_repository_path(raw: str) -> Path:
    path = Path(raw)
    _require(not path.is_absolute(), "input path must be repository-relative")
    resolved = (ROOT / path).resolve()
    _require(resolved == ROOT or ROOT in resolved.parents, "input path escapes repository")
    _require(resolved.is_file(), f"input path does not exist: {raw}")
    return resolved


def load_manifest(path: Path = DEFAULT_MANIFEST) -> Mapping[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    _require(isinstance(value, dict), "manifest must be a mapping")
    return value


def validate_manifest(manifest: Mapping[str, Any]) -> Mapping[str, Any]:
    fields = set(manifest)
    _require(fields == TOP_LEVEL_FIELDS, f"top-level fields drift: {sorted(fields ^ TOP_LEVEL_FIELDS)}")
    _require(manifest["schema"] == SCHEMA, "schema drift")
    _require(manifest["issue"] == 370, "issue drift")
    _require(isinstance(manifest["certificate_id"], str) and manifest["certificate_id"], "missing certificate id")
    _require(manifest["claim_level"] in CLAIM_LEVELS, "invalid claim level")
    _require(manifest["status"] == "verified", "unverified envelope")

    inputs = manifest["inputs"]
    _require(isinstance(inputs, list) and inputs, "at least one immutable input is required")
    seen_paths = set()
    for item in inputs:
        _require(set(item) == {"path", "sha256", "descriptor", "dependency_group", "chronology"}, "input descriptor fields drift")
        _require(item["path"] not in seen_paths, "duplicate input path")
        seen_paths.add(item["path"])
        path = _safe_repository_path(item["path"])
        _require(len(item["sha256"]) == 64 and item["sha256"] == _sha256_file(path), "input digest mismatch")
        for field in ("descriptor", "dependency_group", "chronology"):
            _require(isinstance(item[field], str) and item[field], f"empty input {field}")

    model = manifest["model_class"]
    _require(set(model) == {"id", "state_dimension", "relations", "parameter_bounds"}, "model fields drift")
    _require(isinstance(model["state_dimension"], int) and model["state_dimension"] >= 1, "invalid state dimension")
    _require(isinstance(model["relations"], list) and model["relations"], "model relations required")
    _require(isinstance(model["parameter_bounds"], list), "parameter bounds must be explicit")

    gauge = manifest["gauge"]
    _require(set(gauge) == {"kind", "coverage", "uncovered_set"}, "gauge fields drift")
    _require(gauge["kind"] in {"invariant_coordinates", "reachable_source", "observable_readout", "hankel_minor"}, "invalid gauge kind")
    _require(gauge["coverage"] in {"complete_for_declared_class", "chart_only"}, "invalid gauge coverage")
    if gauge["coverage"] == "chart_only":
        _require(isinstance(gauge["uncovered_set"], str) and gauge["uncovered_set"], "chart-only gauge must declare uncovered set")
    else:
        _require(gauge["uncovered_set"] == "none", "complete gauge must have no uncovered set")

    confidence = manifest["confidence_set"]
    _require(set(confidence) == {"kind", "covariance", "threshold", "calibration"}, "confidence fields drift")
    _require(confidence["kind"] in {"exact_point", "outer_confidence_set"}, "invalid confidence kind")
    if confidence["kind"] == "exact_point":
        _require(confidence["covariance"] == "not_applicable", "exact point cannot declare covariance")
        _require(confidence["threshold"] == "zero", "exact point threshold must be zero")

    certificate = manifest["certificate"]
    _require(set(certificate) == {"level", "type", "order", "result", "dual_support"}, "certificate fields drift")
    _require(certificate["level"] in {"E", "S"}, "invalid certificate level")
    _require(isinstance(certificate["order"], int) and certificate["order"] >= 0, "invalid certificate order")

    verification = manifest["verification"]
    _require(set(verification) == {"verifier", "arithmetic", "status", "residual"}, "verification fields drift")
    _require(_safe_repository_path(verification["verifier"]).is_file(), "verifier does not exist")
    _require(verification["arithmetic"] in {"rational", "integer", "interval"}, "invalid arithmetic")
    _require(verification["status"] == "verified_exact", "verification did not pass exactly")
    _require(verification["residual"] == "0", "nonzero verification residual")

    boundary = manifest["claim_boundary"]
    _require(set(boundary) == {"included", "excluded", "parent_issue"}, "claim boundary fields drift")
    _require(boundary["parent_issue"] == "remain open", "parent issue boundary drift")
    _require(boundary["included"] and boundary["excluded"], "claim boundary must be explicit")
    return {
        "schema": manifest["schema"],
        "certificate_id": manifest["certificate_id"],
        "claim_level": manifest["claim_level"],
        "input_count": len(inputs),
        "model_class": model["id"],
        "gauge_coverage": gauge["coverage"],
        "status": "valid_fail_closed_envelope",
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", nargs="?", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args(argv)
    print(json.dumps(validate_manifest(load_manifest(args.manifest)), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
