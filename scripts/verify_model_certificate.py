#!/usr/bin/env python3
"""Fail-closed dispatcher for checked-in exact model-certificate types."""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/dispatcher/latest.json"
SCHEMA = "matching-one/model-certificate-dispatcher/v1"
REGISTRY = {
    "matching-one/exact-hankel-minor-certificate/v1": ("exact_hankel_minor_certificate", "results/model-certificates/framework/hankel-minor/latest.json"),
    "matching-one/exact-linear-ideal-certificate/v1": ("exact_linear_ideal_certificate", "results/model-certificates/framework/linear-ideal/latest.json"),
    "matching-one/exact-rational-realization-certificate/v1": ("exact_rational_realization_certificate", "results/model-certificates/framework/rational-realization/latest.json"),
    "matching-one/exact-finite-recurrence-certificate/v1": ("exact_recurrence_certificate", "results/model-certificates/framework/finite-recurrence/latest.json"),
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _safe_path(raw: str) -> Path:
    path = Path(raw)
    _require(not path.is_absolute(), "certificate path must be repository-relative")
    resolved = (ROOT / path).resolve()
    _require(ROOT in resolved.parents, "certificate path escapes repository")
    _require(resolved.is_file(), f"certificate path does not exist: {raw}")
    return resolved


def verify_bundle(entries: Sequence[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    _require(isinstance(entries, list) and entries, "certificate bundle cannot be empty")
    seen_paths = set()
    summaries = []
    for entry in entries:
        _require(set(entry) == {"path", "sha256", "schema"}, "bundle entry fields drift")
        raw_path = entry["path"]
        _require(raw_path not in seen_paths, "duplicate certificate path")
        seen_paths.add(raw_path)
        path = _safe_path(raw_path)
        _require(entry["schema"] in REGISTRY, "unregistered certificate schema")
        module_name, canonical_path = REGISTRY[entry["schema"]]
        _require(raw_path == canonical_path, "certificate path is not canonical for schema")
        _require(entry["sha256"] == _sha256_file(path), "certificate digest mismatch")
        value = json.loads(path.read_text(encoding="utf-8"))
        _require(value.get("schema") == entry["schema"], "certificate schema/content mismatch")
        module = importlib.import_module(module_name)
        summary = module.validate_result(value)
        summaries.append({"path": raw_path, "schema": entry["schema"], "verifier": module_name, "summary": summary})
    return summaries


def frozen_entries() -> list[dict[str, str]]:
    entries = []
    for schema, (_, raw_path) in sorted(REGISTRY.items()):
        path = _safe_path(raw_path)
        entries.append({"path": raw_path, "sha256": _sha256_file(path), "schema": schema})
    return entries


def build_result() -> dict[str, Any]:
    entries = frozen_entries()
    summaries = verify_bundle(entries)
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact_verification_dispatch",
        "entries": entries,
        "verification": {
            "certificate_count": len(summaries),
            "all_registered": True,
            "all_digests_match": True,
            "all_type_specific_verifiers_pass": True,
            "summaries": summaries,
            "status": "exact_certificate_bundle_verified",
        },
        "claim_boundary": {
            "included": "fail-closed path, digest, schema whitelist, canonical-location, and type-specific verification of four checked-in exact framework certificates",
            "excluded": "arbitrary plugin loading, unregistered schemas, certificate discovery, SOS solving, statistical calibration, or validation of any physical model class",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    _require(result == expected, "dispatcher artifact does not exactly reproduce")
    return {
        "schema": result["schema"],
        "status": "valid_fail_closed_model_certificate_dispatcher",
        "certificate_count": result["verification"]["certificate_count"],
        "schemas": [entry["schema"] for entry in result["entries"]],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        value = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_result(value), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_result(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
