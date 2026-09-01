#!/usr/bin/env python3
"""Independently replay the arithmetic in the production E_top certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from certify_etop_production_models import (  # noqa: E402
    SCHEMA,
    canonical_bytes,
    certificate_from_projection,
    close,
)


def verify_certificate(
    certificate: dict, manifest: dict, *, verify_source: bool = False,
) -> dict[str, object]:
    if certificate.get("schema") != SCHEMA:
        raise ValueError("unexpected E_top certificate schema")
    if certificate.get("source") != manifest.get("source"):
        raise ValueError("certificate source does not match the frozen manifest")
    projection = certificate.get("source_projection")
    if not isinstance(projection, list) or not projection:
        raise ValueError("certificate source projection is absent")
    expected_projection_hash = hashlib.sha256(canonical_bytes(projection)).hexdigest()
    if certificate.get("projection_sha256") != expected_projection_hash:
        raise ValueError("embedded source projection hash changed")

    replay = certificate_from_projection(projection, manifest)
    for key in ("fixed_ratio_models", "free_ratio_model", "decision_summary"):
        if not close(certificate.get(key), replay[key], tolerance=1e-11):
            raise ValueError(f"E_top certificate arithmetic mismatch in {key}")

    source_verified = False
    if verify_source:
        source = manifest["source"]
        payload = subprocess.run(
            ["git", "show", f"{source['commit']}:{source['path']}"],
            cwd=ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
        if hashlib.sha256(payload).hexdigest() != source["sha256"]:
            raise ValueError("pinned production-derived source hash changed")
        source_verified = True

    return {
        "verified": True,
        "source_verified": source_verified,
        "projection_sha256": expected_projection_hash,
        "decisions": certificate["decision_summary"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificate", type=Path,
        default=ROOT / "results/etop-production-model-certificate/latest.json",
    )
    parser.add_argument(
        "--manifest", type=Path,
        default=ROOT / "analysis/etop_production_model_certificate_manifest.json",
    )
    parser.add_argument("--verify-source", action="store_true")
    args = parser.parse_args()
    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    print(json.dumps(
        verify_certificate(certificate, manifest, verify_source=args.verify_source),
        indent=2,
        sort_keys=True,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
