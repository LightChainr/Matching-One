#!/usr/bin/env python3
"""Independently verify exact interval-comparison certificates for Issue 14."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping, Optional, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CERTIFICATE = ROOT / "analysis" / "synthetic_gadget_comparison_certificate.json"
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
FORBIDDEN_KEYS = frozenset(
    {
        "enumeration_results",
        "gadget_graph",
        "optimized_parameter",
        "production_data",
        "theorem_statement",
    }
)
Interval = Tuple[Fraction, Fraction]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _walk_json(value: Any, path: str = "$") -> None:
    _require(not isinstance(value, float), "%s contains floating point" % path)
    if isinstance(value, Mapping):
        bad = sorted(FORBIDDEN_KEYS.intersection(value))
        _require(not bad, "%s contains out-of-scope fields: %s" % (path, ",".join(bad)))
        for key, child in value.items():
            _walk_json(child, "%s.%s" % (path, key))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _walk_json(child, "%s[%d]" % (path, index))


def canonical_fraction(value: Any, label: str) -> Fraction:
    _require(type(value) is str, "%s must be an exact rational string" % label)
    try:
        parsed = Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError("%s is not a rational" % label) from error
    _require(str(parsed) == value, "%s is not canonically encoded" % label)
    return parsed


def probability_interval(value: Any, label: str) -> Tuple[str, Interval]:
    _require(isinstance(value, Mapping), "%s must be an interval object" % label)
    _require(set(value) == {"artifact_id", "lower", "upper"}, "%s interval fields drift" % label)
    artifact_id = value.get("artifact_id")
    _require(isinstance(artifact_id, str) and artifact_id, "%s artifact id missing" % label)
    lower = canonical_fraction(value.get("lower"), "%s.lower" % label)
    upper = canonical_fraction(value.get("upper"), "%s.upper" % label)
    _require(Fraction(0) <= lower <= upper <= Fraction(1), "%s is not a probability interval" % label)
    return artifact_id, (lower, upper)


def _validate_artifacts(value: Any) -> Mapping[str, Mapping[str, Any]]:
    _require(isinstance(value, list) and len(value) == 2, "exactly two synthetic artifacts required")
    artifacts = {}
    for index, artifact in enumerate(value):
        label = "artifact[%d]" % index
        _require(isinstance(artifact, Mapping), "%s must be an object" % label)
        _require(set(artifact) == {"id", "content", "sha256"}, "%s fields drift" % label)
        artifact_id = artifact.get("id")
        content = artifact.get("content")
        digest = artifact.get("sha256")
        _require(isinstance(artifact_id, str) and artifact_id.startswith("synthetic-"), "%s is not synthetic" % label)
        _require(artifact_id not in artifacts, "duplicate artifact id")
        _require(isinstance(content, str) and content.startswith("synthetic-"), "%s content is not synthetic" % label)
        _require(isinstance(digest, str) and SHA256_RE.fullmatch(digest) is not None, "%s SHA-256 format invalid" % label)
        actual_digest = hashlib.sha256(content.encode("utf-8")).hexdigest()
        _require(digest == actual_digest, "%s SHA-256 mismatch" % label)
        artifacts[artifact_id] = artifact
    _require(set(artifacts) == {"synthetic-baseline", "synthetic-candidate"}, "synthetic artifact roles drift")
    return artifacts


def _validate_definitions(value: Any) -> Mapping[str, Mapping[str, Any]]:
    _require(isinstance(value, list) and bool(value), "definitions must not be empty")
    definitions = {}
    for index, definition in enumerate(value):
        label = "definition[%d]" % index
        _require(isinstance(definition, Mapping), "%s must be an object" % label)
        _require(set(definition) == {"id", "parameter", "observable"}, "%s fields drift" % label)
        definition_id = definition.get("id")
        _require(isinstance(definition_id, str) and definition_id.startswith("synthetic-"), "%s is not synthetic" % label)
        _require(definition_id not in definitions, "duplicate definition id")
        observable = definition.get("observable")
        _require(isinstance(observable, str) and observable.startswith("SYNTHETIC ONLY:"), "%s observable is not synthetic" % label)
        parameter = definition.get("parameter")
        _require(isinstance(parameter, Mapping) and set(parameter) == {"id", "value"}, "%s parameter fields drift" % label)
        parameter_id = parameter.get("id")
        _require(isinstance(parameter_id, str) and parameter_id.startswith("synthetic-"), "%s parameter is not synthetic" % label)
        parameter_value = canonical_fraction(parameter.get("value"), "%s parameter value" % label)
        _require(Fraction(0) <= parameter_value <= Fraction(1), "%s parameter is outside [0,1]" % label)
        definitions[definition_id] = definition
    return definitions


def verify_certificate(certificate: Mapping[str, Any]) -> Mapping[str, Any]:
    _walk_json(certificate)
    _require(certificate.get("schema") == "matching-one/gadget-comparison-certificate/v1", "unknown schema")
    _require(certificate.get("issue") == 14, "wrong issue")
    _require(certificate.get("status") == "synthetic_verifier_fixture_no_bound", "certificate is not the synthetic fixture")
    arithmetic = certificate.get("arithmetic", {})
    _require(arithmetic.get("implementation") == "python_stdlib_fractions.Fraction", "arithmetic implementation drift")
    _require(
        arithmetic.get("endpoint_encoding") == "canonical_integer_or_numerator/denominator_string",
        "endpoint encoding drift",
    )
    _require(arithmetic.get("intervals") == "closed", "interval convention drift")
    _require(arithmetic.get("floating_point_allowed") is False, "floating point cannot be allowed")

    artifacts = _validate_artifacts(certificate.get("artifacts"))
    definitions = _validate_definitions(certificate.get("definitions"))
    claims = certificate.get("claims")
    _require(isinstance(claims, list) and bool(claims), "claims must not be empty")
    audited = []
    claim_ids = set()
    for index, claim in enumerate(claims):
        label = "claim[%d]" % index
        _require(isinstance(claim, Mapping), "%s must be an object" % label)
        _require(
            set(claim) == {"id", "definition_id", "parameter_id", "relation", "strict", "lhs", "rhs"},
            "%s fields drift" % label,
        )
        claim_id = claim.get("id")
        _require(isinstance(claim_id, str) and claim_id.startswith("synthetic-"), "%s id is not synthetic" % label)
        _require(claim_id not in claim_ids, "duplicate claim id")
        claim_ids.add(claim_id)
        definition_id = claim.get("definition_id")
        _require(definition_id in definitions, "%s references an unknown definition" % label)
        definition = definitions[definition_id]
        _require(claim.get("parameter_id") == definition["parameter"]["id"], "%s parameter definition mismatch" % label)
        lhs_artifact, lhs = probability_interval(claim.get("lhs"), "%s.lhs" % label)
        rhs_artifact, rhs = probability_interval(claim.get("rhs"), "%s.rhs" % label)
        _require(lhs_artifact in artifacts and rhs_artifact in artifacts, "%s references an unknown artifact" % label)
        _require(lhs_artifact != rhs_artifact, "%s must compare different artifacts" % label)
        relation = claim.get("relation")
        _require(relation in {"lhs_ge_rhs", "lhs_le_rhs"}, "%s relation unsupported" % label)
        strict = claim.get("strict")
        _require(type(strict) is bool, "%s strict flag must be boolean" % label)
        if relation == "lhs_ge_rhs":
            gap = lhs[0] - rhs[1]
        else:
            gap = rhs[0] - lhs[1]
        _require(gap > 0 if strict else gap >= 0, "%s interval comparison is not certified" % label)
        audited.append({
            "id": claim_id,
            "definition_id": definition_id,
            "relation": relation,
            "strict": strict,
            "separation": str(gap),
        })

    boundary = certificate.get("claim_boundary", {})
    _require(boundary.get("synthetic_fixture") is True, "synthetic boundary changed")
    _require(boundary.get("theorem_claim") is False, "fixture cannot claim a theorem")
    _require(boundary.get("new_bound_claim") is False, "fixture cannot claim a new bound")
    _require(boundary.get("parent_issue") == "remain open", "parent boundary changed")
    return {
        "schema": certificate["schema"],
        "status": "valid_synthetic_fixture_no_bound",
        "claims": audited,
        "floating_point_used": False,
        "proves_new_bound": False,
        "parent_issue": "remain open",
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    args = parser.parse_args(argv)
    certificate = json.loads(args.certificate.read_text(encoding="utf-8"))
    print(json.dumps(verify_certificate(certificate), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
