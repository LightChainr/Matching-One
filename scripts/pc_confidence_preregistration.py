#!/usr/bin/env python3
"""Content-addressed pre-registration and trial-ledger validator.

The validator freezes the statistical constants already derived for Issue 112.
It audits metadata and exact binomial decisions; it does not generate random
fields or certify that a declared source supplied independent trials.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Sequence

import rigorous_pc_confidence_gate as gate


HEX_DIGITS = frozenset("0123456789abcdef")
SIDES = ("upper", "lower")
EXPECTED_GRAPH = {"upper": "square", "lower": "matching"}


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def plan_digest(plan: Mapping[str, Any]) -> str:
    return sha256_text(canonical_json(plan))


def is_sha256(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and set(value) <= HEX_DIGITS


def build_plan() -> dict[str, Any]:
    cutoff, tail = gate.minimal_successes(400, gate.MODERN_P0, gate.PER_RUN_ALPHA)
    assert cutoff == 373
    return {
        "schema": "matching-one/pc-confidence-preregistration-plan/v1",
        "issue": 112,
        "status": "protocol_only_no_trial_data",
        "statistical_contract": {
            "familywise_alpha": gate.fraction_text(gate.FAMILYWISE_ALPHA),
            "sides": list(SIDES),
            "max_attempts_per_side": gate.ATTEMPTS_PER_SIDE,
            "trials_per_attempt": 400,
            "null_block_probability": gate.fraction_text(gate.MODERN_P0),
            "per_run_alpha": gate.fraction_text(gate.PER_RUN_ALPHA),
            "minimum_successes": cutoff,
            "tail_at_cutoff_decimal": gate.decimal_text(tail),
            "decision_rule": "accept one run iff exact Bin(N,p0) upper tail <= per_run_alpha",
        },
        "side_contract": {
            "upper": {
                "graph": "square",
                "conclusion_shape": "pc(square)<=tested_parameter",
            },
            "lower": {
                "graph": "matching",
                "conclusion_shape": "pc(square)>=1-tested_parameter",
            },
        },
        "record_contract": {
            "phase": "final",
            "required_fields": [
                "record_id",
                "plan_digest",
                "phase",
                "side",
                "attempt",
                "graph",
                "tested_parameter",
                "trials",
                "successes",
                "stream_domain",
                "data_digest",
                "independence_attestation",
            ],
            "unique_fields": ["record_id", "stream_domain", "data_digest"],
            "unique_pair": ["side", "attempt"],
            "independence_attestation": (
                "declarative provenance only; the validator cannot prove genuine independence"
            ),
        },
        "forbidden_exploration_data_digests": [],
        "stopping_rule": (
            "at most three pre-indexed attempts per side; stopping after an accepted attempt is allowed "
            "under the frozen Bonferroni allocation"
        ),
    }


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def audit_records(plan: Mapping[str, Any], records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    contract = plan["statistical_contract"]
    digest = plan_digest(plan)
    required = set(plan["record_contract"]["required_fields"])
    trials = int(contract["trials_per_attempt"])
    cutoff = int(contract["minimum_successes"])
    max_attempts = int(contract["max_attempts_per_side"])
    null_probability = gate.MODERN_P0
    per_run_alpha = gate.PER_RUN_ALPHA
    forbidden = set(plan["forbidden_exploration_data_digests"])
    seen_ids = set()
    seen_domains = set()
    seen_data = set()
    seen_attempts = set()
    decisions = []

    for position, record in enumerate(records):
        missing = sorted(required - set(record))
        _require(not missing, "record %d missing fields: %s" % (position, ",".join(missing)))
        _require(record["plan_digest"] == digest, "record %d plan digest mismatch" % position)
        _require(record["phase"] == "final", "record %d is not a final trial" % position)
        side = record["side"]
        _require(side in SIDES, "record %d has unknown side" % position)
        _require(record["graph"] == EXPECTED_GRAPH[side], "record %d graph/side mismatch" % position)
        attempt = record["attempt"]
        _require(isinstance(attempt, int) and 1 <= attempt <= max_attempts, "invalid attempt index")
        key = (side, attempt)
        _require(key not in seen_attempts, "duplicate side/attempt")
        seen_attempts.add(key)
        _require(record["trials"] == trials, "record %d trial count differs from plan" % position)
        successes = record["successes"]
        _require(isinstance(successes, int) and 0 <= successes <= trials, "invalid success count")
        _require(
            isinstance(record["tested_parameter"], str) and record["tested_parameter"],
            "tested_parameter must be a nonempty exact string",
        )
        _require(
            isinstance(record["independence_attestation"], str)
            and bool(record["independence_attestation"].strip()),
            "independence attestation is required",
        )

        record_id = record["record_id"]
        stream_domain = record["stream_domain"]
        data_digest = record["data_digest"]
        _require(isinstance(record_id, str) and record_id, "record_id must be nonempty")
        _require(isinstance(stream_domain, str) and stream_domain, "stream_domain must be nonempty")
        _require(is_sha256(data_digest), "data_digest must be lowercase SHA-256")
        _require(record_id not in seen_ids, "duplicate record_id")
        _require(stream_domain not in seen_domains, "duplicate stream_domain")
        _require(data_digest not in seen_data, "duplicate data_digest")
        _require(data_digest not in forbidden, "final data reuses an exploration digest")
        seen_ids.add(record_id)
        seen_domains.add(stream_domain)
        seen_data.add(data_digest)

        tail = gate.binomial_tail(trials, successes, null_probability)
        accepted = tail <= per_run_alpha
        _require(accepted == (successes >= cutoff), "cutoff and exact-tail decisions disagree")
        decisions.append(
            {
                "record_id": record_id,
                "side": side,
                "attempt": attempt,
                "successes": successes,
                "accepted": accepted,
                "null_tail_decimal": gate.decimal_text(tail),
                "tail_arithmetic": "exact Fraction",
            }
        )

    return {
        "plan_digest": digest,
        "record_count": len(records),
        "unique_side_attempts": len(seen_attempts),
        "accepted_records": sum(1 for row in decisions if row["accepted"]),
        "decisions": decisions,
        "audit_boundary": (
            "metadata and arithmetic validated; the independence attestation remains an external claim"
        ),
    }


def synthetic_fixture(plan: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Return visibly synthetic records for validator regression tests only."""

    digest = plan_digest(plan)
    successes = {"upper": (372, 373, 371), "lower": (373, 370, 369)}
    records = []
    for side in SIDES:
        for attempt, count in enumerate(successes[side], start=1):
            label = "synthetic-%s-%d" % (side, attempt)
            records.append(
                {
                    "record_id": label,
                    "plan_digest": digest,
                    "phase": "final",
                    "side": side,
                    "attempt": attempt,
                    "graph": EXPECTED_GRAPH[side],
                    "tested_parameter": "synthetic-only",
                    "trials": 400,
                    "successes": count,
                    "stream_domain": "matching-one/issue112/%s/%d" % (side, attempt),
                    "data_digest": sha256_text(label + "/not-empirical-data"),
                    "independence_attestation": "synthetic validator fixture; not an empirical claim",
                }
            )
    return records


def build_artifact() -> dict[str, Any]:
    plan = build_plan()
    records = synthetic_fixture(plan)
    audit = audit_records(plan, records)
    assert audit["record_count"] == 6
    assert audit["accepted_records"] == 2
    assert plan["statistical_contract"]["minimum_successes"] == 373
    return {
        "schema": "matching-one/pc-confidence-preregistration-audit/v1",
        "issue": 112,
        "status": "protocol_validator_no_empirical_data",
        "plan": plan,
        "plan_sha256": plan_digest(plan),
        "synthetic_fixture_audit": audit,
        "claim_boundary": {
            "proved": "plan integrity, ledger uniqueness, count ranges, and exact decision arithmetic",
            "not_proved": (
                "genuine randomness, independent Bernoulli trials, an observed event probability, "
                "or a critical-probability bound"
            ),
        },
    }


def render_markdown(artifact: Mapping[str, Any]) -> str:
    contract = artifact["plan"]["statistical_contract"]
    audit = artifact["synthetic_fixture_audit"]
    lines = [
        "# Pre-registration and trial-ledger audit gate",
        "",
        "This artifact freezes a content-addressed statistical plan. It contains only synthetic",
        "validator fixtures and makes no empirical percolation claim.",
        "",
        "## Frozen plan",
        "",
        "- plan SHA-256: `%s`;" % artifact["plan_sha256"],
        "- familywise alpha: `%s`;" % contract["familywise_alpha"],
        "- two sides, at most `%d` attempts per side;" % contract["max_attempts_per_side"],
        "- `%d` trials per attempt; per-run alpha `%s`;" % (
            contract["trials_per_attempt"], contract["per_run_alpha"]
        ),
        "- exact acceptance cutoff: `%d` successes." % contract["minimum_successes"],
        "",
        "## Validator controls",
        "",
        "The synthetic six-record fixture has `%d` accepted records and exercises both sides. Tests"
        % audit["accepted_records"],
        "also require rejection of plan-digest tampering, duplicate side/attempt pairs, duplicate",
        "domains or data digests, exploration-data reuse, graph/side mismatches, and invalid counts.",
        "",
        "## Boundary",
        "",
        "A nonempty independence attestation is provenance, not proof. This validator checks metadata",
        "and exact arithmetic only; genuine randomness and IID Bernoulli sampling remain external",
        "conditions for any future confidence statement.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = build_artifact()
    rendered = (
        json.dumps(artifact, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(artifact)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
