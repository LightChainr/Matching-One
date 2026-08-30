#!/usr/bin/env python3
"""Validate the data-free exact commuting-square contract for Issue 158."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Tuple

from gaussian_harmonic_arithmetic import gmul, harmonic, norm


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "gaussian_commuting_square_root_contract.json"
EXPECTED_FACTORS = {"norm2": Fraction(-1, 4), "norm5": Fraction(-14, 625), "norm10": Fraction(7, 1250)}
FORBIDDEN_DATA_KEYS = frozenset(
    {"samples", "seed", "counter", "roots", "observed_ratio", "covariance", "score", "chi2", "p_value"}
)
Pair = Tuple[int, int]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _walk_forbidden(value: Any, path: str = "$") -> None:
    if isinstance(value, Mapping):
        bad = sorted(FORBIDDEN_DATA_KEYS.intersection(value))
        _require(not bad, "%s contains target-data fields: %s" % (path, ",".join(bad)))
        for key, child in value.items():
            _walk_forbidden(child, "%s.%s" % (path, key))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _walk_forbidden(child, "%s[%d]" % (path, index))


def _pair(value: Any, label: str) -> Pair:
    _require(isinstance(value, list) and len(value) == 2, "%s must be a pair" % label)
    _require(all(isinstance(component, int) for component in value), "%s must contain integers" % label)
    pair = value[0], value[1]
    _require(norm(pair) > 0, "%s must be nonzero" % label)
    return pair


def canonical_d4(value: Pair) -> Pair:
    """Return the deterministic nonnegative, descending D4 representative."""

    return tuple(sorted((abs(value[0]), abs(value[1])), reverse=True))  # type: ignore[return-value]


def delta_cos4(pair: Sequence[Pair]) -> Fraction:
    _require(len(pair) == 2, "an ordered orientation pair is required")
    return harmonic(pair[0], 1)[0] - harmonic(pair[1], 1)[0]


def root_factor(parent: Sequence[Pair], multiplier: Pair) -> Fraction:
    parent_delta = delta_cos4(parent)
    _require(parent_delta != 0, "parent H4 contrast vanishes")
    child = [gmul(parent[0], multiplier), gmul(parent[1], multiplier)]
    angular_ratio = delta_cos4(child) / parent_delta
    return angular_ratio / (norm(multiplier) ** 2)


def validate_contract(contract: Mapping[str, Any], arithmetic_bytes: bytes) -> dict[str, Any]:
    _walk_forbidden(contract)
    _require(contract.get("schema") == "matching-one/gaussian-commuting-square-root-contract/v1", "unknown schema")
    _require(contract.get("issue") == 158, "wrong issue")
    _require(contract.get("status") == "exact_design_only_no_target_data", "contract is not data-free")
    source = contract.get("arithmetic_source", {})
    _require(source.get("path") == "scripts/gaussian_harmonic_arithmetic.py", "wrong arithmetic source")
    digest = hashlib.sha256(arithmetic_bytes).hexdigest()
    _require(source.get("sha256") == digest, "arithmetic source SHA-256 mismatch")

    character = contract.get("character_contract", {})
    _require(character.get("level") == 4, "character level drift")
    _require(character.get("fit_parameters") == 0, "fit parameter introduced")
    for key, expected in EXPECTED_FACTORS.items():
        _require(Fraction(character.get("%s_factor" % key)) == expected, "%s factor drift" % key)

    lineages = contract.get("lineages")
    _require(isinstance(lineages, list) and len(lineages) == 2, "exactly two lineages required")
    audited = []
    ids = []
    for index, lineage in enumerate(lineages):
        label = "lineage[%d]" % index
        lineage_id = lineage.get("id")
        _require(isinstance(lineage_id, str) and lineage_id, "%s id missing" % label)
        ids.append(lineage_id)
        parent = [_pair(value, "%s parent" % label) for value in lineage.get("parent", [])]
        _require(len(parent) == 2, "%s must have two ordered parents" % label)
        parent_n = lineage.get("parent_N")
        _require(norm(parent[0]) == parent_n and norm(parent[1]) == parent_n, "%s parent norm drift" % label)
        multipliers = {name: _pair(value, "%s %s" % (label, name)) for name, value in lineage.get("multipliers", {}).items()}
        _require(set(multipliers) == {"norm2", "norm5", "norm10"}, "%s multiplier set drift" % label)
        _require(norm(multipliers["norm2"]) == 2 and norm(multipliers["norm5"]) == 5, "%s generator norm drift" % label)
        product_25 = gmul(multipliers["norm2"], multipliers["norm5"])
        product_52 = gmul(multipliers["norm5"], multipliers["norm2"])
        _require(product_25 == product_52 == multipliers["norm10"], "%s commuting product drift" % label)
        _require(norm(multipliers["norm10"]) == 10, "%s direct multiplier norm drift" % label)

        children = lineage.get("canonical_children", {})
        factors = {}
        for name in ("norm2", "norm5", "norm10"):
            actual_children = [canonical_d4(gmul(value, multipliers[name])) for value in parent]
            stored_children = [_pair(value, "%s %s child" % (label, name)) for value in children.get(name, [])]
            _require(actual_children == stored_children, "%s %s canonical child drift" % (label, name))
            factors[name] = root_factor(parent, multipliers[name])
            _require(factors[name] == EXPECTED_FACTORS[name], "%s %s root factor drift" % (label, name))
        direct_raw = [gmul(value, multipliers["norm10"]) for value in parent]
        path_25 = [gmul(gmul(value, multipliers["norm2"]), multipliers["norm5"]) for value in parent]
        path_52 = [gmul(gmul(value, multipliers["norm5"]), multipliers["norm2"]) for value in parent]
        _require(direct_raw == path_25 == path_52, "%s descendant paths do not close" % label)
        _require(factors["norm2"] * factors["norm5"] == factors["norm10"], "%s character factors do not compose" % label)
        _require(Fraction(lineage.get("direct_root_factor")) == factors["norm10"], "%s stored direct factor drift" % label)
        _require(norm(direct_raw[0]) == lineage.get("target_N") and norm(direct_raw[1]) == lineage.get("target_N"), "%s target norm drift" % label)
        audited.append({"id": lineage_id, "target_N": lineage["target_N"], "direct_root_factor": str(factors["norm10"]), "paths_close": True})
    _require(len(ids) == len(set(ids)), "duplicate lineage ids")

    protocol = contract.get("scoring_protocol", {})
    for key in (
        "score_fixed_ratio_before_finite_size_corrections",
        "direct_and_two_step_paths_form_one_correlated_block",
        "freeze_corrections_from_source_data_before_target_read",
    ):
        _require(protocol.get(key) is True, "scoring protocol drift: %s" % key)
    _require(protocol.get("claim_depends_on_crn") is False, "claim cannot depend on CRN")
    _require(contract.get("claim_boundary", {}).get("parent_issue") == "remain open", "parent boundary changed")
    return {
        "schema": contract["schema"],
        "status": "valid_exact_design_only",
        "arithmetic_source_sha256": digest,
        "lineages": audited,
        "common_root_factor": "7/1250",
        "contains_target_data": False,
        "parent_issue": "remain open",
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    args = parser.parse_args(argv)
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    source_path = ROOT / contract["arithmetic_source"]["path"]
    result = validate_contract(contract, source_path.read_bytes())
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
