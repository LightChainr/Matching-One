#!/usr/bin/env python3
"""Reproduce the frozen Issue #212 matching-odd post-reveal synthesis.

This scorer deliberately reads only two already-scored primary blocks and two
registered models per block.  It does not create a new prequential evidence
row: the output is a derived, post-reveal summary of independent raw-data
groups that are already represented in the evidence ledger.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


SCHEMA = "matching-one/issue212-matching-odd-synthesis/v1"
BLOCK_SPECS = (
    {
        "id": "issue43_n185_n265_deltaM",
        "raw_data_group": "issue43_n185_n265_500m_histograms",
        "h4_score": "H4_x21_over_4",
    },
    {
        "id": "issue57_norm5",
        "raw_data_group": "issue57_norm5_production",
        "h4_score": "H4_norm5",
    },
)


def chi_square_survival_even(chi_square: float, dof: int) -> float:
    """Return P(ChiSquare_dof >= chi_square) for positive even ``dof``."""
    if not math.isfinite(chi_square) or chi_square < 0.0:
        raise ValueError("chi-square must be finite and nonnegative")
    if not isinstance(dof, int) or dof <= 0 or dof % 2:
        raise ValueError("this exact scorer requires a positive even number of degrees of freedom")
    half_x = 0.5 * chi_square
    return math.exp(-half_x) * sum(
        half_x**order / math.factorial(order) for order in range(dof // 2)
    )


def _validated_score(block: dict, score_name: str) -> dict:
    try:
        score = block["scores"][score_name]
    except KeyError as exc:
        raise ValueError(f"block {block.get('id')} lacks registered score {score_name}") from exc
    required = ("chi_square", "dimension", "nlpd")
    if any(field not in score for field in required):
        raise ValueError(f"score {block['id']}/{score_name} lacks required fields")
    chi_square = float(score["chi_square"])
    nlpd = float(score["nlpd"])
    dimension = score["dimension"]
    if not math.isfinite(chi_square) or chi_square < 0.0:
        raise ValueError(f"score {block['id']}/{score_name} has invalid chi-square")
    if not isinstance(dimension, int) or isinstance(dimension, bool) or dimension <= 0:
        raise ValueError(f"score {block['id']}/{score_name} has invalid dimension")
    if not math.isfinite(nlpd):
        raise ValueError(f"score {block['id']}/{score_name} has invalid NLPD")
    return {"score_name": score_name, "chi_square": chi_square, "dof": dimension, "nlpd": nlpd}


def synthesize(ledger: dict, *, source_sha256: str | None = None) -> dict:
    """Validate and combine the two Issue #212 ledger blocks."""
    blocks = ledger.get("blocks")
    if not isinstance(blocks, list):
        raise ValueError("ledger blocks must be a list")
    by_id = {}
    for block in blocks:
        block_id = block.get("id")
        if block_id in by_id:
            raise ValueError(f"duplicate block id in ledger: {block_id}")
        by_id[block_id] = block

    selected = []
    raw_data_groups = []
    for spec in BLOCK_SPECS:
        if spec["id"] not in by_id:
            raise ValueError(f"ledger lacks frozen Issue #212 block {spec['id']}")
        block = by_id[spec["id"]]
        if block.get("role") != "primary":
            raise ValueError(f"block {spec['id']} is not primary")
        if block.get("status") != "SCORED":
            raise ValueError(f"block {spec['id']} is not SCORED")
        channel = block.get("channel", {})
        if channel.get("source") != "matching_odd" or channel.get("target") != "matching_odd":
            raise ValueError(f"block {spec['id']} is not matching_odd -> matching_odd")
        if block.get("raw_data_group") != spec["raw_data_group"]:
            raise ValueError(f"block {spec['id']} raw_data_group differs from the frozen contract")
        raw_data_groups.append(block["raw_data_group"])
        zero_score = _validated_score(block, "zero_effect")
        h4_score = _validated_score(block, spec["h4_score"])
        if zero_score["dof"] != h4_score["dof"]:
            raise ValueError(f"block {spec['id']} compares scores with unequal dimensions")
        selected.append({
            "id": spec["id"],
            "role": block["role"],
            "status": block["status"],
            "channel": {"source": channel["source"], "target": channel["target"]},
            "raw_data_group": block["raw_data_group"],
            "zero_effect": zero_score,
            "fixed_H4": h4_score,
        })

    if len(set(raw_data_groups)) != len(raw_data_groups):
        raise ValueError("frozen blocks do not have distinct raw_data_group values")

    joint = {}
    for model in ("zero_effect", "fixed_H4"):
        chi_square = sum(block[model]["chi_square"] for block in selected)
        dof = sum(block[model]["dof"] for block in selected)
        nlpd = sum(block[model]["nlpd"] for block in selected)
        survival_p = chi_square_survival_even(chi_square, dof)
        joint[model] = {
            "chi_square": chi_square,
            "dof": dof,
            "chi_square_survival_p": survival_p,
            "reject_at_alpha_0_05": survival_p < 0.05,
            "nlpd": nlpd,
        }

    delta_nlpd = joint["fixed_H4"]["nlpd"] - joint["zero_effect"]["nlpd"]
    return {
        "schema": SCHEMA,
        "issue": 212,
        "analysis_class": "post_reveal_synthesis",
        "source_ledger": {
            "schema_version": ledger.get("schema_version"),
            "manifest_version": ledger.get("manifest_version"),
            "sha256": source_sha256,
        },
        "selection_contract": {
            "block_ids": [spec["id"] for spec in BLOCK_SPECS],
            "required_role": "primary",
            "required_channel": {"source": "matching_odd", "target": "matching_odd"},
            "read_scores_only": {
                spec["id"]: ["zero_effect", spec["h4_score"]] for spec in BLOCK_SPECS
            },
        },
        "blocks": selected,
        "independence_contract": {
            "raw_data_groups": raw_data_groups,
            "raw_data_groups_distinct": True,
            "joint_covariance": "block_diagonal",
            "scope": "the distinct target streams/seeds frozen by Issue #212; no other ledger views enter",
        },
        "joint_scores": joint,
        "predictive_comparison": {
            "delta_nlpd_fixed_H4_minus_zero_effect": delta_nlpd,
            "preferred_lower_nlpd": "fixed_H4" if delta_nlpd < 0.0 else "zero_effect",
        },
        "governance": {
            "adds_new_primary_evidence": False,
            "may_be_appended_as_primary_ledger_row": False,
            "interpretation": "derived synthesis of two already-primary independent target blocks",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    raw = args.ledger.read_bytes()
    result = synthesize(json.loads(raw), source_sha256=hashlib.sha256(raw).hexdigest())
    result["source_ledger"]["path"] = args.ledger.as_posix()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
