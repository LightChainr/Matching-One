#!/usr/bin/env python3
"""Reclassify the P250 reveal support-first and audit normalization fields."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from score_z5_charged_threepoint import zero_score


def header(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        return next(csv.reader(handle))


def audit(score_path: Path, cubic_batches: Path, charged_archive_batches: Path) -> dict:
    score = json.loads(score_path.read_text())
    support = zero_score(score["primary_point"], score["primary_covariance_8x8"])
    cubic_fields = header(cubic_batches)
    charged_fields = header(charged_archive_batches)
    required = {
        "same_replica_local_pair_products": False,
        "separation_label": False,
        "local_cubic_pair_cross_covariance": False,
    }
    return {
        "schema": "matching-one/p250-near-zero-reveal-audit/v1",
        "status": "existing_archives_reanalyzed",
        "existing_1m_cubic_support": {
            **support,
            "decision_at_0.05": (
                "detected" if support["survival_p"] < 0.05 else "not_detected"
            ),
        },
        "existing_phase_closure": {
            **score["cross_product_closure"],
            "reclassified_status": "not_interpretable_until_nonzero_support",
        },
        "field_inventory": {
            "be80f25_threepoint_batch_fields": cubic_fields,
            "P226_charged_batch_fields": charged_fields,
            "P226_semantics": (
                "global marked-row one-point response under a different root schedule; "
                "not the local O_r(x) used by the cubic"
            ),
        },
        "normalization_requirements_present": required,
        "can_construct_same_operator_separation_ratio": all(required.values()),
        "exact_reason": (
            "The retained cubic batches contain only products after the three local rows "
            "were multiplied. P226 contains a different global charged response. Neither "
            "archive retains O_r(x)O_-r(y), a separation, or its joint covariance with C_rst."
        ),
        "consequence": (
            "Do not add replicas to the compact cubic. Acquire pair products and cubics "
            "jointly at multiple separations, and gate phase closure behind cubic support."
        ),
    }


def render(result: dict) -> str:
    support = result["existing_1m_cubic_support"]
    closure = result["existing_phase_closure"]
    return "\n".join(
        [
            "# P250 existing-archive near-zero audit",
            "",
            f"- cubic support: `{support['chi_square']}/{support['degrees_of_freedom']}`, p `{support['survival_p']}` -> `{support['decision_at_0.05']}`",
            f"- phase closure: `{closure['chi_square']}/{closure['degrees_of_freedom']}`, p `{closure['survival_p']}` -> `{closure['reclassified_status']}`",
            f"- same-operator separation ratio recoverable: `{result['can_construct_same_operator_separation_ratio']}`",
            "",
            result["exact_reason"],
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("score", type=Path)
    parser.add_argument("cubic_batches", type=Path)
    parser.add_argument("charged_archive_batches", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.score, args.cubic_batches, args.charged_archive_batches)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
