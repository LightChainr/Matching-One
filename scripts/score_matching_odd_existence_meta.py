#!/usr/bin/env python3
"""Post-reveal target-only existence synthesis for independent matching-odd blocks.

This script does not create a new prequential evidence row. It reads two
already-scored primary blocks, verifies that their raw-data groups are distinct
and their target channel is matching-odd, then combines the registered
zero-effect and fixed-H4 goodness-of-fit statistics using block-diagonal
independence.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


BLOCKS = (
    ("issue43_n185_n265_deltaM", "H4_x21_over_4"),
    ("issue57_norm5", "H4_norm5"),
)


def chi_square_survival_even(chi_square: float, degrees: int) -> float:
    """Exact chi-square survival function for positive even degrees of freedom."""
    if chi_square < 0 or degrees <= 0 or degrees % 2:
        raise ValueError("requires chi_square>=0 and positive even degrees")
    x = chi_square / 2.0
    return math.exp(-x) * sum(
        x**order / math.factorial(order) for order in range(degrees // 2)
    )


def index_blocks(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    rows = payload.get("blocks")
    if not isinstance(rows, list):
        raise ValueError("evidence ledger has no block list")
    indexed = {row["id"]: row for row in rows}
    if len(indexed) != len(rows):
        raise ValueError("duplicate evidence block id")
    return indexed


def combine(payload: dict[str, Any]) -> dict[str, Any]:
    indexed = index_blocks(payload)
    selected = []
    raw_groups = set()

    for block_id, h4_model in BLOCKS:
        block = indexed.get(block_id)
        if block is None:
            raise ValueError(f"missing block {block_id}")
        if block.get("role") != "primary" or block.get("status") != "SCORED":
            raise ValueError(f"{block_id} is not a scored primary block")
        channel = block.get("channel", {})
        if channel.get("target") != "matching_odd":
            raise ValueError(f"{block_id} target channel is not matching_odd")
        raw_group = block.get("raw_data_group")
        if not raw_group or raw_group in raw_groups:
            raise ValueError("selected blocks do not have distinct raw-data groups")
        raw_groups.add(raw_group)
        scores = block.get("scores", {})
        if "zero_effect" not in scores or h4_model not in scores:
            raise ValueError(f"{block_id} lacks the frozen zero/H4 score")
        selected.append((block_id, raw_group, h4_model, scores))

    result: dict[str, Any] = {
        "schema": "matching-one/independent-target-existence-meta/v1",
        "status": "post-reveal synthesis; not an additive prequential evidence row",
        "blocks": [],
        "models": {},
        "governance": {
            "target_channel": "matching_odd",
            "independence_rule": "distinct raw_data_group; cross-block target covariance set to zero",
            "anti_double_counting": "use exactly the two registered primary target blocks; add no derived views",
        },
    }

    for block_id, raw_group, h4_model, scores in selected:
        result["blocks"].append(
            {
                "id": block_id,
                "raw_data_group": raw_group,
                "h4_model": h4_model,
                "zero": scores["zero_effect"],
                "h4": scores[h4_model],
            }
        )

    for label, key in (("zero_effect", "zero_effect"), ("H4_fixed_predictions", None)):
        chi_square = 0.0
        dimensions = 0
        nlpd = 0.0
        for block_id, _raw_group, h4_model, scores in selected:
            model_key = key if key is not None else h4_model
            row = scores[model_key]
            chi_square += float(row["chi_square"])
            dimensions += int(row["dimension"])
            nlpd += float(row["nlpd"])
        result["models"][label] = {
            "chi_square": chi_square,
            "degrees_of_freedom": dimensions,
            "chi_square_survival": chi_square_survival_even(chi_square, dimensions),
            "nlpd": nlpd,
        }

    h4 = result["models"]["H4_fixed_predictions"]
    zero = result["models"]["zero_effect"]
    result["comparison"] = {
        "delta_nlpd_H4_minus_zero": h4["nlpd"] - zero["nlpd"],
        "interpretation": (
            "negative favors the block-specific fixed H4 predictions; "
            "do not treat the chi-square difference as a calibrated likelihood-ratio test"
        ),
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "ledger", nargs="?", type=Path, default=Path("results/evidence-ledger/latest.json")
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.ledger.read_text(encoding="utf-8"))
    result = combine(payload)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
