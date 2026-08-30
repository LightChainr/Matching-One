#!/usr/bin/env python3
"""Tests for the Issue #40 joint-Gram production scorer."""

from __future__ import annotations

import sys
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_p40_production_motif_projection as scorer  # noqa: E402


NAMES = [
    "q", "C_black", "C_white", "V", "E", "F0", "nnn_pos", "nnn_neg",
    "path3_x", "path3_y", "corners", "E_mc", "F0_mc", "nnn_pos_mc",
    "nnn_neg_mc", "path3_x_mc", "path3_y_mc", "corners_mc", "right_angle",
    "right_angle_mc",
]


def geometry_payload(a: int, b: int, vectors: list[list[float]]) -> dict:
    width = len(NAMES)
    return {
        "a": a,
        "b": b,
        "sum": [sum(row[i] for row in vectors) for i in range(width)],
        "gram": [[sum(row[i] * row[j] for row in vectors) for j in range(width)] for i in range(width)],
        "wrapping_l1": 0,
        "identity_l1": 0,
    }


def synthetic_row(batch: int) -> dict:
    vectors = []
    for replica in range(20):
        t = batch * 20 + replica
        controls = [
            (t % 7) - 3,
            ((3 * t + 1) % 11) - 5,
            ((5 * t + 2) % 13) - 6,
            ((7 * t + 4) % 17) - 8,
        ]
        noise = ((11 * t + 3) % 19 - 9) * 0.03
        row = [0.0] * len(NAMES)
        row[0] = 1.4 * controls[0] - 0.8 * controls[1] + 0.5 * controls[2] + 0.2 * controls[3] + noise
        row[11] = controls[0]
        row[13] = controls[1]
        row[14] = 0.0
        row[12] = controls[2]
        row[19] = controls[3]
        vectors.append(row)
    zeros = [[0.0] * len(NAMES) for _ in vectors]
    return {
        "n": 65,
        "batch": batch,
        "samples": len(vectors),
        "p_ref": 0.592746050790,
        "names": NAMES,
        "first": geometry_payload(8, 1, vectors),
        "second": geometry_payload(7, 4, zeros),
        "cross_gram_semantics": "sum first[i]*second[j] over same replicas",
        "cross_gram": [[0.0] * len(NAMES) for _ in NAMES],
    }


class ProductionMotifProjectionTest(unittest.TestCase):
    def test_joint_gram_recovers_delta_and_heldout_reduction(self) -> None:
        moments = [scorer.batch_moments(synthetic_row(batch)) for batch in range(10)]
        result = scorer.score_size(moments)
        self.assertEqual(result["labels"], list(scorer.LABELS))
        self.assertEqual(result["batches"], 10)
        self.assertGreater(result["cross_fit"]["oof_variance_reduction"], 100.0)
        self.assertEqual(result["exact_gates"]["identity_l1"], 0.0)

    def test_missing_joint_gram_is_rejected(self) -> None:
        row = synthetic_row(0)
        del row["cross_gram"]
        with self.assertRaisesRegex(ValueError, "joint first-second"):
            scorer.batch_moments(row)


if __name__ == "__main__":
    unittest.main()
