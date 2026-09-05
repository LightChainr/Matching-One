
from __future__ import annotations
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import digital_alexander_rank_oracle as oracle  # noqa: E402


class DigitalAlexanderRankOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        config = json.loads(
            (ROOT / "analysis" / "digital_alexander_rank_manifest.json").read_text()
        )
        cls.result = oracle.analyze(config)

    def test_elementary_rank_lemma_is_exhaustive(self) -> None:
        lemma = self.result["elementary_rank_lemma"]
        self.assertEqual(len(lemma["rows"]), 9)
        self.assertTrue(lemma["premise_implies_weak_identity"])
        self.assertEqual(
            lemma["premise_rank_pairs"],
            [[0, 0], [0, 2], [1, 1], [2, 0], [2, 2]],
        )

    def test_every_declared_exhaustive_geometry_has_no_weak_failure(self) -> None:
        rows = self.result["exhaustive_geometries"]
        self.assertEqual(len(rows), 5)
        for row in rows:
            with self.subTest(geometry=row["id"]):
                self.assertFalse(row["common_channel_counterexamples"])
                self.assertFalse(row["weak_rank_counterexamples"])
                self.assertFalse(row["strong_rank_sum_counterexamples"])
                self.assertEqual(
                    row["raw_rank_equality_count"], row["q_zero_count"]
                )
                self.assertEqual(
                    row["raw_rank_equality_outside_q_zero_count"], 0
                )

    def test_expected_exact_configuration_counts(self) -> None:
        counts = {
            row["id"]: row["configurations"]
            for row in self.result["exhaustive_geometries"]
        }
        self.assertEqual(
            counts,
            {
                "axis-L2": 16,
                "axis-L3": 512,
                "gaussian-2-1": 32,
                "diamond-L2": 256,
                "c4-self-matching-3-1": 1024,
            },
        )

    def test_empty_and_full_configuration_normalization(self) -> None:
        geometry = oracle.axis_integer_torus(3)
        empty = oracle.rank_identity_record(geometry, (False,) * geometry.n)
        full = oracle.rank_identity_record(geometry, (True,) * geometry.n)
        self.assertEqual(
            (empty["q_either"], empty["rank_black"], empty["rank_white"]),
            (-1, 0, 2),
        )
        self.assertEqual(
            (full["q_either"], full["rank_black"], full["rank_white"]),
            (1, 2, 0),
        )

    def test_counterexample_search_is_frozen_and_weak_null(self) -> None:
        search = self.result["counterexample_search"]
        self.assertEqual(search["seed"], 2690829)
        self.assertEqual(search["matrix_count"], 160)
        self.assertGreaterEqual(search["configurations_evaluated"], 50000)
        self.assertFalse(search["common_channel_counterexamples"])
        self.assertFalse(search["weak_rank_counterexamples"])
        self.assertFalse(search["strong_rank_sum_counterexamples"])

    def test_c4_control_uses_self_matching_edges(self) -> None:
        row = next(
            row
            for row in self.result["exhaustive_geometries"]
            if row["id"] == "c4-self-matching-3-1"
        )
        self.assertTrue(row["matching_edges_equal_primal_edges"])


if __name__ == "__main__":
    unittest.main()
