from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_fk_q_score_oracle import render  # noqa: E402


class ExactFKQScoreOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = render(2)

    def test_exponential_family_statistic_is_exact(self) -> None:
        self.assertEqual(self.payload["geometry"]["configurations"], 256)
        self.assertTrue(self.payload["exact_means"]["T_equals_k_plus_b_over_2"])
        self.assertEqual(sum(row["count"] for row in self.payload["sqrtQ_histogram"]), 256)

    def test_three_score_orders_match_direct_ratio_derivatives(self) -> None:
        for observable in self.payload["observable_derivatives"].values():
            self.assertTrue(observable["direct_equals_score"])

    def test_dual_wrap_baseline_and_nontrivial_q_tangent(self) -> None:
        derivatives = self.payload["observable_derivatives"]
        self.assertEqual(
            derivatives["open_wrap"]["expectation_at_Q1"],
            derivatives["closed_dual_wrap"]["expectation_at_Q1"],
        )
        self.assertNotEqual(
            derivatives["wrap_difference"]["t_derivatives_direct_ratio"][0],
            "0/1",
        )

    def test_wrap_difference_polynomial_has_simple_q1_zero(self) -> None:
        polynomial = self.payload["wrap_difference_polynomial"]
        self.assertTrue(polynomial["Q1_factorization_exact"])
        self.assertTrue(polynomial["Q1_zero_is_simple"])
        self.assertTrue(polynomial["P_L_coefficients_strictly_positive"])
        self.assertEqual(polynomial["tangent_from_simple_zero"], "69/256")

    def test_committed_l3_has_same_simple_positive_quotient_structure(self) -> None:
        payload = json.loads(
            (ROOT / "results/fk-q-score/L3.json").read_text(encoding="utf-8")
        )
        polynomial = payload["wrap_difference_polynomial"]
        self.assertEqual(payload["geometry"]["configurations"], 262144)
        self.assertTrue(polynomial["Q1_factorization_exact"])
        self.assertTrue(polynomial["Q1_zero_is_simple"])
        self.assertTrue(polynomial["P_L_coefficients_strictly_positive"])
        self.assertEqual(polynomial["tangent_from_simple_zero"], "18865/65536")

    def test_committed_artifact_is_byte_reproducible(self) -> None:
        artifact = json.loads(
            (ROOT / "results/fk-q-score/latest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(artifact, self.payload)


if __name__ == "__main__":
    unittest.main()
