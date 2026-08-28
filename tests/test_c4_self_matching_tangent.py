from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from c4_self_matching_tangent import (  # noqa: E402
    center_score_function_response,
    exact_tangent_report,
)


class C4SelfMatchingTangentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.report = exact_tangent_report(3, 1)

    def test_family_and_selection_rule(self) -> None:
        self.assertEqual(
            self.report["family"]["matching_exchange"],
            "(t,lambda)->(-t,-lambda)",
        )
        self.assertIn("m+n is even", self.report["selection_rule"])

    def test_response_matrix_is_channel_independent(self) -> None:
        for row in self.report["channels"].values():
            self.assertEqual(
                row["response_matrix_rows_Rplus_Rminus_columns_t_lambda"],
                [["0", "0"], ["15/8", "5/4"]],
            )

    def test_score_function_oracle_reproduces_exact_n10_derivatives(self) -> None:
        oracle = self.report["score_function_oracle"]
        self.assertEqual(oracle["identity"], "d_u <O> = <O*S_u>")
        self.assertTrue(oracle["verified_against_bivariate_polynomial"])
        for row in self.report["channels"].values():
            score = row["score_function_center_response_Rminus"]
            self.assertEqual(score, {"d_t": "15/8", "d_lambda": "5/4"})
            self.assertEqual(
                [score["d_t"], score["d_lambda"]],
                row["response_matrix_rows_Rplus_Rminus_columns_t_lambda"][1],
            )

    def test_score_function_signs_on_two_site_synthetic_observable(self) -> None:
        # One even and one odd site.  O=1 exactly when the even site is occupied,
        # so d_t <O> = d_lambda <O> = 1 at the center.
        counts = [[0, 0], [1, 1]]
        self.assertEqual(
            center_score_function_response(counts),
            {"d_t": Fraction(1), "d_lambda": Fraction(1)},
        )

    def test_only_physical_odd_root_is_self_matching_center(self) -> None:
        gate = self.report["exact_odd_root_gate"]
        self.assertEqual(gate["factorization"], "lambda*(5/4 - 4*lambda^4)")
        self.assertEqual(gate["only_legal_root"], "0")
        self.assertGreater(gate["other_real_root_magnitudes"], 0.5)
        self.assertFalse(gate["other_roots_are_legal"])

    def test_even_h4_scan_is_frozen_and_has_holdout(self) -> None:
        protocol = self.report["nontrivial_improved_action_protocol"]
        self.assertEqual(protocol["minimum_orientation_design"]["N"], 130)
        self.assertEqual(protocol["frozen_nonnegative_scan"], ["0", "1/8", "1/4", "3/8"])
        self.assertIn("lack-of-fit", protocol["fit"])


if __name__ == "__main__":
    unittest.main()
