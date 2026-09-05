
from __future__ import annotations
import math
import sys
import unittest
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from intrinsic_quantile_center import (  # noqa: E402
    DOUBLING_RATIO,
    FROZEN_U,
    Q_EXPONENT,
    WIDTH_EXPONENT,
    descriptive_p49,
    frozen_monomials,
    midpoint_difference,
    n10_matching,
    n10_oracle,
    quantile_levels,
)


class IntrinsicQuantileCenterTests(unittest.TestCase):
    def test_frozen_levels_and_exponents_are_locked(self) -> None:
        self.assertEqual(FROZEN_U, (0.025, 0.05))
        self.assertEqual(Q_EXPONENT, Fraction(-3, 4))
        self.assertEqual(WIDTH_EXPONENT, Fraction(3, 8))
        self.assertAlmostEqual(DOUBLING_RATIO, 2 ** (-0.75), places=15)
        freeze = frozen_monomials()
        self.assertTrue(freeze["do_not_add_levels_after_looking"])
        self.assertTrue(freeze["not_a_numerical_target_for_p43_or_issue_57"])
        self.assertTrue(freeze["p49_130_170_are_not_a_doubling_pair"])
        with self.assertRaisesRegex(ValueError, "do not add levels"):
            quantile_levels(n10_matching, levels=(0.025, 0.05, 0.1))

    def test_n10_beta33_midpoints_are_exactly_one_half(self) -> None:
        oracle = n10_oracle()
        self.assertTrue(oracle["c_u_equals_half"])
        self.assertTrue(oracle["Q_vanishes_by_oddness"])
        self.assertEqual(oracle["N"], 10)
        self.assertAlmostEqual(n10_matching(0.5), 0.0, places=15)
        levels = quantile_levels(n10_matching)
        for u, level in levels.items():
            self.assertAlmostEqual(n10_matching(level.p_minus), -u, places=12)
            self.assertAlmostEqual(n10_matching(level.p_plus), u, places=12)
            self.assertAlmostEqual(level.c, 0.5, places=15)
        self.assertEqual(midpoint_difference(levels), oracle["Q"])
        self.assertAlmostEqual(oracle["Q"], 0.0, places=15)

    def test_prediction_yaml_forbids_p43_and_new_levels(self) -> None:
        text = (
            ROOT / "predictions" / "intrinsic_quantile_center_20260829.yaml"
        ).read_text(encoding="utf-8")
        self.assertIn("0.025", text)
        self.assertIn("0.05", text)
        self.assertIn("not_a_numerical_target_for_p43_or_issue_57", text)
        self.assertIn("2^{-3/4}", text)
        self.assertIn("do_not_add_levels_after_looking", text)
        self.assertNotIn("0.075", text)
        self.assertNotIn("P43 target", text)

    def test_descriptive_p49_is_labeled_and_solves_frozen_levels(self) -> None:
        payload = descriptive_p49()
        self.assertTrue(payload["not_a_numerical_target_for_p43_or_issue_57"])
        self.assertTrue(payload["p43_was_not_read_as_a_target"])
        self.assertEqual(payload["p49_role"], "development_descriptive_only")
        for n, expected_q, expected_scaled in (
            ("130", -4.035e-6, -1.554e-4),
            ("170", -3.355e-6, -1.580e-4),
        ):
            row = payload["p49_descriptive"][n]
            self.assertTrue(row["not_a_numerical_target_for_p43_or_issue_57"])
            self.assertEqual(row["role"], "development_descriptive_only")
            self.assertEqual(row["channel"], "rank-2 cross")
            self.assertEqual(row["batches"], 100)
            self.assertEqual(row["samples"], 100_000_000)
            self.assertAlmostEqual(row["Q"], expected_q, delta=2e-8)
            self.assertAlmostEqual(
                row["Q_times_N_to_3_over_4"], expected_scaled, delta=1e-6
            )
            for u in ("0.025", "0.05"):
                level = row["levels"][u]
                self.assertAlmostEqual(level["M_minus"], -float(u), places=10)
                self.assertAlmostEqual(level["M_plus"], float(u), places=10)
            # Thermal widths scale more cleanly than Q_N on these two children.
            self.assertAlmostEqual(
                row["w_0.025_times_N_to_3_over_8"], 0.01431, places=4
            )
            self.assertAlmostEqual(
                row["w_0.05_times_N_to_3_over_8"], 0.02864, places=4
            )
            self.assertAlmostEqual(
                row["w_0.05_times_N_to_3_over_8"]
                / row["w_0.025_times_N_to_3_over_8"],
                2.0,
                places=2,
            )
            # 130 and 170 are not a doubling pair, so do not test 2^{-3/4} here.
            self.assertNotAlmostEqual(
                abs(payload["p49_descriptive"]["170"]["Q"]
                    / payload["p49_descriptive"]["130"]["Q"]),
                DOUBLING_RATIO,
                places=2,
            )
        self.assertAlmostEqual(
            payload["p49_descriptive"]["130"]["Q_times_N_to_3_over_4"]
            / payload["p49_descriptive"]["170"]["Q_times_N_to_3_over_4"],
            1.0,
            delta=0.03,
        )
        self.assertTrue(math.isfinite(DOUBLING_RATIO))


if __name__ == "__main__":
    unittest.main()
