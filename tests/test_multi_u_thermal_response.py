from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from multi_u_thermal_response import (  # noqa: E402
    DOUBLING_EVEN,
    DOUBLING_ODD,
    EVEN_N_EXPONENT,
    FROZEN_U,
    ODD_N_EXPONENT,
    frozen_templates,
    n10_oracle,
    p49_descriptive,
    run_suite,
)


class MultiUThermalResponseTests(unittest.TestCase):
    def test_frozen_grid_and_monomials(self) -> None:
        self.assertEqual(FROZEN_U, (0.0, 0.025, 0.05))
        self.assertEqual(EVEN_N_EXPONENT, Fraction(-3, 4))
        self.assertEqual(ODD_N_EXPONENT, Fraction(-3, 8))
        self.assertEqual(DOUBLING_EVEN, Fraction(1, 2) ** Fraction(3, 4))
        self.assertEqual(DOUBLING_ODD, Fraction(1, 2) ** Fraction(3, 8))
        templates = frozen_templates()
        self.assertTrue(templates["do_not_add_levels_after_looking"])
        self.assertTrue(templates["not_a_target_for_p43_or_issue_57"])
        self.assertEqual(templates["u_grid"], [0.0, 0.025, 0.05])

    def test_n10_even_shape_vanishes(self) -> None:
        payload = n10_oracle()
        self.assertTrue(payload["even_shape_vanishes"])
        self.assertEqual(payload["Q"], 0.0)
        for u, row in payload["levels"].items():
            self.assertAlmostEqual(row["c"], 0.5, places=15, msg=str(u))

    def test_p49_descriptive_obeys_leading_u_monomials(self) -> None:
        for n in (130, 170):
            row = p49_descriptive(n)
            self.assertTrue(row["not_a_target_for_p43_or_issue_57"])
            self.assertEqual(row["role"], "development_descriptive_only")
            shape = row["shape"]
            self.assertAlmostEqual(shape["odd_u_ratio"], 1.0, places=2)
            self.assertAlmostEqual(shape["even_u_ratio"], 1.0, places=2)
            self.assertAlmostEqual(shape["w_over_u_0.025"] / shape["w_over_u_0.05"], 1.0, places=2)

    def test_suite_refuses_p43(self) -> None:
        payload = run_suite()
        self.assertTrue(payload["not_a_target_for_p43_or_issue_57"])
        freeze = (ROOT / "predictions" / "multi_u_thermal_response_20260829.yaml").read_text(
            encoding="utf-8"
        )
        self.assertIn("not_a_target_for_p43_or_issue_57", freeze)
        self.assertIn("do_not_add_levels_after_looking", freeze)
        self.assertIn("Do not score P43", freeze)


if __name__ == "__main__":
    unittest.main()
