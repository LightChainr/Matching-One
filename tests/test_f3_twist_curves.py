from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_f3_twist_curves import (  # noqa: E402
    character_coefficients,
    derivative,
    evaluate,
    score,
    twist_coefficients,
)


class F3TwistCurveTests(unittest.TestCase):
    def test_constraint_and_character_polynomials(self) -> None:
        base = {
            "P0": [1.0, 0.0, 0.0],
            "P2": [0.0, 0.0, 1.0],
            "L_1_0": [0.0, 0.4, 0.0],
            "L_1_1": [0.0, 0.1, 0.0],
            "L_1_2": [0.0, 0.2, 0.0],
            "L_0_1": [0.0, 0.3, 0.0],
        }
        twists = twist_coefficients(base)
        self.assertEqual(twists["T_0_0"], [1.0, 1.0, 1.0])
        self.assertEqual(twists["T_1_0"], twists["T_2_0"])
        self.assertTrue(all(
            curve[k + 1] <= curve[k]
            for name, curve in twists.items() if name != "T_0_0"
            for k in range(len(curve) - 1)
        ))
        characters = character_coefficients(base)
        self.assertAlmostEqual(characters["H4_axis_diag"][1], 0.2)
        self.assertAlmostEqual(evaluate(twists["T_0_0"], 0.4), 1.0)
        self.assertAlmostEqual(derivative(twists["T_0_0"], 0.4), 0.0)

    def test_n65_full_curve_and_crossing_selector(self) -> None:
        source = ROOT / "results/local-20260830/P334-projective-birth-N65-smoke"
        result = score(
            source / "n65_20k.births.csv",
            source / "n65_20k.metadata.json",
            0.592746050790,
        )
        self.assertTrue(result["exact_gates"]["passed"])
        self.assertEqual(result["exact_gates"]["maximum_zero_twist_residual"], 0.0)
        self.assertEqual(
            result["exact_gates"]["maximum_nonzero_twist_monotonicity_violation"],
            0.0,
        )
        coefficients = result["bernstein_batch_coefficients"]
        self.assertEqual(len(coefficients["values"]), 20)
        self.assertEqual(len(coefficients["values"][0]["orientations"]["first"]["P0"]), 66)
        grid_row = result["evaluated_grid"]["joint_estimates"][0]
        self.assertEqual(len(grid_row["order"]), 48)
        self.assertEqual(len(grid_row["covariance"]), 48)
        self.assertTrue(all(len(row) == 48 for row in grid_row["covariance"]))
        selected = result["crossing_selector"]["selected"]
        self.assertEqual(selected["character"], "H4_axis_diag")
        self.assertAlmostEqual(selected["root"], 0.5736333260803443)
        self.assertEqual(selected["leave_one_root_fraction"], 1.0)
        self.assertGreater(selected["jackknife_standard_error"], 0.02)


if __name__ == "__main__":
    unittest.main()
