from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p250_cross_scale_power_freeze import freeze  # noqa: E402
from score_z5_projective_leg_cross_scale import scalar_gls  # noqa: E402
from z5_projective_leg_cross_scale_mc import FIELD_ORDER, exact_gate, run  # noqa: E402


class Z5ProjectiveLegCrossScaleTests(unittest.TestCase):
    def test_n101_is_first_five_separation_oblique_parent(self) -> None:
        gate = exact_gate()
        self.assertTrue(gate["passed"])
        self.assertEqual(gate["parent_order"], 101)
        self.assertEqual(gate["child_order"], 505)
        self.assertEqual(gate["distinct_signed_axis_displacements"], 20)
        self.assertEqual(gate["candidate_audit"]["candidates"][-1]["norm"], 101)

    def test_tiny_stream_has_80_coordinate_covariance(self) -> None:
        rows, analysis = run(4, 2, 1, 0.59274605079, 25050510120260830, 0, cap=2_000)
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(FIELD_ORDER), 80)
        self.assertEqual(len(analysis["covariance_of_mean"]), 80)

    def test_scalar_gls_recovers_exact_coefficient(self) -> None:
        design = [1.0, 2.0, 4.0]
        point = [1.25 * value for value in design]
        covariance = [[0.01 if i == j else 0.0 for j in range(3)] for i in range(3)]
        fit = scalar_gls(point, covariance, design)
        self.assertAlmostEqual(fit["coefficient"], 1.25, places=14)
        self.assertAlmostEqual(fit["chi_square"], 0.0, places=14)

    def test_existing_variance_selects_40k(self) -> None:
        path = ROOT / "results/huawei-20260830/P250-z5-projective-leg-pair-transfer-40k/response_40k.batches.csv"
        result = freeze(path)
        self.assertEqual(result["selected_samples"], 40_000)
        selected = next(row for row in result["grid"] if row["samples"] == 40_000)
        self.assertGreater(selected["projected_minimum_d1_d5_real_abs_z"], 6.0)


if __name__ == "__main__":
    unittest.main()
