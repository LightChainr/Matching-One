from __future__ import annotations

import sys
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p50_fullcurve_n290 import independent_residual_covariance, scalar_score  # noqa: E402


class P50FullcurveExecutionTests(unittest.TestCase):
    def test_execution_plan_is_balanced_and_independent(self) -> None:
        plan = yaml.safe_load(
            (ROOT / "experiments/p50_n145_n290_fullcurve_20260829.yaml").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(plan["status"], "production_ready_before_N290_fullcurve_reveal")
        self.assertEqual(plan["sampling"]["samples_per_size"], 100000000)
        self.assertEqual(plan["sampling"]["batches"], 100)
        self.assertEqual(plan["sampling"]["cross_size_rng"], "independent")
        self.assertNotEqual(plan["parent_N145"]["seed"], plan["child_N290"]["seed"])
        self.assertEqual(plan["parent_N145"]["representations"], [[12, 1], [9, 8]])
        self.assertEqual(plan["child_N290"]["lineage_order"], [[13, 11], [17, 1]])

    def test_n290_period_matrices_have_order_290(self) -> None:
        plan = yaml.safe_load(
            (ROOT / "experiments/p50_n145_n290_fullcurve_20260829.yaml").read_text(
                encoding="utf-8"
            )
        )
        for key in ("first_period_matrix", "second_period_matrix"):
            matrix = plan["child_N290"][key]
            determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
            self.assertEqual(determinant, 290)
        self.assertEqual(plan["child_N290"]["determinant_order"], 290)

    def test_independent_residual_covariance_adds_size_local_terms(self) -> None:
        parent = [[4.0, 1.0], [1.0, 9.0]]
        child = [[16.0, 2.0], [2.0, 25.0]]
        result = independent_residual_covariance(parent, child, -0.5)
        self.assertEqual(result, [[17.0, 2.25], [2.25, 27.25]])

    def test_scalar_score_includes_frozen_ratio_uncertainty(self) -> None:
        without = scalar_score(10.0, 13.0, 1.0, 4.0, 1.2)
        with_source = scalar_score(10.0, 13.0, 1.0, 4.0, 1.2, ratio_se=0.1)
        self.assertGreater(with_source["variance"], without["variance"])
        self.assertLess(abs(with_source["signed_z"]), abs(without["signed_z"]))

    def test_primary_score_order_remains_frozen(self) -> None:
        plan = yaml.safe_load(
            (ROOT / "experiments/p50_n145_n290_fullcurve_20260829.yaml").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            plan["frozen_scoring"]["order"],
            [
                "intrinsic_and_thermal_even_DeltaM_transfer",
                "raw_asymptotic_slope_baseline",
                "frozen_scalar_plus_H4_slope_correction",
                "raw_root_ratio_baseline",
                "frozen_induced_root_ratio",
            ],
        )


if __name__ == "__main__":
    unittest.main()
