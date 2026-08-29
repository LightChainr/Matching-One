import math
import sys
import unittest
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p50_matrix_state_transfer import (  # noqa: E402
    STATE_ORDER,
    score_model,
    state_from_statistics,
)


class P50MatrixStateTransferTests(unittest.TestCase):
    def test_state_definition(self) -> None:
        stat = {
            "mean_slope": 2.0,
            "P4_S": 3.0,
            "P4_D_prime": 5.0,
            "P4_D": 7.0,
            "P4_S_prime": 11.0,
        }
        state = state_from_statistics(stat, 4)
        self.assertEqual(state[0], 12.0)
        self.assertEqual(state[1], 10.0)
        self.assertAlmostEqual(state[2], (4.0 ** (13.0 / 8.0)) * 7.0)
        self.assertAlmostEqual(state[3], (4.0 ** (13.0 / 8.0)) * 11.0 / 2.0)

    def test_frozen_increment_arithmetic(self) -> None:
        payload = yaml.safe_load(
            (ROOT / "predictions" / "p50_matrix_state_transfer_20260829.yaml").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(tuple(payload["state_definition"]["state_order"]), STATE_ORDER)
        ordinary = payload["models"]["analytic_inverse_N"]
        jordan = payload["models"]["rank2_Jordan_log"]

        c = float(ordinary["source_correction_coefficient_C"])
        expected_ordinary = c * (1.0 / 290.0 - 1.0 / 145.0)
        self.assertAlmostEqual(
            float(ordinary["expected_child_minus_parent"]["T_Su"]), expected_ordinary, places=14
        )
        self.assertAlmostEqual(
            float(ordinary["source_se_T_Su_increment"]),
            math.sqrt(float(ordinary["source_variance_C"])) / 290.0,
            places=14,
        )

        b = float(jordan["source_log_coefficient_B"])
        expected_jordan = b * math.log(2.0)
        self.assertAlmostEqual(
            float(jordan["expected_child_minus_parent"]["T_Su"]), expected_jordan, places=14
        )
        self.assertAlmostEqual(
            float(jordan["source_se_T_Su_increment"]),
            math.sqrt(float(jordan["source_variance_B"])) * math.log(2.0),
            places=14,
        )

    def test_zero_residual_scores_zero(self) -> None:
        payload = yaml.safe_load(
            (ROOT / "predictions" / "p50_matrix_state_transfer_20260829.yaml").read_text(
                encoding="utf-8"
            )
        )
        model = payload["models"]["rank2_Jordan_log"]
        expected = [float(model["expected_child_minus_parent"][name]) for name in STATE_ORDER]
        covariance = [
            [1.0 if i == j else 0.0 for j in range(4)]
            for i in range(4)
        ]
        score = score_model(expected, covariance, model)
        self.assertAlmostEqual(float(score["joint_chi_square"]), 0.0, places=14)
        self.assertAlmostEqual(float(score["marginal_T_Su_signed_z"]), 0.0, places=14)
        self.assertAlmostEqual(float(score["chi_square_survival_df4"]), 1.0, places=14)


if __name__ == "__main__":
    unittest.main()
