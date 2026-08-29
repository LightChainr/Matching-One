#!/usr/bin/env python3

from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from q1_spin4_velocity_oracle import (  # noqa: E402
    REQUIRED_DERIVATIVE_SEMANTICS,
    exact_targets,
    score_input,
    synthetic_oracle,
    velocity_estimate,
)


class Q1Spin4VelocityOracleTests(unittest.TestCase):
    def test_exact_targets_and_gap(self) -> None:
        targets = exact_targets()
        fields = targets["fields"]
        self.assertEqual(
            fields["four_leg_V_2_2"]["dx_dQ_at_Q1"]["coefficient_sqrt3_over_pi"],
            "-5/16",
        )
        self.assertEqual(
            fields["thermal_Q4_epsilon"]["dx_dQ_at_Q1"]["coefficient_sqrt3_over_pi"],
            "-9/16",
        )
        self.assertEqual(
            targets["velocity_gap_four_leg_minus_thermal_Q4"]["coefficient_sqrt3_over_pi"],
            "1/4",
        )

    def test_synthetic_oracles_cancel_amplitude_velocity(self) -> None:
        payload = synthetic_oracle(exact_targets())
        for row in payload.values():
            self.assertLess(row["absolute_error"], 1e-14)
            self.assertEqual(row["normalization_velocity_cancelled"], -0.43)

    def test_full_covariance_changes_velocity_variance(self) -> None:
        point = [0.2, 0.03, 0.1, 0.01]
        diagonal = [[0.0] * 4 for _ in range(4)]
        for i, value in enumerate((1e-4, 2e-4, 1.5e-4, 2.5e-4)):
            diagonal[i][i] = value
        correlated = [row[:] for row in diagonal]
        correlated[1][3] = correlated[3][1] = 5e-5
        first = velocity_estimate([8, 16], point, diagonal)
        second = velocity_estimate([8, 16], point, correlated)
        self.assertNotEqual(first["variance"], second["variance"])

    def test_measure_only_input_is_not_scoreable(self) -> None:
        payload = {
            "schema": "matching-one.q-velocity-two-size-input.v1",
            "sizes": [8, 16],
            "point": [1.0, 0.0, 0.5, 0.0],
            "covariance": [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)],
            "derivative_semantics": "measure_score_only",
            "explicit_field_definition_derivative_included": False,
        }
        result = score_input(payload, exact_targets())
        self.assertEqual(result["status"], "NOT_SCOREABLE")
        payload["derivative_semantics"] = REQUIRED_DERIVATIVE_SEMANTICS
        result = score_input(payload, exact_targets())
        self.assertEqual(result["status"], "NOT_SCOREABLE")


if __name__ == "__main__":
    unittest.main()
