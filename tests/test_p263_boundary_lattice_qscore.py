from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p263_boundary_lattice_qscore import (  # noqa: E402
    LINK_PATTERNS,
    direct_polynomial_tangent,
    measure_tangent_from_sufficient,
    render,
    square_cycle_rows,
    sufficient_statistics,
)


class BoundaryLatticeQScoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = render()
        cls.frozen = json.loads(
            (ROOT / "predictions" / "p263_boundary_lattice_qscore_20260829.json").read_text(
                encoding="utf-8"
            )
        )

    def test_tiny_link_probabilities_and_tangents(self) -> None:
        derivatives = self.payload["tiny_square_bond_regression"]["derivatives"]
        self.assertEqual(derivatives["1234"]["probability_Q1"], "5/16")
        self.assertEqual(derivatives["12|34"]["probability_Q1"], "1/16")
        self.assertEqual(derivatives["14|23"]["probability_Q1"], "1/16")
        self.assertEqual(derivatives["1234"]["measure_score_covariance"], "-37/256")
        self.assertEqual(derivatives["12|34"]["measure_score_covariance"], "-1/256")
        self.assertEqual(derivatives["14|23"]["measure_score_covariance"], "-1/256")

    def test_measure_score_matches_direct_ratio_differentiation(self) -> None:
        rows = square_cycle_rows()
        statistics = sufficient_statistics(rows)
        for pattern in LINK_PATTERNS:
            self.assertEqual(
                measure_tangent_from_sufficient(statistics, pattern),
                direct_polynomial_tangent(rows, pattern),
            )

    def test_frozen_geometries_have_declared_cross_ratios(self) -> None:
        for geometry in self.payload["frozen_geometries"]:
            x1, x2, x3, x4 = [
                float(Fraction(value)) for value in geometry["normalized_boundary_points"]
            ]
            lam = (x2 - x1) * (x4 - x3) / ((x4 - x2) * (x3 - x1))
            self.assertAlmostEqual(lam, float(Fraction(geometry["lambda"])), places=15)

    def test_amplitude_projected_high_branch_target_is_frozen(self) -> None:
        target = self.payload["amplitude_gauge_and_score"]["frozen_target"]
        self.assertEqual(target["lambda_order"], ["1/4", "1/3", "2/3", "3/4"])
        self.assertEqual(target["anchor_lambda"], "1/3")
        self.assertAlmostEqual(target["anchored_dQ_logU"][1], 0.0, places=15)
        self.assertLess(target["max_abs_order_difference"], 1e-12)
        self.assertAlmostEqual(target["anchored_dQ_logU"][0], -0.14829994193484564, places=13)
        self.assertAlmostEqual(target["anchored_dQ_logU"][2], 0.3249693720763415, places=13)
        self.assertAlmostEqual(target["anchored_dQ_logU"][3], 0.36786256566279774, places=13)

    def test_derivative_ledger_keeps_terms_separate(self) -> None:
        ledger = self.payload["q_score_decomposition"]
        self.assertIn("Cov(O,J/2)", ledger["measure_score"])
        self.assertIn("zero", ledger["explicit_projector_derivative"])
        self.assertIn("4*h_prime*log(L)", ledger["explicit_field_derivative"])
        self.assertIn("2*h_prime*log(K)", ledger["conformal_prefactor_derivative"])

    def test_frozen_artifact_matches(self) -> None:
        self.assertEqual(self.payload, self.frozen)


if __name__ == "__main__":
    unittest.main()
