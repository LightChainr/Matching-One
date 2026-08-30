from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
import tempfile
import unittest

from scripts.gaussian_mixed_curvature_oracle import (
    build_contract,
    mixed_log_curvature,
    ordinary_power_mixed_response,
    validate_contract,
)


class GaussianMixedCurvatureOracleTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_constant_and_affine_log_laws_are_annihilated(self) -> None:
        for constant, linear in ((0, 0), (7, 0), (0, 11), (-5, 13)):
            self.assertEqual(
                mixed_log_curvature(Fraction(constant), Fraction(linear), Fraction(0)),
                {},
            )

    def test_quadratic_log_keeps_only_the_cross_curvature(self) -> None:
        for constant, linear, quadratic in ((0, 0, 1), (7, 11, 13), (-3, 5, -2)):
            self.assertEqual(
                mixed_log_curvature(
                    Fraction(constant), Fraction(linear), Fraction(quadratic)
                ),
                {"log2_log5": Fraction(2 * quadratic)},
            )

    def test_ordinary_power_response_factorizes_exactly(self) -> None:
        expected = {1: Fraction(2, 5), 2: Fraction(18, 25), 3: Fraction(217, 250)}
        for beta, target in expected.items():
            self.assertEqual(ordinary_power_mixed_response(beta), target)
            self.assertGreater(target, 0)

    def test_contract_drift_fails_closed(self) -> None:
        frozen = build_contract()
        frozen["ordinary_power_responses"]["2"] = "17/25"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "drift.json"
            path.write_text(json.dumps(frozen), encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "contract drifted"):
                validate_contract(path)


if __name__ == "__main__":
    unittest.main()
