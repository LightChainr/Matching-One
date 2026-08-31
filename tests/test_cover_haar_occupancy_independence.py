from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest

from scripts.cover_haar_occupancy_independence import (
    bounded_box_cdf,
    build_contract,
    fractional_window_volume,
    joint_parent_child,
    occupancy_covariance,
    validate_contract,
)


class CoverHaarOccupancyIndependenceTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_h0_and_h1_joint_laws_factor_exactly(self) -> None:
        for degree in range(2, 7):
            for probability in (Fraction(1, 5), Fraction(2, 5), Fraction(1, 2), Fraction(5, 7)):
                product = probability * probability
                self.assertEqual(joint_parent_child(degree, probability), product)
                self.assertEqual(
                    joint_parent_child(degree, probability, antithetic=True), product
                )

    def test_covariance_with_the_fiber_mean_is_zero(self) -> None:
        for degree in (2, 3, 5):
            for probability in (Fraction(1, 3), Fraction(2, 5), Fraction(3, 5)):
                self.assertEqual(occupancy_covariance(degree, probability), 0)
                self.assertEqual(
                    occupancy_covariance(degree, probability, antithetic=True), 0
                )

    def test_rational_volume_oracle_handles_box_and_fractional_windows(self) -> None:
        self.assertEqual(
            bounded_box_cdf(Fraction(2), [Fraction(1), Fraction(1)]), 1
        )
        self.assertEqual(
            fractional_window_volume([Fraction(2, 5), Fraction(1)], Fraction(3, 5)),
            Fraction(6, 25),
        )
        with self.assertRaisesRegex(ValueError, "degree at least two"):
            joint_parent_child(1, Fraction(1, 2))

    def test_contract_drift_fails_closed(self) -> None:
        frozen = build_contract()
        frozen["frozen_cases"][0]["h0_covariance_with_fiber_mean"] = "1"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "drift.json"
            path.write_text(json.dumps(frozen), encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "contract drifted"):
                validate_contract(path)


if __name__ == "__main__":
    unittest.main()
