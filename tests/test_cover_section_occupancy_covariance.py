from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import tempfile
import unittest

from scripts.cover_section_occupancy_covariance import (
    antithetic_covariance,
    build_contract,
    coupling_case,
    section_covariance,
    squared_correlation,
    validate_contract,
)


class CoverSectionOccupancyCovarianceTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_section_squared_correlation_is_one_over_degree(self) -> None:
        for degree in (1, 2, 5, 11):
            for probability in (Fraction(1, 7), Fraction(2, 5), Fraction(3, 5)):
                covariance = section_covariance(degree, probability)
                self.assertEqual(
                    squared_correlation(covariance, degree, probability),
                    Fraction(1, degree),
                )

    def test_antithetic_formula_is_symmetric_and_negative(self) -> None:
        for degree in (2, 5):
            left = antithetic_covariance(degree, Fraction(2, 5))
            right = antithetic_covariance(degree, Fraction(3, 5))
            self.assertEqual(left, right)
            self.assertLess(left, 0)

    def test_frozen_norm_two_and_five_values_are_exact(self) -> None:
        self.assertEqual(coupling_case(2, Fraction(2, 5))["section_covariance"], "3/25")
        self.assertEqual(coupling_case(2, Fraction(2, 5))["antithetic_squared_correlation"], "2/9")
        self.assertEqual(coupling_case(5, Fraction(3, 5))["section_covariance"], "6/125")
        self.assertEqual(coupling_case(5, Fraction(3, 5))["antithetic_squared_correlation"], "4/45")

    def test_contract_drift_and_invalid_inputs_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "strictly between"):
            section_covariance(2, Fraction())
        frozen = build_contract()
        frozen["frozen_cases"][0]["section_covariance"] = "0"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "drift.json"
            path.write_text(json.dumps(frozen), encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "contract drifted"):
                validate_contract(path)


if __name__ == "__main__":
    unittest.main()
