from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_boolean_noise_semigroup import render  # noqa: E402


class ExactBooleanNoiseSemigroupTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = render()
        cls.artifact = json.loads(
            (
                ROOT
                / "results"
                / "exact-boolean-noise-semigroup"
                / "oracle.json"
            ).read_text(encoding="utf-8")
        )

    def test_direct_noisy_pairs_equal_fourier_generating_functions(self) -> None:
        for geometry in self.payload["geometries"]:
            checks = geometry["direct_noisy_pair_cross_checks"]
            self.assertTrue(checks["all_exact"])
            self.assertEqual(checks["rho_grid"], ["0", "1/4", "1/2", "3/4", "1"])
            for rows in checks["pairs"].values():
                for row in rows:
                    self.assertEqual(row["difference"], "0")
                    self.assertEqual(
                        row["direct_noisy_pair_covariance"],
                        row["fourier_generating_value"],
                    )

    def test_autocorrelation_coefficients_are_nonnegative(self) -> None:
        for geometry in self.payload["geometries"]:
            for spectrum in geometry["autocorrelation_generating_functions"].values():
                self.assertTrue(
                    all(Fraction(value) >= 0 for value in spectrum["coefficients"])
                )

    def test_rho_one_derivative_is_p_times_q_pivotal_mass(self) -> None:
        pq = Fraction(2, 5) * Fraction(3, 5)
        for geometry in self.payload["geometries"]:
            for row in geometry["rho_one_pivotal_checks"].values():
                derivative = Fraction(row["raw_autocorrelation_derivative"])
                mass = Fraction(row["total_unsigned_pivotal_mass"])
                self.assertEqual(derivative, pq * mass)
                self.assertEqual(
                    Fraction(row["derivative_divided_by_p_times_q"]), mass
                )

    def test_orientation_and_matching_cross_spectra_are_signed(self) -> None:
        for geometry in self.payload["geometries"]:
            spectra = geometry["cross_spectra"]
            for name in ("orientation_cross_spectrum", "matching_cross_spectrum"):
                coefficients = [Fraction(value) for value in spectra[name]["coefficients"]]
                self.assertTrue(any(value > 0 for value in coefficients))
                self.assertTrue(any(value < 0 for value in coefficients))
            self.assertTrue(
                all(
                    Fraction(value) == 0
                    for value in spectra["H4_matching_odd_cross_spectrum"]["coefficients"]
                )
            )

    def test_n10_level_support_is_frozen(self) -> None:
        n10 = next(row for row in self.payload["geometries"] if row["N"] == 10)
        orientation = [
            Fraction(value)
            for value in n10["autocorrelation_generating_functions"][
                "orientation_difference"
            ]["coefficients"]
        ]
        matching_odd = [
            Fraction(value)
            for value in n10["autocorrelation_generating_functions"][
                "matching_odd_cross"
            ]["coefficients"]
        ]
        self.assertEqual(
            {index for index, value in enumerate(orientation) if value}, {2, 3, 4}
        )
        self.assertEqual(
            {index for index, value in enumerate(matching_odd) if value},
            {1, 2, 3, 4, 5},
        )

    def test_frozen_artifact_matches_oracle(self) -> None:
        self.assertEqual(self.artifact, self.payload)


if __name__ == "__main__":
    unittest.main()
