from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_gaussian_matching_polynomial import enumerate_gaussian  # noqa: E402


class ExactGaussianMatchingPolynomialTests(unittest.TestCase):
    def test_gaussian_3_1_matches_independent_frozen_vector(self) -> None:
        expected = json.loads(
            (
                ROOT
                / "results"
                / "exact-axis-l5-frontier"
                / "gaussian_3_1_target.json"
            ).read_text(encoding="utf-8")
        )
        observed = enumerate_gaussian(3, 1, "either", max_n=10)
        self.assertEqual(observed["geometry"]["N"], 10)
        self.assertEqual(
            observed["primal_wrap_counts_by_occupancy"],
            expected["primal_wrap_counts_by_occupancy"],
        )
        self.assertEqual(
            observed["matching_complement_wrap_counts_by_occupancy"],
            expected["matching_complement_wrap_counts_by_occupancy"],
        )
        self.assertEqual(observed["bernstein_sums"], expected["bernstein_sums"])
        self.assertEqual(
            observed["power_coefficients_ascending"],
            expected["power_coefficients_ascending"],
        )
        self.assertEqual(observed["degree"], 10)


if __name__ == "__main__":
    unittest.main()
