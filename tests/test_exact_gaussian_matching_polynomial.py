from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_gaussian_matching_polynomial import enumerate_gaussian  # noqa: E402


CASES = [
    (3, 1, 10, "gaussian_3_1_target.json"),
    (3, 2, 13, "gaussian_3_2_target.json"),
    (4, 1, 17, "gaussian_4_1_target.json"),
]


class ExactGaussianMatchingPolynomialTests(unittest.TestCase):
    def test_independent_frozen_vectors_match_canonical_topology_engine(self) -> None:
        for a, b, n, filename in CASES:
            with self.subTest(a=a, b=b, N=n):
                expected = json.loads(
                    (
                        ROOT
                        / "results"
                        / "exact-axis-l5-frontier"
                        / filename
                    ).read_text(encoding="utf-8")
                )
                observed = enumerate_gaussian(a, b, "either", max_n=n)
                self.assertEqual(observed["geometry"]["N"], n)
                self.assertEqual(
                    observed["primal_wrap_counts_by_occupancy"],
                    expected["primal_wrap_counts_by_occupancy"],
                )
                self.assertEqual(
                    observed["matching_complement_wrap_counts_by_occupancy"],
                    expected["matching_complement_wrap_counts_by_occupancy"],
                )
                self.assertEqual(
                    observed["bernstein_sums"], expected["bernstein_sums"]
                )
                self.assertEqual(
                    observed["power_coefficients_ascending"],
                    expected["power_coefficients_ascending"],
                )
                self.assertEqual(observed["degree"], n)


if __name__ == "__main__":
    unittest.main()
