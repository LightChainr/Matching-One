from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_matching_polynomial import bernstein_to_power  # noqa: E402
from matching_defect_polynomial import (  # noqa: E402
    AXIS_BERNSTEIN,
    N10_BERNSTEIN,
    N10_POWER,
    run_suite,
)


class MatchingDefectPolynomialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = run_suite()

    def test_suite_passes_and_rejects_tutte(self) -> None:
        self.assertTrue(self.payload["passed"])
        self.assertFalse(self.payload["tutte_specialization"])
        self.assertEqual(
            self.payload["identification"],
            "M(p) is the Bernoulli generating function of the wrapping-event q",
        )

    def test_axis_bernstein_reproduces_committed_matching_polynomials(self) -> None:
        for row in self.payload["axis"]:
            self.assertEqual(row["bernstein"], AXIS_BERNSTEIN[row["L"]], row["name"])
            self.assertEqual(row["identity_failures"], 0, row["name"])
            self.assertFalse(row["obstruction"]["tutte_specialization"])
            self.assertFalse(row["obstruction"]["primal_equals_matching"])
            self.assertEqual(row["power"], bernstein_to_power(row["bernstein"]))

    def test_axis_l2_power_basis_is_the_committed_biquadratic(self) -> None:
        row = next(item for item in self.payload["axis"] if item["L"] == 2)
        self.assertEqual(row["power"], [-1, 0, 4, 0, -2])
        # Not antisymmetric under k -> N-k, so complement is not q -> -q.
        self.assertNotEqual(row["complement_involution_failures"], 0)

    def test_n10_is_beta33_and_complement_is_the_duality_involution(self) -> None:
        n10 = self.payload["n10"]
        self.assertEqual(n10["bernstein"], N10_BERNSTEIN)
        self.assertEqual(n10["power"], N10_POWER)
        self.assertEqual(n10["complement_involution_failures"], 0)
        self.assertTrue(n10["obstruction"]["primal_equals_matching"])
        self.assertFalse(n10["obstruction"]["tutte_specialization"])
        self.assertTrue(n10["beta33"])


if __name__ == "__main__":
    unittest.main()
