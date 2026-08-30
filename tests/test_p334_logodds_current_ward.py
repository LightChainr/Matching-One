from __future__ import annotations

from fractions import Fraction
import sys
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p334_logodds_current_ward import build_certificate  # noqa: E402


class ProjectiveCurrentWardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_certificate()

    def test_both_exact_quotients_close(self) -> None:
        self.assertTrue(self.payload["all_exact_gates_pass"])
        self.assertEqual([row["N"] for row in self.payload["quotients"]], [13, 17])

    def test_each_orbit_closes_polynomial_and_eta_jet(self) -> None:
        for quotient in self.payload["quotients"]:
            for row in quotient["orbit_rows"]:
                self.assertTrue(row["ward_polynomial_exact"])
                self.assertEqual(row["eta_jet_shift_exact_through_order"], 6)

    def test_shares_are_coordinate_free_and_complete(self) -> None:
        for quotient in self.payload["quotients"]:
            shares = quotient["coordinate_free_orbit_shares"]
            self.assertTrue(all(row["same_in_p_and_eta"] for row in shares))
            self.assertEqual(sum(Fraction(row["exact_share"]) for row in shares), 1)

    def test_current_sum_rule(self) -> None:
        for quotient in self.payload["quotients"]:
            self.assertTrue(quotient["empty_and_full_state_zero"])
            self.assertTrue(quotient["integrated_net_current_zero"])

    def test_each_orbit_and_total_has_one_exact_stationary_point(self) -> None:
        for quotient in self.payload["quotients"]:
            certificate = quotient["unique_stationary_point_certificate"]
            self.assertTrue(certificate["all_unique"])
            self.assertEqual(certificate["total_net_bernstein_sign_variations"], 1)
            self.assertTrue(all(certificate["orbit_and_total"].values()))


if __name__ == "__main__":
    unittest.main()
