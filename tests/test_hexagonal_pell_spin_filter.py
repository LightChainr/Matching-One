import sys
import unittest
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from hexagonal_pell_spin_filter import build, eisenstein_e4, pell_plus  # noqa: E402


class HexagonalPellSpinFilterTests(unittest.TestCase):
    def test_first_nontrivial_pell_solutions(self) -> None:
        rows = pell_plus(3)
        self.assertEqual(rows, [(7, 4), (26, 15), (97, 56)])
        for x, m in rows:
            self.assertEqual(x * x - 3 * m * m, 1)

    def test_period_determinants(self) -> None:
        payload = build(3, 60, 90)
        self.assertEqual(
            [row["site_count_det"] for row in payload["records"]],
            [56, 780, 10864],
        )

    def test_E4_hexagonal_zero_and_frozen_ratios(self) -> None:
        mp.mp.dps = 60
        rho = mp.mpf(1) / 2 + 1j * mp.sqrt(3) / 2
        self.assertLess(abs(eisenstein_e4(rho, 100)), mp.mpf("1e-45"))
        payload = build(3, 60, 100)
        ratios = [mp.mpf(row["E4_over_E4_i"]) for row in payload["records"]]
        expected = [
            mp.mpf("0.0362565052970333734263688162329262"),
            mp.mpf("0.00265353663269522063293766164637801"),
            mp.mpf("0.000190777377235422695593606819608911"),
        ]
        for actual, target in zip(ratios, expected):
            self.assertLess(abs(actual - target), mp.mpf("1e-34"))

    def test_Pell_E4_ratio_has_m_inverse_two_limit(self) -> None:
        payload = build(4, 60, 100)
        scaled = [mp.mpf(row["m2_times_E4_ratio"]) for row in payload["records"]]
        self.assertLess(abs(scaled[-1] - scaled[-2]), mp.mpf("2e-4"))
        self.assertGreater(scaled[-1], mp.mpf("0.598"))
        self.assertLess(scaled[-1], mp.mpf("0.599"))


if __name__ == "__main__":
    unittest.main()
