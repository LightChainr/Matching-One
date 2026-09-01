import math
import sys
import unittest
from pathlib import Path

import mpmath as mp
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from hexagonal_pell_spin_filter import build, eisenstein_e4, pell_plus  # noqa: E402


class HexagonalPellSpinFilterTests(unittest.TestCase):
    def test_first_pell_targets_and_determinants(self) -> None:
        rows = pell_plus(3)
        self.assertEqual(rows, [(7, 4), (26, 15), (97, 56)])
        self.assertEqual([2 * x * m for x, m in rows], [56, 780, 10864])
        for x, m in rows:
            self.assertEqual(x * x - 3 * m * m, 1)

    def test_hexagonal_E4_zero_and_pell_ratios(self) -> None:
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

    def test_pell_shape_error_is_inverse_square(self) -> None:
        scaled = []
        for x, m in pell_plus(4):
            error = abs(x / (2 * m) - math.sqrt(3) / 2)
            scaled.append(m * m * error)
        self.assertAlmostEqual(scaled[-1], 1 / (4 * math.sqrt(3)), places=5)

    def test_frozen_prediction_artifact(self) -> None:
        path = ROOT / "predictions" / "hexagonal_pell_spin_filter_20260828.yaml"
        payload = yaml.safe_load(path.read_text(encoding="utf-8"))
        h4 = payload["frozen_hypotheses"]["H4_thermal_level4"]
        h12 = payload["frozen_hypotheses"]["H12_same_radial_alias"]
        self.assertEqual(h4["pell_root_bias_length_exponent"], "6")
        self.assertEqual(h12["pell_root_bias_length_exponent"], "4")
        self.assertEqual(h4["pell_matching_residual_length_exponent"], "21/4")
        self.assertEqual(h12["pell_matching_residual_length_exponent"], "13/4")


if __name__ == "__main__":
    unittest.main()
