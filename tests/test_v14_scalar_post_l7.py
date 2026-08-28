from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from v14_scalar_post_l7 import (  # noqa: E402
    conditional_interchiral_parity,
    critical_potts_h,
    diagonal_x,
    mechanism_table,
)


class V14ScalarPostL7Tests(unittest.TestCase):
    def test_v14_spectrum_and_annihilator_power(self) -> None:
        self.assertEqual(critical_potts_h(1, 4), F(33, 8))
        self.assertEqual(diagonal_x(4), F(33, 4))
        v14 = mechanism_table()["V14_scalar_H0"]
        self.assertEqual(v14["residual_power_L"], F(25, 4))
        self.assertEqual(v14["residual_power_N"], F(25, 8))
        self.assertEqual(v14["relative_q_L"], F(3))
        self.assertEqual(v14["root_power_L"], F(7))

    def test_next_thermal_h4_is_q6_w10(self) -> None:
        nxt = mechanism_table()["thermal_next_H4"]
        self.assertEqual(nxt["x"], F(45, 4))
        self.assertEqual(nxt["residual_power_L"], F(37, 4))
        self.assertEqual(nxt["relative_q_L"], F(6))
        self.assertEqual(nxt["root_power_L"], F(10))

    def test_parity_helper_is_explicitly_conditional(self) -> None:
        self.assertEqual(conditional_interchiral_parity(2), -1)
        self.assertEqual(conditional_interchiral_parity(3), 1)
        self.assertEqual(conditional_interchiral_parity(4), -1)


if __name__ == "__main__":
    unittest.main()
