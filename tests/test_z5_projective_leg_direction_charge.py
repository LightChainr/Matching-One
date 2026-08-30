from __future__ import annotations

import cmath
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_projective_leg_direction_charge import (  # noqa: E402
    ZETA5,
    direction_factor,
    wrap_phase,
)


class Z5ProjectiveLegDirectionChargeTests(unittest.TestCase):
    def test_internal_direction_factor_is_zero(self) -> None:
        self.assertEqual(direction_factor(1.0 + 0j), 0j)

    def test_minus_i_direction_factor_is_plus_i(self) -> None:
        self.assertAlmostEqual(abs(direction_factor(-1j) - 1j), 0.0, places=14)

    def test_c4_minus_i_is_not_in_z5_alphabet(self) -> None:
        distances = [abs(wrap_phase(-math.pi / 2.0 - cmath.phase(ZETA5**index))) for index in range(5)]
        self.assertGreater(min(distances), 0.3)


if __name__ == "__main__":
    unittest.main()
