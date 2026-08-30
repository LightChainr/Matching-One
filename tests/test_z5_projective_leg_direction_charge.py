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
    fit_difference_amplitudes,
    unrestricted_shared_roots_heldout_residual,
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

    def test_difference_amplitudes_recover_two_modes(self) -> None:
        roots = (0.6 + 0.1j, -0.2 - 0.3j)
        amplitudes = {
            ("plus", 1): (1.0 + 0.2j, -0.3 + 0.1j),
            ("plus", 2): (0.5 - 0.4j, 0.2 + 0.3j),
            ("minus", 1): (-0.7 + 0.1j, 0.4 - 0.2j),
            ("minus", 2): (0.9 + 0.5j, -0.1 - 0.6j),
        }
        values = {}
        for channel, row in amplitudes.items():
            hand, charge = channel
            for distance in range(1, 5):
                value = sum(row[index] * roots[index] ** (distance - 1) for index in range(2))
                values[f"d{distance}_T{charge}_{hand}_re"] = 0.0
                values[f"d{distance}_T{charge}_{hand}_im"] = 0.0
                values[f"d{distance}_A{charge}_{hand}_re"] = value.real
                values[f"d{distance}_A{charge}_{hand}_im"] = value.imag
        fitted = fit_difference_amplitudes(values, roots)
        for channel, expected in amplitudes.items():
            for observed, target in zip(fitted[channel], expected):
                self.assertAlmostEqual(abs(observed - target), 0.0, places=12)

    def test_unrestricted_shared_roots_predict_fifth_distance(self) -> None:
        # Supply both T and A from the same exact rank-two recurrence.  The T
        # row determines the roots while A has independent amplitudes.
        roots = (0.55 + 0.05j, -0.1 - 0.35j)
        values = {}
        for channel_index, (hand, charge) in enumerate((("plus", 1), ("plus", 2), ("minus", 1), ("minus", 2))):
            t_amplitudes = (1 + 0.1j * channel_index, 0.2 - 0.05j * channel_index)
            a_amplitudes = (0.3 + 0.02j * channel_index, -0.1 + 0.03j * channel_index)
            for distance in range(1, 6):
                tvalue = sum(t_amplitudes[index] * roots[index] ** (distance - 1) for index in range(2))
                avalue = sum(a_amplitudes[index] * roots[index] ** (distance - 1) for index in range(2))
                values[f"d{distance}_T{charge}_{hand}_re"] = tvalue.real
                values[f"d{distance}_T{charge}_{hand}_im"] = tvalue.imag
                values[f"d{distance}_A{charge}_{hand}_re"] = avalue.real
                values[f"d{distance}_A{charge}_{hand}_im"] = avalue.imag
        residual = unrestricted_shared_roots_heldout_residual(values)
        self.assertLess(max(abs(value) for value in residual), 1e-12)


if __name__ == "__main__":
    unittest.main()
