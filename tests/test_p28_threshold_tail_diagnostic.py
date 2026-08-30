import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from diagnose_p28_threshold_profile_tail_rejection import (  # noqa: E402
    canonical_modes,
    local_effective_exponents,
)


class P28TailDiagnosticTests(unittest.TestCase):
    def test_modes_annihilate_frozen_model(self):
        z = [2.5, 2.75, 3.0, 3.25, 3.5]
        modes = canonical_modes(z, "stretched_4_over_3")
        self.assertEqual(len(modes), 3)
        for row in modes:
            self.assertAlmostEqual(sum(row), 0.0, places=11)
            self.assertAlmostEqual(sum(a * b ** (4 / 3) for a, b in zip(row, z)), 0.0, places=11)

    def test_local_exponent_recovers_power(self):
        z = [2.5, 2.75, 3.0, 3.25, 3.5]
        beta = 4 / 3
        # The secant-based diagnostic converges to beta; at this grid it is close.
        values = [1.7 - 0.8 * value**beta for value in z]
        observed = local_effective_exponents(values, z)
        for value in observed:
            self.assertAlmostEqual(value, beta, delta=0.01)

    def test_nested_modes_have_two_dimensions(self):
        z = [2.5, 2.75, 3.0, 3.25, 3.5]
        self.assertEqual(len(canonical_modes(z, "post_reveal_4_over_3_plus_2_over_3")), 2)


if __name__ == "__main__":
    unittest.main()
