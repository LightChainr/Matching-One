import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from hexagonal_pell_spin_filter import pell_plus  # noqa: E402


class HexagonalPellPeriodMatrixTests(unittest.TestCase):
    def test_shape_error_is_quadratically_small(self) -> None:
        errors = []
        for x, m in pell_plus(4):
            error = abs(x / (2 * m) - math.sqrt(3) / 2)
            errors.append((m, error))
            self.assertEqual(x * x - 3 * m * m, 1)
        scaled = [m * m * error for m, error in errors]
        self.assertLess(max(scaled[-2:]) - min(scaled[-2:]), 5e-5)
        self.assertAlmostEqual(scaled[-1], 1 / (4 * math.sqrt(3)), places=5)


if __name__ == "__main__":
    unittest.main()
