from __future__ import annotations

import cmath
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_projective_leg_pair_transfer import CHANNELS  # noqa: E402
from score_z5_projective_leg_state_dimension import (  # noqa: E402
    coefficients_from_roots,
    fit_recurrence,
    image_kernel,
    recurrence_residual,
    recurrence_roots,
    transport_coefficients,
)


class Z5ProjectiveLegStateDimensionTests(unittest.TestCase):
    def synthetic_rank2(self):
        roots = [0.72 * cmath.exp(0.11j), 0.39 * cmath.exp(-0.21j)]
        rows = {}
        for index, channel in enumerate(CHANNELS):
            first = complex(1.0 + index, 0.2 * index)
            second = complex(-0.3 * index, 0.7 + index)
            rows[channel] = [first * roots[0] ** d + second * roots[1] ** d for d in range(5)]
        return roots, rows

    def test_rank2_recurrence_recovers_two_states(self) -> None:
        roots, rows = self.synthetic_rank2()
        expected = coefficients_from_roots(roots)
        fitted = fit_recurrence(rows, 2, 4)
        for observed, target in zip(fitted, expected):
            self.assertAlmostEqual(abs(observed - target), 0.0, places=12)
        self.assertLess(max(abs(value) for value in recurrence_residual(rows, fitted, 5)), 1e-12)

    def test_rank1_does_not_close_rank2_sequence(self) -> None:
        _, rows = self.synthetic_rank2()
        fitted = fit_recurrence(rows, 1, 4)
        self.assertGreater(max(abs(value) for value in recurrence_residual(rows, fitted, 5)), 1e-3)

    def test_conformal_transport_scales_complex_logs(self) -> None:
        roots = [0.7 * cmath.exp(0.1j), 0.4 * cmath.exp(-0.2j)]
        coefficients = coefficients_from_roots(roots)
        transported = recurrence_roots(transport_coefficients(coefficients, 1.5))
        expected = [cmath.exp(1.5 * cmath.log(root)) for root in roots]
        self.assertAlmostEqual(sorted(abs(value) for value in transported)[0], sorted(abs(value) for value in expected)[0], places=12)

    def test_nearest_image_kernel_is_positive(self) -> None:
        self.assertGreater(image_kernel(5, (10, 1), 1.25), 0.0)


if __name__ == "__main__":
    unittest.main()
