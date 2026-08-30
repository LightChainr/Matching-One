from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p321_ttbar_shape_operator import (  # noqa: E402
    build_artifact,
    first_order_mode_multiplier,
    normalized_deformed_mode,
    product_rule_obstruction,
    undeformed_mode,
)


class P321TtbarShapeOperatorTests(unittest.TestCase):
    def test_exact_small_alpha_mode_coefficient(self) -> None:
        k, x, p, y, theta = 1.25, 0.75, 2.0, 1.3, 0.17
        base = undeformed_mode(x, p, y, theta)
        expected = first_order_mode_multiplier(k, x, p, y)
        estimates = []
        for alpha in (1.0e-5, 5.0e-6):
            deformed = normalized_deformed_mode(k, x, p, y, theta, alpha)
            estimates.append(((deformed / base - 1.0) / alpha).real)
        richardson = 2.0 * estimates[-1] - estimates[0]
        self.assertLess(abs(richardson - expected), 2.0e-3)
        self.assertLess(abs(estimates[-1] - expected), abs(estimates[0] - expected))

    def test_zero_frequency_is_unchanged(self) -> None:
        self.assertEqual(first_order_mode_multiplier(3.0, 0.0, 0.0, 2.0), 0.0)
        self.assertAlmostEqual(
            normalized_deformed_mode(3.0, 0.0, 0.0, 2.0, 0.0, 0.01).real,
            1.0,
        )

    def test_angular_frequency_enters_with_opposite_sign(self) -> None:
        radial = first_order_mode_multiplier(0.0, 2.0, 0.0, 1.0)
        angular = first_order_mode_multiplier(0.0, 2.0, 3.0, 1.0)
        self.assertAlmostEqual(radial - angular, 4.0 * math.pi**2 * 9.0)

    def test_product_identity_names_thermal_derivative(self) -> None:
        identity = product_rule_obstruction()
        self.assertIn("grad(log f)", identity["ratio_identity"])
        self.assertIn("F_t", identity["p321_substitution"]["f"])

    def test_artifact_keeps_spin_boundary(self) -> None:
        artifact = build_artifact()
        self.assertIn("spin four", artifact["spin_boundary"])
        self.assertIn("does not follow", artifact["decision"]["closed"])


if __name__ == "__main__":
    unittest.main()
