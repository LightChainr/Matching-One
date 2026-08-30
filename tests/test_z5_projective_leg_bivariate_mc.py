from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from z5_projective_leg_bivariate_mc import (  # noqa: E402
    FIELD_ORDER,
    FIRST_QUADRANT,
    GRID,
    exact_gate,
    label,
    rotate,
    rotation_gauges,
)


class Z5ProjectiveLegBivariateTests(unittest.TestCase):
    def test_grid_is_c4_closed_degree_four_diamond(self) -> None:
        self.assertEqual(len(GRID), 41)
        self.assertEqual(len(FIRST_QUADRANT), 15)
        self.assertEqual(set(map(rotate, GRID)), set(GRID))
        self.assertTrue({(1, 1), (2, 1), (1, 2)}.issubset(GRID))
        self.assertEqual(len(FIELD_ORDER), 41 * 2 * 2 * 2)

    def test_labels_are_sign_explicit(self) -> None:
        self.assertEqual(label(2, -1), "ap2_bm1")
        self.assertEqual(label(-2, 1), "am2_bp1")

    def test_exact_gate(self) -> None:
        gate = exact_gate()
        self.assertTrue(gate["passed"])
        self.assertEqual(gate["rotation_fiber_gate"]["fiber_multipliers_mod5"], {"plus": 3, "minus": 2})
        self.assertEqual(gate["rotation_fiber_gate"]["channel_map_fourth_power"], "identity")
        gauges = rotation_gauges()
        self.assertEqual(set(gauges), {"plus", "minus"})
        self.assertTrue(all(len(row) == 101 for row in gauges.values()))


if __name__ == "__main__":
    unittest.main()
