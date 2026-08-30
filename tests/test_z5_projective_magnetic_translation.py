from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_projective_magnetic_translation import (  # noqa: E402
    canonical_gate,
    spatial_bundle_flatness_gate,
    weyl_displacement,
)


class Z5ProjectiveMagneticTranslationTests(unittest.TestCase):
    def test_canonical_weyl_and_weil_relations(self) -> None:
        gate = canonical_gate()
        self.assertTrue(gate["passed"])
        for row in gate["models"].values():
            self.assertLess(row["D_center_fifth_power_minus_I"], 1e-12)

    def test_spatial_cover_bundle_is_flat(self) -> None:
        gate = spatial_bundle_flatness_gate()
        self.assertTrue(gate["passed"])
        self.assertEqual(gate["counts"]["plus"]["0"], 101)
        self.assertEqual(gate["counts"]["minus"]["0"], 101)

    def test_symmetric_weyl_adjoint(self) -> None:
        for m in range(1, 5):
            for point in ((1, 1), (2, 1), (-1, 2)):
                opposite = (-point[0], -point[1])
                self.assertLess(np.max(np.abs(weyl_displacement(point, m).conjugate().T - weyl_displacement(opposite, m))), 1e-12)


if __name__ == "__main__":
    unittest.main()
