#!/usr/bin/env python3
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p250_adaptive_order_pilot as pilot  # noqa: E402


class AdaptiveOrderPilotTests(unittest.TestCase):
    def test_frozen_geometries(self):
        self.assertEqual(pilot.MediumOracle("N325").geometry.n, 325)
        self.assertEqual(pilot.MediumOracle("N425").geometry.n, 425)

    def test_dihedral_stencil_has_eight_covariant_images(self):
        images = {
            (
                pilot.transform_offset(orientation, (1, 1)),
                pilot.transform_offset(orientation, (0, 1)),
            )
            for orientation in range(8)
        }
        self.assertEqual(len(images), 8)

    def test_counter_field_and_marks_are_reproducible(self):
        oracle = pilot.MediumOracle("N325")
        first = oracle.marks(123, 456)
        second = oracle.marks(123, 456)
        self.assertEqual(first, second)
        self.assertEqual(
            pilot.counter_uniform(123, 456, 7, 325),
            pilot.counter_uniform(123, 456, 7, 325),
        )

    def test_tiny_medium_smoke_keeps_exact_controls(self):
        row = pilot.run_batch(("N325", 0, 10_250_000_000, 40, pilot.P_REF, 12345))
        self.assertEqual(row["typed_defined_mismatch"], 0)
        self.assertEqual(row["typed_support_mismatch"], 0)
        self.assertEqual(row["typed_Rminus_residual_max"], 0)
        self.assertEqual(row["typed_Rplus_sum_max"], 0)
        self.assertEqual(row["fixed_support_order_null_failures"], 0)


if __name__ == "__main__":
    unittest.main()
