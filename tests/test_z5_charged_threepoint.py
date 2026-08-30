from __future__ import annotations

import sys
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from norm5_chiral_hecke_phase import gaussian_ratio_power  # noqa: E402
from score_z5_charged_threepoint import gls_model  # noqa: E402
from z5_charged_threepoint_mc import (  # noqa: E402
    JOINT_REAL_ORDER,
    PRIMARY_REAL_ORDER,
    dft_charges,
    exact_mapping_gate,
    parent_anchor_indices,
    run,
)


class Z5ChargedThreePointTests(unittest.TestCase):
    def test_exact_mapping_and_anchor_gate(self):
        gate = exact_mapping_gate()
        self.assertTrue(gate["passed"])
        self.assertEqual(gate["parent_translations_checked"], 65)
        self.assertEqual(gate["root_labels_checked"], 1950)
        self.assertEqual(gate["projection_failures"], 0)
        self.assertEqual(gate["anchor_failures"], 0)

    def test_fiber_dft_conjugacy_and_charge_labels(self):
        values = (0, 1, -2, 3, -4)
        rows = dft_charges(values)
        self.assertAlmostEqual(rows[4], rows[1].conjugate())
        self.assertAlmostEqual(rows[3], rows[2].conjugate())

    def test_counter_anchor_schedule_is_deterministic_and_nondegenerate(self):
        first = parent_anchor_indices(25011312220260830, 17)
        self.assertEqual(first, parent_anchor_indices(25011312220260830, 17))
        self.assertEqual(len(set(first[1])), 3)
        self.assertNotEqual(first, parent_anchor_indices(25011312220260830, 18))

    def test_tiny_stream_retains_full_covariance_and_controls(self):
        batches, analysis = run(4, 2, 1, 0.592746050790, 25011312220260830, 1, 0)
        self.assertEqual(len(batches), 2)
        self.assertEqual(analysis["primary_order"], list(PRIMARY_REAL_ORDER))
        self.assertEqual(len(analysis["primary_covariance_of_mean"]), 8)
        self.assertEqual(len(analysis["joint_order"]), len(JOINT_REAL_ORDER))
        self.assertEqual(len(analysis["joint_covariance_of_mean"]), 24)
        self.assertLess(analysis["conjugacy_max_abs"], 1e-12)
        self.assertEqual(analysis["closure"]["delete_one_replicates"], 2)

    def test_synthetic_joint_gls_recovers_each_frozen_phase(self):
        covariance = [[1e-4 if i == j else 0.0 for j in range(8)] for i in range(8)]
        amplitudes = (complex(0.3, -0.2), complex(-0.1, 0.25))
        for spin in (4, 8, 12):
            real, imag, denominator = gaussian_ratio_power(3 * spin)
            q = complex(real / denominator, imag / denominator)
            mean = []
            for amplitude in amplitudes:
                mean.extend(((q * amplitude).real, (q * amplitude).imag, amplitude.real, amplitude.imag))
            own = gls_model(mean, covariance, q)
            self.assertLess(own["chi_square"], 1e-20)


if __name__ == "__main__":
    unittest.main()
