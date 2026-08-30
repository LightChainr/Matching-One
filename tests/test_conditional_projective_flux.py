from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_conditional_projective_flux import (  # noqa: E402
    align_spin4,
    conditional_values,
    score,
)


class ConditionalProjectiveFluxTests(unittest.TestCase):
    def test_factorized_line_timing_has_zero_discriminator(self) -> None:
        values = conditional_values({
            "birth0": 4.0,
            "exit0": 2.0,
            "direct": 0.5,
            "birth4_re": 1.0,
            "birth4_im": -2.0,
            "exit4_re": 0.5,
            "exit4_im": -1.0,
        })
        self.assertEqual(values["delta_mu_re"], 0.0)
        self.assertEqual(values["delta_mu_im"], 0.0)
        self.assertEqual(values["direct_flux"], 0.5)

    def test_gaussian_phase_alignment_is_exact(self) -> None:
        phase = complex(8, 1) ** 4 / 65**2
        aligned = align_spin4(
            {"delta_mu_re": phase.real * 0.125,
             "delta_mu_im": phase.imag * 0.125},
            8, 1,
        )
        self.assertAlmostEqual(aligned["delta_parallel"], 0.125, places=14)
        self.assertAlmostEqual(aligned["delta_perpendicular"], 0.0, places=14)

    def test_n65_archive_regression_and_full_covariance(self) -> None:
        result_root = ROOT / "results/local-20260830/P334-projective-birth-N65-smoke"
        result = score(
            result_root / "n65_20k.births.csv",
            result_root / "n65_20k.metadata.json",
            0.592746050790,
        )
        contrast = result["mechanism_gate"]["same_modulus_orientation_contrast"]
        self.assertAlmostEqual(contrast["value"], 0.011715329989237022)
        self.assertGreater(contrast["quadratic_against_equal_aligned_sorting"], 7.0)
        order = result["joint_estimate"]["order"]
        covariance = result["joint_estimate"]["covariance"]
        self.assertEqual(len(covariance), len(order))
        self.assertTrue(all(len(row) == len(order) for row in covariance))
        self.assertTrue(all(math.isfinite(value) for row in covariance for value in row))
        self.assertEqual(
            result["q_lift_boundary"]["classification"],
            "intrinsic horizontal Q=1 p-response",
        )


if __name__ == "__main__":
    unittest.main()

