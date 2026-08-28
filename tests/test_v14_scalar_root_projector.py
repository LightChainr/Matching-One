from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_v14_scalar_root_projector import lineage_stat  # noqa: E402
from v14_scalar_post_l7 import (  # noqa: E402
    conditional_interchiral_parity,
    critical_potts_h,
    diagonal_x,
    root_bias_exponent_in_L,
)


def scalar_project(c1: float, c2: float, p1: float, p2: float) -> float:
    return (c1 * p2 - c2 * p1) / (c1 - c2)


class V14ScalarRootProjectorTests(unittest.TestCase):
    def test_potts_arithmetic_and_l7_exponent(self) -> None:
        self.assertEqual(critical_potts_h(1, 2), F(5, 8))
        self.assertEqual(diagonal_x(4), F(33, 4))
        self.assertEqual(conditional_interchiral_parity(4), -1)
        self.assertEqual(root_bias_exponent_in_L(F(33, 4)), F(7, 1))
        self.assertEqual(diagonal_x(4) - F(21, 4), F(3, 1))

    def test_two_angle_projector_exactly_cancels_h4_root_term(self) -> None:
        pc = 0.59274605079
        n = 65
        beta0 = 3.5
        a0 = 3.75
        a4 = -0.41
        c1, c2 = 0.87, -0.49
        common = pc + a0 * n ** (-beta0)
        p1 = common + a4 * c1 * n ** (-2.0)
        p2 = common + a4 * c2 * n ** (-2.0)
        projected = scalar_project(c1, c2, p1, p2)
        self.assertAlmostEqual(projected, common, places=15)
        self.assertAlmostEqual((p1 - p2) / (c1 - c2), a4 * n ** (-2.0), places=15)

    def test_norm2_lineage_recovers_pc_and_scalar_amplitude_without_external_pc(self) -> None:
        pc = 0.59274605079
        beta = 3.5
        amplitude = -2.75
        sample = {
            65: {"p_scalar": pc + amplitude * 65 ** (-beta)},
            130: {"p_scalar": pc + amplitude * 130 ** (-beta)},
            85: {"p_scalar": pc + amplitude * 85 ** (-beta)},
            170: {"p_scalar": pc + amplitude * 170 ** (-beta)},
        }
        first = lineage_stat(sample, 65, 130, beta)
        second = lineage_stat(sample, 85, 170, beta)
        self.assertAlmostEqual(first["pc_hat"], pc, places=15)
        self.assertAlmostEqual(second["pc_hat"], pc, places=15)
        self.assertAlmostEqual(first["amplitude"], amplitude, places=9)
        self.assertAlmostEqual(second["amplitude"], amplitude, places=9)

    def test_wrong_beta_breaks_cross_lineage_pc_consistency(self) -> None:
        pc = 0.59274605079
        true_beta = 3.5
        amplitude = 5.0
        sample = {
            65: {"p_scalar": pc + amplitude * 65 ** (-true_beta)},
            130: {"p_scalar": pc + amplitude * 130 ** (-true_beta)},
            85: {"p_scalar": pc + amplitude * 85 ** (-true_beta)},
            170: {"p_scalar": pc + amplitude * 170 ** (-true_beta)},
        }
        first = lineage_stat(sample, 65, 130, 3.0)
        second = lineage_stat(sample, 85, 170, 3.0)
        self.assertGreater(abs(first["pc_hat"] - second["pc_hat"]), 1e-10)


if __name__ == "__main__":
    unittest.main()
