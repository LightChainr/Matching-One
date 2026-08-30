from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_projective_leg_radius6_flat import (  # noqa: E402
    ideal_bridge_residual,
    kernel_projector,
    schur,
)


class Z5ProjectiveLegRadius6FlatScoreTests(unittest.TestCase):
    def test_fixed_chart_schur_vanishes_for_rank_five_matrix(self) -> None:
        rng = np.random.default_rng(250606)
        left = rng.normal(size=(20, 5)) + 1j * rng.normal(size=(20, 5))
        right = rng.normal(size=(5, 10)) + 1j * rng.normal(size=(5, 10))
        matrix = left @ right
        residual = schur(matrix, tuple(range(5)), tuple(range(5)))
        self.assertLess(float(np.max(np.abs(residual))), 1e-10)

    def test_conjugate_kernel_projector_bridge_is_basis_independent(self) -> None:
        rng = np.random.default_rng(250607)
        matrix = rng.normal(size=(20, 10)) + 1j * rng.normal(size=(20, 10))
        plus = matrix[:, :5] @ matrix[:5, :]
        minus = plus.conjugate()
        self.assertLess(max(map(abs, ideal_bridge_residual(plus, minus))), 1e-10)
        self.assertAlmostEqual(float(np.trace(kernel_projector(plus)).real), 5.0, places=10)


if __name__ == "__main__":
    unittest.main()
