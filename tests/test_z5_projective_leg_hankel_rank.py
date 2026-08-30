from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_projective_leg_hankel_rank import (  # noqa: E402
    MONOMIALS_2,
    maximum_volume_pivot,
    schur_complement,
)


def moment_matrix(tx: np.ndarray, ty: np.ndarray, left: np.ndarray, right: np.ndarray) -> np.ndarray:
    def moment(point):
        return left @ np.linalg.matrix_power(tx, point[0]) @ np.linalg.matrix_power(ty, point[1]) @ right

    return np.asarray([[moment((u[0] + v[0], u[1] + v[1])) for v in MONOMIALS_2] for u in MONOMIALS_2])


class Z5ProjectiveLegHankelRankTests(unittest.TestCase):
    def test_rank_chart_includes_commuting_jordan_state(self) -> None:
        nilpotent = np.zeros((3, 3), dtype=complex)
        nilpotent[0, 1] = 1.0
        nilpotent[1, 2] = 1.0
        tx = 0.61 * np.eye(3) + nilpotent
        ty = -0.22 * np.eye(3) + 0.37 * nilpotent + 0.11 * nilpotent @ nilpotent
        self.assertLess(np.max(np.abs(tx @ ty - ty @ tx)), 1e-14)
        left = np.asarray([1.0 + 0.2j, -0.4j, 0.7])
        right = np.asarray([0.3, 1.1 - 0.1j, -0.8j])
        matrix = moment_matrix(tx, ty, left, right)
        pivot = maximum_volume_pivot(matrix, 3)
        self.assertLess(np.max(np.abs(schur_complement(matrix, pivot))), 1e-9)

    def test_generic_six_state_rejects_rank_five_algebraically(self) -> None:
        rng = np.random.default_rng(0)
        xroots = np.exp(2j * np.pi * rng.random(6))
        yroots = np.exp(2j * np.pi * rng.random(6))
        amplitudes = np.ones(6, dtype=complex)
        matrix = np.asarray([
            [sum(amplitudes * xroots ** (u[0] + v[0]) * yroots ** (u[1] + v[1])) for v in MONOMIALS_2]
            for u in MONOMIALS_2
        ])
        self.assertGreater(abs(np.linalg.det(matrix)), 1e-12)
        pivot = maximum_volume_pivot(matrix, 5)
        self.assertGreater(np.max(np.abs(schur_complement(matrix, pivot))), 1e-10)


if __name__ == "__main__":
    unittest.main()
