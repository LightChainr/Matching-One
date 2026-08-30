from __future__ import annotations

from itertools import combinations
from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_projective_leg_augmented_map import (  # noqa: E402
    centered_influences,
    maximum_volume_pivot,
    score_from_influences,
)
from score_z5_projective_leg_radius5_morphism import transformed_rows  # noqa: E402
from score_z5_projective_leg_annihilator_bridge import transformed_basis  # noqa: E402


class P250AugmentedMapTests(unittest.TestCase):
    def test_vectorized_pivot_matches_brute_force(self) -> None:
        rng = np.random.default_rng(250)
        matrix = rng.normal(size=(8, 6)) + 1j * rng.normal(size=(8, 6))
        observed = maximum_volume_pivot(matrix)
        best = None
        for columns in combinations(range(6), 5):
            for rows in combinations(range(8), 5):
                volume = float(abs(np.linalg.det(matrix[np.ix_(rows, columns)])))
                if best is None or volume > best[0]:
                    best = (volume, rows, columns)
        self.assertAlmostEqual(observed["abs_determinant"], best[0], places=12)
        self.assertEqual(observed["rows"], best[1])
        self.assertEqual(observed["columns"], best[2])

    def test_same_d4_action_reaches_basis_and_left_shifts(self) -> None:
        for power in range(4):
            basis = transformed_basis(alexander_reflection=True, rotation_power=power)
            rows = transformed_rows(True, power)
            ximage, yimage = basis[1], basis[2]
            expected = tuple(
                (u[0] * ximage[0] + u[1] * yimage[0], u[0] * ximage[1] + u[1] * yimage[1])
                for u in ((3, 0), (2, 1), (1, 2), (0, 3))
            )
            self.assertEqual(rows, expected)

    def test_centered_source_influences_add_not_pair(self) -> None:
        old_deleted = [[1.0, 0.0], [3.0, 2.0], [2.0, 1.0]]
        fresh_deleted = [[5.0, 2.0], [4.0, 4.0], [6.0, 3.0]]
        old = centered_influences(old_deleted)
        fresh = centered_influences(fresh_deleted)
        result = score_from_influences([0.2, -0.1], old, fresh, effective_batches=20)
        expected = old.T @ old + fresh.T @ fresh
        np.testing.assert_allclose(result["covariance"], expected)
        self.assertEqual(result["covariance_addition_max_abs_error"], 0.0)


if __name__ == "__main__":
    unittest.main()
