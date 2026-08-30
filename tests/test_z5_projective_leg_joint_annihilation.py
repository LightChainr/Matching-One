from __future__ import annotations

from pathlib import Path
import sys
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_projective_leg_joint_annihilation import (  # noqa: E402
    MONOMIALS_2,
    joint_stack,
    maximum_volume_pivot,
    projective_q,
    schur_complement,
    transformed_basis,
)


def matrix_with_kernel(rng: np.random.Generator, kernel: np.ndarray) -> np.ndarray:
    if abs(kernel[-1]) < 1e-8:
        raise ValueError("synthetic normalization coefficient vanished")
    matrix = rng.normal(size=(12, 6)) + 1j * rng.normal(size=(12, 6))
    matrix[:, -1] = -(matrix[:, :5] @ kernel[:5]) / kernel[-1]
    return matrix


class Z5ProjectiveLegJointAnnihilationTests(unittest.TestCase):
    def test_conjugating_bridge_is_one_complex_rank_five_null(self) -> None:
        rng = np.random.default_rng(250505)
        plus_kernel = rng.normal(size=6) + 1j * rng.normal(size=6)
        plus = matrix_with_kernel(rng, plus_kernel)
        minus = matrix_with_kernel(rng, plus_kernel.conjugate())
        joint = joint_stack(plus, minus, coefficient_conjugation=True)
        pivot = maximum_volume_pivot(joint, 5)
        self.assertLess(np.max(np.abs(schur_complement(joint, pivot))), 1e-10)
        fitted = projective_q(joint, pivot)
        self.assertLess(np.max(np.abs(joint @ fitted)), 1e-10)

    def test_direction_mismatch_is_visible_after_individual_rank_five(self) -> None:
        rng = np.random.default_rng(250506)
        plus_kernel = rng.normal(size=6) + 1j * rng.normal(size=6)
        minus_kernel = rng.normal(size=6) + 1j * rng.normal(size=6)
        plus = matrix_with_kernel(rng, plus_kernel)
        minus = matrix_with_kernel(rng, minus_kernel)
        self.assertEqual(np.linalg.matrix_rank(plus), 5)
        self.assertEqual(np.linalg.matrix_rank(minus), 5)
        joint = joint_stack(plus, minus, coefficient_conjugation=False)
        self.assertEqual(np.linalg.matrix_rank(joint), 6)
        pivot = maximum_volume_pivot(joint, 5)
        self.assertGreater(np.max(np.abs(schur_complement(joint, pivot))), 1e-8)

    def test_identity_linear_stack_is_literal_vertical_stack(self) -> None:
        rng = np.random.default_rng(250507)
        plus = rng.normal(size=(12, 6)) + 1j * rng.normal(size=(12, 6))
        minus = rng.normal(size=(12, 6)) + 1j * rng.normal(size=(12, 6))
        observed = joint_stack(plus, minus, coefficient_conjugation=False)
        self.assertTrue(np.array_equal(observed, np.vstack((plus, minus))))

    def test_alexander_basis_transforms_rows_and_columns_in_one_order(self) -> None:
        expected = ((0, 0), (-1, 0), (0, 1), (-2, 0), (-1, 1), (0, 2))
        observed = transformed_basis(alexander_reflection=True, rotation_power=2)
        self.assertEqual(observed, expected)
        self.assertEqual(len(observed), len(MONOMIALS_2))


if __name__ == "__main__":
    unittest.main()
