from __future__ import annotations

from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from covariance_nullspace import (  # noqa: E402
    CovarianceNullspaceViolation,
    covariance_spectral_diagnostics,
)


class CovarianceNullspaceTests(unittest.TestCase):
    def test_exact_null_compatible_residual(self) -> None:
        score = covariance_spectral_diagnostics(
            [1, 1], [[1, 1], [1, 1]], nullspace_policy="structural"
        )
        self.assertEqual(score["nullspace_status"], "structural_null_compatible")
        self.assertTrue(score["nullspace_compatible"])
        self.assertAlmostEqual(float(score["chi_square"]), 1.0)

    def test_exact_null_violation_fails_closed(self) -> None:
        with self.assertRaisesRegex(
            CovarianceNullspaceViolation, "structural covariance null"
        ):
            covariance_spectral_diagnostics(
                [1, -1], [[1, 1], [1, 1]], nullspace_policy="structural"
            )

    def test_estimated_null_violation_is_explicit(self) -> None:
        score = covariance_spectral_diagnostics([1, -1], [[1, 1], [1, 1]])
        self.assertEqual(
            score["nullspace_status"], "estimated_near_null_incompatibility"
        )
        self.assertFalse(score["nullspace_compatible"])
        self.assertEqual(len(score["discarded_eigendirections"]), 1)
        self.assertAlmostEqual(
            float(score["max_abs_discarded_standardized_residual"]),
            2 ** 0.5,
        )

    def test_full_rank_statistic_is_unchanged(self) -> None:
        score = covariance_spectral_diagnostics([1, 2], [[1, 0], [0, 4]])
        self.assertEqual(score["nullspace_status"], "full_rank")
        self.assertAlmostEqual(float(score["chi_square"]), 2.0)
        self.assertEqual(score["degrees_of_freedom"], 2)

    def test_near_singular_cutoff_sensitivity_is_frozen(self) -> None:
        epsilon = mp.mpf("1e-12")
        score = covariance_spectral_diagnostics(
            [1, -1], [[1, 1 - epsilon], [1 - epsilon, 1]]
        )
        rows = {
            mp.mpf(row["relative_eigenvalue_cutoff"]): row
            for row in score["cutoff_sensitivity"]
        }
        self.assertEqual(rows[mp.mpf("1e-14")]["numerical_rank"], 2)
        self.assertEqual(rows[mp.mpf("1e-10")]["numerical_rank"], 1)
        self.assertEqual(rows[mp.mpf("1e-6")]["numerical_rank"], 1)
        self.assertGreater(
            rows[mp.mpf("1e-14")]["chi_square"], mp.mpf("1e10")
        )

    def test_materially_indefinite_covariance_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "materially indefinite"):
            covariance_spectral_diagnostics([1, 1], [[1, 2], [2, 1]])

    def test_raw_covariance_allows_a_structural_zero_variance_coordinate(self) -> None:
        score = covariance_spectral_diagnostics(
            [1, 0], [[1, 0], [0, 0]], standardize=False
        )
        self.assertEqual(score["spectral_basis"], "raw_covariance")
        self.assertEqual(score["numerical_rank"], 1)
        self.assertAlmostEqual(float(score["chi_square"]), 1.0)
        self.assertEqual(score["nullspace_status"], "estimated_near_null_compatible")

    def test_raw_nullspace_compatibility_is_scale_invariant(self) -> None:
        maxima = []
        for scale in (mp.mpf(1), mp.mpf("1e-13")):
            score = covariance_spectral_diagnostics(
                [scale, -scale],
                [
                    [scale ** 2, scale ** 2],
                    [scale ** 2, scale ** 2],
                ],
                standardize=False,
            )
            self.assertEqual(
                score["nullspace_status"], "estimated_near_null_incompatibility"
            )
            maxima.append(score["max_abs_discarded_nullspace_projection"])
        self.assertAlmostEqual(float(maxima[0]), float(maxima[1]), places=12)


if __name__ == "__main__":
    unittest.main()
