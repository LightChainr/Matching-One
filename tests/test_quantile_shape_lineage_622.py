#!/usr/bin/env python3
"""Lock the #622 reflection-residual pipeline: anchors, redundancies, numerics.

These are tiny exact-algebra checks, not the heavy lineage run.  They stop the
classes of mistake that would make ``A_N`` look like a shape signal when it is
an anchor artefact, a wrong coordinate, or an inverted near-singular covariance.
"""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import threshold_quantile_lineage as lineage  # noqa: E402

_SPEC = importlib.util.spec_from_file_location(
    "quantile_shape_lineage_622",
    ROOT / "scripts" / "probe_invariant_shape" / "quantile_shape_lineage_622.py")
assert _SPEC is not None and _SPEC.loader is not None
qsl = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(qsl)


class AnchorsAndRedundancy(unittest.TestCase):
    def test_the_anchor_values_are_exact(self) -> None:
        """Stops A being read at a coordinate that is not a hard zero."""
        levels = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
        Q = [0.5, 0.55, 0.56, 0.58, 0.593, 0.605, 0.618, 0.625, 0.65]
        Z, width = qsl.z_chart(Q, levels, 0.2, 0.8)
        self.assertAlmostEqual(Z[1], 0.0, places=15)
        self.assertAlmostEqual(Z[7], 1.0, places=15)
        self.assertGreater(width, 0.0)

    def test_reflected_pairs_repeat_and_anchor_is_zero(self) -> None:
        """Stops a redundant A coordinate being counted as a new measurement."""
        levels = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
        Z = [-0.27, 0.0, 0.19, 0.35, 0.50, 0.65, 0.81, 1.0, 1.25]
        A = qsl.reflection_residual(Z, levels)
        self.assertAlmostEqual(A[1], 0.0, places=15)
        self.assertAlmostEqual(A[7], 0.0, places=15)
        for i in range(9):
            self.assertAlmostEqual(A[i], A[8 - i], places=15)

    def test_the_pairing_is_the_declared_one(self) -> None:
        """Stops the raw-p reflection M(p)+M(1-p) being used instead of A."""
        levels = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
        Z = [-0.27, 0.0, 0.19, 0.35, 0.50, 0.65, 0.81, 1.0, 1.25]
        A = qsl.reflection_residual(Z, levels)
        # A(u) = Z(u) + Z(1-u) - 1, on the normalized chart.
        for i, u in enumerate(levels):
            self.assertAlmostEqual(A[i], Z[i] + Z[8 - i] - 1.0, places=15)


class PureNumerics(unittest.TestCase):
    def test_jacobi_matches_a_known_spectrum(self) -> None:
        eigenvalues = qsl.jacobi_eigenvalues([[2.0, 1.0], [1.0, 2.0]])
        self.assertAlmostEqual(eigenvalues[0], 1.0, places=12)
        self.assertAlmostEqual(eigenvalues[1], 3.0, places=12)

    def test_cholesky_solves_the_identity_case(self) -> None:
        solution = qsl.cholesky_solve([[2.0, 1.0], [1.0, 2.0]], [3.0, 3.0])
        self.assertAlmostEqual(solution[0], 1.0, places=12)
        self.assertAlmostEqual(solution[1], 1.0, places=12)

    def test_chi_square_survival_matches_published_values(self) -> None:
        self.assertAlmostEqual(qsl.chi2_survival(1.0, 4), 0.9097959895689501,
                               places=12)
        self.assertAlmostEqual(qsl.chi2_survival(20.0, 4), 0.0004993992273873336,
                               places=12)

    def test_a_singular_covariance_is_not_pseudoinverted(self) -> None:
        """Stops an exact redundancy being discarded by a pseudoinverse.

        The wrong number is a finite chi-square formed from a rank-deficient
        covariance, which hides the residual that contradicts the model.
        """
        levels = (0.1, 0.2, 0.3, 0.4, 0.5)
        covariance = [[1.0, 0.0, 1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0, 0.0, 0.0],
                      [1.0, 0.0, 1.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0, 1.0, 0.0],
                      [0.0, 0.0, 0.0, 0.0, 1.0]]
        vector = [0.3, 0.0, 0.3, 0.0, 0.1]
        diagnostic = qsl._quadratic_diagnostic(
            vector, covariance, (0.1, 0.3, 0.4, 0.5), levels)
        self.assertIsNone(diagnostic["chi2"])
        self.assertIn("no pseudoinverse", diagnostic["reason"])


class FiniteMovement(unittest.TestCase):
    def test_a_resolved_monotone_decrease_is_named(self) -> None:
        sizes = {
            "145": {"analysis": {"spin0": {"A_norm_on_independent_coordinates": {
                "value": 0.030, "se": 1e-5}}}},
            "290": {"analysis": {"spin0": {"A_norm_on_independent_coordinates": {
                "value": 0.020, "se": 1e-5}}}},
            "725": {"analysis": {"spin0": {"A_norm_on_independent_coordinates": {
                "value": 0.010, "se": 1e-5}}}},
        }
        adjacent = {
            "145->290": {"spin0": {"Z_norm": {"value": 0.005, "se": 1e-5},
                                   "A_norm": {"value": 0.007, "se": 1e-5}}},
            "290->725": {"spin0": {"Z_norm": {"value": 0.004, "se": 1e-5},
                                   "A_norm": {"value": 0.006, "se": 1e-5}}},
        }
        move = qsl.finite_movement(sizes, adjacent, "spin0")
        self.assertEqual(move["A_norm_trend"], "decreasing")
        self.assertTrue(move["adjacent_displacement_interval_change"]
                        ["delta_Z_norm"]["resolved_change"])

    def test_an_unresolved_step_is_not_called_a_trend(self) -> None:
        sizes = {
            "145": {"analysis": {"spin0": {"A_norm_on_independent_coordinates": {
                "value": 0.0300, "se": 1e-5}}}},
            "290": {"analysis": {"spin0": {"A_norm_on_independent_coordinates": {
                "value": 0.02998, "se": 1e-5}}}},
            "725": {"analysis": {"spin0": {"A_norm_on_independent_coordinates": {
                "value": 0.0100, "se": 1e-5}}}},
        }
        adjacent = {
            "145->290": {"spin0": {"Z_norm": {"value": 0.005, "se": 1e-5},
                                   "A_norm": {"value": 0.007, "se": 1e-5}}},
            "290->725": {"spin0": {"Z_norm": {"value": 0.005, "se": 1e-5},
                                   "A_norm": {"value": 0.007, "se": 1e-5}}},
        }
        move = qsl.finite_movement(sizes, adjacent, "spin0")
        self.assertEqual(move["A_norm_trend"], "unresolved_on_this_lineage")


class LineageWeights(unittest.TestCase):
    def test_the_three_committed_pairs_are_interpolating_spin_zero(self) -> None:
        """Stops the corrected cos 4theta from being replaced by cos 2theta."""
        pairs = {145: ((12, 1), (9, 8)),
                 290: ((13, 11), (17, 1)),
                 725: ((26, 7), (23, 14))}
        for size, (first, second) in pairs.items():
            cos4 = {"first": lineage.cos_four_theta(*first),
                    "second": lineage.cos_four_theta(*second)}
            weights = lineage.spin_zero_weights(cos4)
            self.assertAlmostEqual(weights["first"] + weights["second"], 1.0,
                                   places=12)
            self.assertAlmostEqual(weights["first"] * cos4["first"]
                                   + weights["second"] * cos4["second"], 0.0,
                                   places=12)
            self.assertTrue(lineage.is_interpolation(weights), size)


if __name__ == "__main__":
    unittest.main()
