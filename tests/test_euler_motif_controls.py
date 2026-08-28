from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from euler_motif_controls import (  # noqa: E402
    canonical_channels,
    centered_cluster_estimator,
    configuration_observables,
    exact_validation_suite,
    microcanonical_control_means,
    pilot_frozen_monte_carlo,
    prepare_convention,
)
from integer_period_torus import axis_integer_torus  # noqa: E402


class EulerExactIdentityTests(unittest.TestCase):
    def test_axis_diamond_and_gaussian_exhaustive_suite(self) -> None:
        result = exact_validation_suite()
        summaries = {
            row["geometry"]: (
                row["N"],
                row["euler_identity_max_abs_residual"],
                row["all_fixed_K_centered_sums_zero"],
            )
            for row in result["geometries"]
        }
        self.assertEqual(
            summaries,
            {
                "axis": (9, 0, True),
                "diamond": (8, 0, True),
                "gaussian-2-1": (5, 0, True),
            },
        )

    def test_fixed_k_hypergeometric_means(self) -> None:
        means = microcanonical_control_means(9, 4)
        self.assertEqual(means["V"], Fraction(4, 1))
        self.assertEqual(means["E"], Fraction(3, 1))
        self.assertEqual(means["F0"], Fraction(1, 14))
        self.assertEqual(means["diag_11"], Fraction(3, 2))
        self.assertEqual(means["corner_3"], Fraction(3, 7))

    def test_centered_canonical_identity(self) -> None:
        geometry = axis_integer_torus(3)
        convention = prepare_convention(geometry)
        active = [index in (0, 1, 3, 4) for index in range(geometry.n)]
        observed = configuration_observables(convention, active)
        p = 0.59274605079210
        q, controls = canonical_channels(observed, geometry.n, p)
        cluster = centered_cluster_estimator(observed, geometry.n, p)
        self.assertAlmostEqual(cluster, q + controls[0] - controls[1] + controls[2])


class EulerPilotFrozenTests(unittest.TestCase):
    def test_reproducible_fresh_evaluation_variance_reduction(self) -> None:
        geometry = axis_integer_torus(4)
        first = pilot_frozen_monte_carlo(
            geometry, 0.59274605079210, 1000, 4000, 17
        )
        second = pilot_frozen_monte_carlo(
            geometry, 0.59274605079210, 1000, 4000, 17
        )
        self.assertEqual(first, second)
        self.assertEqual(first["best_single_estimator"], "D_cluster")
        self.assertGreater(first["variance_reduction_vs_best_single"], 2.0)
        self.assertEqual(
            list(first["predeclared_control_hierarchy"]),
            ["euler", "euler_plus_diagonal", "euler_plus_local_motifs"],
        )
        self.assertEqual(
            first["estimators"]["pilot_frozen_control_variate"]["samples"],
            4000,
        )


if __name__ == "__main__":
    unittest.main()
