
from __future__ import annotations
import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from control_variate_estimator import (  # noqa: E402
    DuplicateChannelError,
    FrozenControlVariate,
    FrozenEstimator,
)
from euler_motif_controls import (  # noqa: E402
    ALL_MOTIFS,
    WRAPPING_CHANNELS,
    configuration_identity,
    exhaustive_conditional_means,
    exhaustive_identity,
    microcanonical_motif_mean,
    named_tiny_geometries,
    run_identity_suite,
)
from integer_period_torus import axis_integer_torus, gaussian_integer_torus  # noqa: E402


class EulerIdentityTests(unittest.TestCase):
    def test_prescribed_tiny_identities(self) -> None:
        payload = run_identity_suite()
        self.assertTrue(payload["passed"])
        names = {row["name"] for row in payload["exhaustive_identities"]}
        self.assertIn("axis", names)
        self.assertTrue(any(row["name"].startswith("gaussian-2-1") for row in payload["exhaustive_identities"]))
        for row in payload["exhaustive_identities"]:
            self.assertEqual(row["identity_failures"], 0)
            self.assertEqual(row["wrapping_not_identical"], 0)
            self.assertEqual(row["q_values"], [-1, 0, 1])

    def test_empty_and_full_axis_l2_boundary(self) -> None:
        geometry = axis_integer_torus(2)
        empty = configuration_identity(geometry, [False] * 4, 0)
        self.assertEqual(empty.q, -1)
        self.assertEqual(empty.cluster_difference, -1)
        self.assertEqual(empty.motifs["V"], 0)
        self.assertEqual(empty.residual, 0)
        full = configuration_identity(geometry, [True] * 4, 15)
        self.assertEqual(full.q, 1)
        self.assertEqual(full.cluster_difference, 1)
        self.assertEqual(full.motifs["V"], 4)
        self.assertEqual(full.motifs["E"], 8)
        self.assertEqual(full.motifs["F0"], 4)
        self.assertEqual(full.residual, 0)

    def test_conditional_means_are_exactly_zero_after_centering(self) -> None:
        for geometry in named_tiny_geometries():
            result = exhaustive_conditional_means(geometry)
            self.assertTrue(result["passed"], result["failures"])
            for row in result["by_K"]:
                for name, payload in row["motifs"].items():
                    self.assertTrue(payload["zero_conditional_mean_of_centered"])
                    self.assertAlmostEqual(
                        payload["expected_mean"],
                        microcanonical_motif_mean(name, geometry.n, row["K"]),
                    )


class DuplicateWrappingRejectionTests(unittest.TestCase):
    def test_exhaustive_tiny_wrapping_channels_are_identical(self) -> None:
        geometry = gaussian_integer_torus(2, 1)
        for mask in range(1 << geometry.n):
            active = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
            wrapping = configuration_identity(geometry, active, mask).wrapping
            self.assertEqual(len(set(wrapping.values())), 1)
            self.assertEqual(set(wrapping), set(WRAPPING_CHANNELS))

    def test_gls_rejects_duplicate_wrapping_channels(self) -> None:
        geometry = axis_integer_torus(3)
        rows = []
        for mask in range(1 << geometry.n):
            active = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
            wrapping = configuration_identity(geometry, active, mask).wrapping
            rows.append([float(wrapping[name]) for name in WRAPPING_CHANNELS])
        with self.assertRaisesRegex(DuplicateChannelError, "duplicate channels"):
            FrozenEstimator.fit(WRAPPING_CHANNELS, rows)

    def test_control_variate_fit_rejects_cloned_wrapping_controls(self) -> None:
        target = [0.0, 1.0, -1.0, 0.0]
        controls = [[1.0, 1.0], [0.0, 0.0], [-1.0, -1.0], [0.0, 0.0]]
        with self.assertRaises(DuplicateChannelError):
            FrozenControlVariate.fit(
                "q", ("q_cross", "q_either"), target, controls, (0.0, 0.0)
            )

    def test_ols_recovers_known_control_weight(self) -> None:
        # Noise [1,1,-1,-1] is orthogonal to z=[1,-1,1,-1] on this four-point design.
        target = [1.0, 1.0, -1.0, -1.0]
        controls = [[1.0], [-1.0], [1.0], [-1.0]]
        estimator = FrozenControlVariate.fit("q", ("Z",), target, controls, (0.0,))
        self.assertAlmostEqual(estimator.weights[0], 0.0)
        result = estimator.evaluate(target, controls, ("Z",))
        self.assertAlmostEqual(result["sample_variance"], 4.0 / 3.0)

        target = [3.0, -1.0, 1.0, -3.0]  # y = 2 z + previous orthogonal noise
        estimator = FrozenControlVariate.fit("q", ("Z",), target, controls, (0.0,))
        self.assertAlmostEqual(estimator.weights[0], 2.0)


if __name__ == "__main__":
    unittest.main()
