from __future__ import annotations

from pathlib import Path
import sys
import unittest

try:
    import numpy as np
    import scipy  # noqa: F401
except ImportError as exc:  # pragma: no cover - lightweight CI contract
    raise unittest.SkipTest(f"optional scientific stack unavailable: {exc}")


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_projective_leg_augmented_joint import (  # noqa: E402
    ROWS_3,
    augmented_residual_from_blocks,
    covariance_score_from_components,
    transformed_rows,
)


class Z5ProjectiveLegAugmentedJointTests(unittest.TestCase):
    def test_all_d4_transforms_remain_on_degree_three(self) -> None:
        for alexander in (False, True):
            for power in range(4):
                candidate = {
                    "alexander_reflection": alexander,
                    "rotation_power": power,
                }
                rows = transformed_rows(candidate)
                self.assertEqual(len(set(rows)), 4)
                self.assertTrue(all(abs(a) + abs(b) == 3 for a, b in rows))

    def test_fixed_old_pivot_appends_every_extension_equation(self) -> None:
        rng = np.random.default_rng(250505)
        q = np.asarray([1 + 0j, 2 - 1j, -1 + 2j, 3 + 0.5j, -2j, 1 - 3j])
        old = rng.normal(size=(24, 6)) + 1j * rng.normal(size=(24, 6))
        extension = rng.normal(size=(16, 6)) + 1j * rng.normal(size=(16, 6))
        # Make every row obey the same q while keeping the first five columns invertible.
        old[:, 5] = -(old[:, :5] @ q[:5]) / q[5]
        extension[:, 5] = -(extension[:, :5] @ q[:5]) / q[5]
        pivot = {"rows": [0, 1, 2, 3, 4], "columns": [0, 1, 2, 3, 4]}
        residual, observed_q = augmented_residual_from_blocks(old, extension, pivot)
        self.assertEqual(residual.shape, (35,))
        self.assertLess(float(np.max(np.abs(residual))), 1e-10)
        self.assertLess(float(np.max(np.abs(old @ observed_q))), 1e-10)
        self.assertLess(float(np.max(np.abs(extension @ observed_q))), 1e-10)

    def test_independent_covariance_components_are_saved_and_added(self) -> None:
        point = [0.2, -0.1]
        old_deleted = [[0.19, -0.11], [0.21, -0.09], [0.20, -0.12], [0.18, -0.10]]
        fresh_deleted = [[0.20, -0.08], [0.22, -0.11], [0.19, -0.09], [0.21, -0.12]]
        result = covariance_score_from_components(
            point, old_deleted, fresh_deleted, eigen_relative_cutoff=1e-10,
        )
        old = np.asarray(result["old_influence_covariance"])
        fresh = np.asarray(result["fresh_influence_covariance"])
        total = np.asarray(result["total_covariance"])
        np.testing.assert_allclose(total, old + fresh, atol=0.0, rtol=0.0)
        self.assertLessEqual(result["covariance_identity_max_abs_error"], 1e-18)


if __name__ == "__main__":
    unittest.main()
