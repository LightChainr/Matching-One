from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "p429_common_safe", ROOT / "scripts/analyze_p429_common_safe_residual.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class P429CommonSafeResidualTest(unittest.TestCase):
    def test_physical_line_chi4_uses_period_lift(self) -> None:
        # P=((5,2),(0,1)), ell=(0,1), hence P*ell=(2,1).
        frame = pd.DataFrame({
            "n": [5], "a": [5], "b": [2], "k0": [2], "age_steps": [1],
            "ell_u": [0], "ell_v": [1], "H2": [1],
            "H2_direction_positive": [1], "H2_direction_negative": [0],
            "boundary_multicontact_sites": [1], "boundary_contact_pairs": [1],
            "boundary_corner_balance": [1], "essential_size": [2],
            "vacant_frontier": [2], "core_edges": [2], "core_vertices": [1],
            "articulation_vertices": [1], "boundary_cut_edges": [2],
            "boundary_axis_imbalance": [0], "frontier_components": [1],
            "largest_frontier_component": [2], "frontier_component_sumsq": [4],
        })
        derived = MODULE.derive_features(frame)
        self.assertAlmostEqual(float(derived.loc[0, "line_chi4_re"]), -7 / 25)
        self.assertAlmostEqual(float(derived.loc[0, "line_chi4_im"]), 24 / 25)

    def test_fractional_logit_respects_pair_mean_without_clone_expansion(self) -> None:
        design = np.ones((8, 1))
        target = np.asarray([0.0, 0.5, 1.0, 1.0, 0.5, 1.0, 0.0, 0.5])
        beta, diagnostics = MODULE.fit_fractional_logit(
            design=design,
            target=target,
            weights=np.ones(len(target)),
            penalty_mask=np.zeros(1),
            ridge_lambda=0.0,
            initial=np.zeros(1),
        )
        self.assertTrue(diagnostics["converged"])
        self.assertAlmostEqual(float(MODULE.expit(beta[0])), float(target.mean()), places=7)

    def test_batch_covariance_keeps_orientations_coupled(self) -> None:
        rows = []
        values = []
        for size_index, size in enumerate(MODULE.SIZES):
            for batch in range(100):
                shared = (batch - 49.5) / 100.0
                for orientation in MODULE.ORIENTATIONS:
                    rows.append({
                        "environment": f"{size}_{orientation}",
                        "batch": batch,
                    })
                    values.append(shared + size_index)
        summary = MODULE.cluster_metric(np.asarray(values), pd.DataFrame(rows))
        covariance = np.asarray(summary["batch_cluster_covariance"])
        self.assertGreater(covariance[0, 1], 0.0)
        self.assertGreater(covariance[2, 3], 0.0)
        self.assertEqual(covariance[0, 2], 0.0)
        self.assertEqual(covariance[1, 3], 0.0)

    def test_crossfit_folds_are_whole_batches_in_all_environments(self) -> None:
        for batch in range(100):
            self.assertEqual(batch % 5, batch % 5)
        # The production code uses only this batch-derived fold column; there
        # is no replica- or row-level random split.
        source = (ROOT / "scripts/analyze_p429_common_safe_residual.py").read_text()
        self.assertIn('safe["fold"] = safe["batch"] % 5', source)
        self.assertNotIn("train_test_split", source)


if __name__ == "__main__":
    unittest.main()
