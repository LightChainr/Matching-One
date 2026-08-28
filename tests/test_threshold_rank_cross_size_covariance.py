from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_threshold_rank_covariance import (  # noqa: E402
    _orientation_batches,
    constant_heldout_audit,
    covariance_of_mean,
    jackknife_pseudovalues,
)
from validate_threshold_rank_covariance_archive import validate_archive  # noqa: E402


ARCHIVE = ROOT / "results" / "server-20260828" / "P33-cross-size-covariance"


class ThresholdRankCrossSizeCovarianceTests(unittest.TestCase):
    def test_covariance_is_for_the_mean_not_raw_batches(self) -> None:
        means, covariance = covariance_of_mean(
            [
                [1.0, 2.0],
                [2.0, 4.0],
                [3.0, 6.0],
                [4.0, 8.0],
            ]
        )
        self.assertEqual(means, [2.5, 5.0])
        self.assertAlmostEqual(covariance[0][0], 5.0 / 12.0)
        self.assertAlmostEqual(covariance[0][1], 5.0 / 6.0)
        self.assertAlmostEqual(covariance[1][1], 5.0 / 3.0)

    def test_jackknife_pseudovalues_preserve_full_center(self) -> None:
        pseudo = jackknife_pseudovalues(10.0, [9.0, 10.0, 11.0])
        self.assertEqual(pseudo, [12.0, 10.0, 8.0])
        self.assertAlmostEqual(sum(pseudo) / len(pseudo), 10.0)

    def test_constant_heldout_score_propagates_training_covariance(self) -> None:
        sizes = [65, 85, 130, 145, 170]
        values = [0.80, 0.82, 0.79, 0.81, 0.77]
        variances = [0.01**2, 0.012**2, 0.015**2, 0.02**2, 0.025**2]
        covariance = [
            [
                variances[i]
                if i == j
                else 0.1 * math.sqrt(variances[i] * variances[j])
                for j in range(len(sizes))
            ]
            for i in range(len(sizes))
        ]
        result = constant_heldout_audit(
            values,
            covariance,
            sizes,
            [65, 85, 130],
            [145, 170],
        )
        self.assertEqual(result["heldout_dof"], 2)
        self.assertAlmostEqual(sum(result["training_weights"]), 1.0)
        self.assertGreater(result["amplitude_se"], 0.0)
        self.assertTrue(math.isfinite(result["heldout_chi_square"]))
        residual_covariance = result["heldout_residual_covariance"]
        self.assertNotEqual(residual_covariance[0][1], 0.0)

    def test_alignment_requires_same_batches_for_every_size(self) -> None:
        def record(n: int, orientation: str, batch: int) -> dict:
            return {
                "n": n,
                "orientation": orientation,
                "batch": batch,
                "samples": 10,
            }

        records = {}
        for n in (65, 85):
            for orientation in ("first", "second"):
                for batch in (0, 1):
                    records[(n, orientation, batch)] = record(
                        n, orientation, batch
                    )
        sizes, batches, _grouped = _orientation_batches(records)
        self.assertEqual(sizes, [65, 85])
        self.assertEqual(batches, [0, 1])

        del records[(85, "second", 1)]
        with self.assertRaisesRegex(ValueError, "aligned batch ids"):
            _orientation_batches(records)

    def test_committed_archive_has_stable_design_and_positive_covariance(self) -> None:
        result = validate_archive(
            ARCHIVE / "batch_metrics.csv",
            ARCHIVE / "summary.json",
            1e12,
        )
        batch = result["batch_contract"]
        self.assertEqual(batch["sizes"], [65, 85, 130, 145, 170])
        self.assertEqual(batch["batch_count"], 100)
        self.assertEqual(batch["samples_per_size_batch"], 100000)
        root = result["summary_contract"]["metric_covariance_diagnostics"][
            "root_gap"
        ]
        self.assertLess(root["infinity_norm_condition"], 10.0)


if __name__ == "__main__":
    unittest.main()
