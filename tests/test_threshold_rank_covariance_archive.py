from __future__ import annotations

import copy
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from validate_threshold_rank_covariance_archive import (  # noqa: E402
    covariance_diagnostics,
    read_batch_rows,
    validate_archive,
    validate_batch_rows,
)


ARCHIVE = ROOT / "results" / "server-20260828" / "P33-cross-size-covariance"


class ThresholdRankCovarianceArchiveTests(unittest.TestCase):
    def test_committed_archive_satisfies_contracts(self) -> None:
        result = validate_archive(
            ARCHIVE / "batch_metrics.csv",
            ARCHIVE / "summary.json",
            1e12,
        )
        batch = result["batch_contract"]
        self.assertEqual(batch["sizes"], [65, 85, 130, 145, 170])
        self.assertEqual(batch["batch_count"], 100)
        self.assertEqual(batch["samples_per_size_batch"], 100000)
        metrics = result["summary_contract"]["metric_covariance_diagnostics"]
        self.assertIn("root_gap", metrics)
        self.assertLess(metrics["root_gap"]["infinity_norm_condition"], 10.0)

    def test_batch_contract_rejects_sample_or_geometry_drift(self) -> None:
        rows = read_batch_rows(ARCHIVE / "batch_metrics.csv")[:10]

        changed_samples = copy.deepcopy(rows)
        changed_samples[1]["samples"] = "99999"
        with self.assertRaisesRegex(ValueError, "one sample count"):
            validate_batch_rows(changed_samples)

        changed_geometry = copy.deepcopy(rows)
        changed_geometry[5]["a1"] = "7"
        changed_geometry[5]["b1"] = "4"
        with self.assertRaisesRegex(ValueError, "Gaussian norm|changes geometry"):
            validate_batch_rows(changed_geometry)

    def test_covariance_contract_rejects_indefinite_or_ill_conditioned_input(self) -> None:
        with self.assertRaisesRegex(ValueError, "positive definite"):
            covariance_diagnostics(
                [[1.0, 2.0], [2.0, 1.0]],
                "indefinite",
                1e12,
            )
        with self.assertRaisesRegex(ValueError, "condition number"):
            covariance_diagnostics(
                [[1.0, 0.0], [0.0, 1e-14]],
                "ill-conditioned",
                1e12,
            )


if __name__ == "__main__":
    unittest.main()
