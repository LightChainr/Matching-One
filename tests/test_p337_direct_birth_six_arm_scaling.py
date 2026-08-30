from __future__ import annotations

import csv
import gzip
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p337_direct_birth_six_arm_scaling import (  # noqa: E402
    BETA_FIXED,
    extract_batches,
    model_scores,
    summarize_size,
)


class P337DirectBirthScalingTests(unittest.TestCase):
    def write_fixture(self, path: Path, *, compressed: bool = False) -> None:
        rows = []
        for orientation in ("first", "second"):
            for batch, direct in enumerate((10, 20)):
                rows += [
                    [10, 1, 1, orientation, batch, 100, 3, 3, "DIRECT_RANK2", 0, 0, direct],
                    [10, 1, 1, orientation, batch, 100, 2, 7, "LINE", 1, 0, 100 - direct],
                ]
        opener = gzip.open if compressed else open
        with opener(path, "wt", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["n", "a", "b", "orientation", "batch", "samples", "tau1", "tau2", "kind", "ell_x", "ell_y", "count"])
            writer.writerows(rows)

    def test_sparse_partition_and_paired_summary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "births.csv"
            self.write_fixture(path)
            spec = {"N": 10, "batches": 2, "samples_per_orientation": 200}
            reduced = extract_batches(path, spec)
            summary = summarize_size(reduced, spec)
            self.assertAlmostEqual(summary["Dbar"], 0.15)
            self.assertEqual(summary["direct_counts"], {"first": 30, "second": 30})
            self.assertEqual(summary["input_audit"]["partition_gate"], True)

    def test_gzip_input_has_identical_semantics(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "births.csv.gz"
            self.write_fixture(path, compressed=True)
            spec = {"N": 10, "batches": 2, "samples_per_orientation": 200, "compression": "gzip"}
            reduced = extract_batches(path, spec)
            self.assertEqual(int(reduced["direct"].sum()), 60)

    def test_exact_fixed_power_is_recovered(self) -> None:
        sizes = []
        for n in (85, 170, 340, 680):
            value = 0.6 * n ** (-BETA_FIXED)
            sizes.append({"N": n, "Dbar": value, "log_Dbar_variance_delta": 1e-6})
        models = model_scores(sizes, 0.01)
        self.assertAlmostEqual(models["fixed_5_6"]["amplitude"], 0.6, places=10)
        self.assertAlmostEqual(models["free_power"]["beta"], BETA_FIXED, places=10)
        self.assertAlmostEqual(models["minimal_log_curvature"]["kappa"], 0.0, places=10)


if __name__ == "__main__":
    unittest.main()
