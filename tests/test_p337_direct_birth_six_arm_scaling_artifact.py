from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/local-20260830/P337-direct-birth-six-arm-scaling"


def jackknife_covariance(values: np.ndarray) -> np.ndarray:
    centered = values - values.mean(axis=0)
    return (len(values) - 1) / len(values) * centered.T @ centered


class P337DirectBirthScalingArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.score = json.loads((RESULT / "score.json").read_text())
        with (RESULT / "batch_direct_rates.csv").open(newline="") as handle:
            cls.rows = list(csv.DictReader(handle))

    def test_sufficient_statistics_hash_and_size(self) -> None:
        path = RESULT / "batch_direct_rates.csv"
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        self.assertEqual(digest, self.score["batch_sufficient_statistics"]["sha256"])
        self.assertEqual(len(self.rows), 260)

    def test_points_and_paired_covariances_reconstruct(self) -> None:
        for name in self.score["size_order"]:
            expected = self.score["sizes"][name]
            rows = [row for row in self.rows if row["size"] == name]
            direct = np.asarray([
                [float(row["first_direct_count"]), float(row["second_direct_count"])]
                for row in rows
            ])
            samples = np.asarray([
                [float(row["first_samples"]), float(row["second_samples"])]
                for row in rows
            ])
            point = direct.sum(axis=0) / samples.sum(axis=0)
            deleted = np.asarray([
                (direct.sum(axis=0) - direct[index]) / (samples.sum(axis=0) - samples[index])
                for index in range(len(rows))
            ])
            np.testing.assert_allclose(
                point,
                [expected["D_by_orientation"]["first"], expected["D_by_orientation"]["second"]],
            )
            np.testing.assert_allclose(
                jackknife_covariance(deleted), expected["orientation_delete_one_covariance"],
                rtol=1e-13, atol=1e-20,
            )

    def test_external_lineage_is_excluded_and_boundary_is_conditional(self) -> None:
        self.assertIn("excluded", self.score["external_N325_N425"]["role"])
        self.assertEqual(self.score["decision"], "conditional_fixed_5_6_line_rejected")
        self.assertLess(
            abs(
                self.score["models"]["doubling_contrasts"]["ratios"][-1]["observed_ratio"]
                - self.score["models"]["doubling_contrasts"]["ratios"][-1]["fixed_ratio"]
            ),
            0.0002,
        )


if __name__ == "__main__":
    unittest.main()
