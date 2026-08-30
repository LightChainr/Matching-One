from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/local-20260830/P334-age-iota-site-control"


def jackknife_covariance(values: np.ndarray) -> np.ndarray:
    centered = values - values.mean(axis=0)
    return (len(values) - 1) / len(values) * centered.T @ centered


class P334AgeIotaSiteArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.score = json.loads((RESULT / "score.json").read_text())
        with (RESULT / "batch_delete_one_slopes.csv").open(newline="") as handle:
            cls.rows = list(csv.DictReader(handle))

    def test_batch_hash_and_rows(self) -> None:
        path = RESULT / "batch_delete_one_slopes.csv"
        self.assertEqual(
            hashlib.sha256(path.read_bytes()).hexdigest(),
            self.score["batch_delete_one_slopes"]["sha256"],
        )
        self.assertEqual(len(self.rows), 200)

    def test_covariances_reconstruct(self) -> None:
        for name, size in self.score["sizes"].items():
            rows = [row for row in self.rows if int(row["N"]) == size["N"]]
            order = size["point_vector_order"]
            values = np.asarray([[float(row[key]) for key in order] for row in rows])
            np.testing.assert_allclose(
                jackknife_covariance(values), size["delete_one_covariance"],
                rtol=1e-13, atol=1e-18,
            )

    def test_iota_is_exact_noop_and_site_control_survives(self) -> None:
        self.assertTrue(self.score["primary_replay"]["passed"])
        for size in self.score["sizes"].values():
            self.assertTrue(size["iota_exactly_saturated"])
            self.assertEqual(
                {(row["iota01"], row["iota12"]) for row in size["iota_pair_support"]},
                {(1, 1)},
            )
            self.assertTrue(size["controls"]["smith_site_pair"]["association_survives"])
            self.assertGreater(
                size["controls"]["smith_site_pair"]["minimum_age_denominator_retention"],
                0.999,
            )


if __name__ == "__main__":
    unittest.main()
