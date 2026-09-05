
from __future__ import annotations
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_rank_gap_thermal_window import (  # noqa: E402
    jackknife_se,
    pooled_statistics,
    read_run,
)


FIELDS = [
    "n", "a", "b", "orientation", "batch", "samples",
    "sum_kminus", "sum_kplus", "sum_kminus2", "sum_kplus2",
    "sum_product", "sum_gap", "sum_gap2",
]


def moment_row(n: int, orientation: str, batch: int, pairs: list[tuple[int, int]]) -> dict:
    minus = [pair[0] for pair in pairs]
    plus = [pair[1] for pair in pairs]
    gaps = [right - left for left, right in pairs]
    return {
        "n": n, "a": 2, "b": 1, "orientation": orientation,
        "batch": batch, "samples": len(pairs),
        "sum_kminus": sum(minus), "sum_kplus": sum(plus),
        "sum_kminus2": sum(value * value for value in minus),
        "sum_kplus2": sum(value * value for value in plus),
        "sum_product": sum(left * right for left, right in pairs),
        "sum_gap": sum(gaps), "sum_gap2": sum(value * value for value in gaps),
    }


class RankGapThermalWindowTests(unittest.TestCase):
    def make_run(self, mutate_gap2: bool = False):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        root = Path(directory.name)
        rows = []
        for orientation in ("first", "second"):
            for batch, gap in enumerate((1, 2, 3)):
                rows.append(moment_row(5, orientation, batch, [(1, 1 + gap), (2, 2 + gap)]))
        if mutate_gap2:
            rows[0]["sum_gap2"] += 1
        moments = root / "moments.csv"
        with moments.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        metadata = root / "metadata.json"
        metadata.write_text(json.dumps({
            "seed": 7,
            "replica_counter_first": 10,
            "replica_counter_last_exclusive": 16,
            "samples_per_pair": 6,
            "batches": 3,
            "designs": [{"N": 5, "first": [2, 1], "second": [2, 1]}],
            "git_commit": "0" * 40,
        }), encoding="utf-8")
        return moments, metadata

    def test_gap_mean_and_delete_one_se(self) -> None:
        moments, metadata = self.make_run()
        run = read_run(5, moments, metadata)
        self.assertAlmostEqual(float(pooled_statistics(run)["gap_mean"]), 2.0)
        deleted = [pooled_statistics(run, batch)["gap_mean"] for batch in range(3)]
        self.assertEqual([float(value) for value in deleted], [2.5, 2.0, 1.5])
        self.assertAlmostEqual(float(jackknife_se(deleted)), 1 / 3**0.5)

    def test_joint_shape_statistics_use_paired_moments(self) -> None:
        moments, metadata = self.make_run()
        statistics = pooled_statistics(read_run(5, moments, metadata))
        self.assertAlmostEqual(float(statistics["gap_variance"]), 2 / 3)
        self.assertGreater(float(statistics["rank_correlation"]), 0.0)
        self.assertGreater(float(statistics["gap_cv"]), 0.0)


if __name__ == "__main__":
    unittest.main()
