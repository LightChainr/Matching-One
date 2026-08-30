from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_marked_birth_path as base  # noqa: E402
from score_p267_density_clock_orthogonal import (  # noqa: E402
    METRICS,
    fixed_k_mu,
    subtract_group,
)


class P267DensityClockTests(unittest.TestCase):
    def test_fixed_k_mu_matches_exact_three_by_three_torus_census(self) -> None:
        side = 3
        n = side * side
        edges = set()
        faces = []
        index = lambda x, y: (x % side) + side * (y % side)
        for y in range(side):
            for x in range(side):
                v = index(x, y)
                for w in (index(x + 1, y), index(x, y + 1)):
                    edges.add(tuple(sorted((v, w))))
                faces.append({v, index(x + 1, y), index(x, y + 1), index(x + 1, y + 1)})
        by_k = {k: [] for k in range(n + 1)}
        for mask in range(1 << n):
            occupied = {v for v in range(n) if mask & (1 << v)}
            value = len(occupied)
            value -= sum(edge[0] in occupied and edge[1] in occupied for edge in edges)
            value += sum(face <= occupied for face in faces)
            by_k[len(occupied)].append(value)
        for k, values in by_k.items():
            self.assertAlmostEqual(
                float(fixed_k_mu(n, k)), sum(values) / len(values), places=13
            )

    def test_gzip_path_reader_is_lossless(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tiny.path.csv.gz"
            with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("n", "a", "b", "orientation", "batch", "samples", "k"),
                )
                writer.writeheader()
                writer.writerow({
                    "n": 1, "a": 1, "b": 0, "orientation": "first",
                    "batch": 0, "samples": 2, "k": 0,
                })
            groups = base.read_path(path)
            self.assertEqual(list(groups), [(1, "first", 0)])
            self.assertEqual(groups[(1, "first", 0)][0].samples, 2)

    def test_total_minus_omitted_preserves_exact_path_aggregates(self) -> None:
        def row(samples: int, value: int) -> base.PathRow:
            values = {name: 0 for name in base.VALUE_COLUMNS}
            values["sum_O_ext"] = value
            return base.PathRow(5, 2, 1, "first", -1, samples, 0, values)

        result = subtract_group([row(30, 17)], [row(10, 4)])
        self.assertEqual(result[0].samples, 20)
        self.assertEqual(result[0].values["sum_O_ext"], 13)

    def test_manifest_forbids_far_observer_inference(self) -> None:
        manifest = json.loads(
            (ROOT / "analysis/p267_density_clock_orthogonal_20260830.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            manifest["status"],
            "retrospective_protocol_locked_before_density_clock_score",
        )
        self.assertIn("infer or score an O_far conditional mean", manifest["forbidden"])
        self.assertEqual(METRICS[:4], ("P4_raw_re", "P4_raw_im", "P4_clock_re", "P4_clock_im"))


if __name__ == "__main__":
    unittest.main()
