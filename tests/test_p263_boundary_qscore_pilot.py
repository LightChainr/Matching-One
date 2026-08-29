from __future__ import annotations

import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p263_boundary_qscore_pilot import read_rows, score  # noqa: E402


class BoundaryQScorePilotTests(unittest.TestCase):
    def test_runner_is_deterministic_and_emits_integer_sufficient_statistics(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / "pilot"
            first = Path(directory) / "first.csv"
            second = Path(directory) / "second.csv"
            subprocess.run(
                [
                    "c++", "-O1", "-std=c++17",
                    str(ROOT / "src" / "p263_boundary_qscore_pilot.cpp"),
                    "-o", str(executable),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            command = [
                str(executable), "--level", "1", "--samples", "8",
                "--batches", "4", "--seed", "263", "--output",
            ]
            subprocess.run(command + [str(first)], check=True, capture_output=True, text=True)
            subprocess.run(command + [str(second)], check=True, capture_output=True, text=True)
            self.assertEqual(first.read_bytes(), second.read_bytes())

            with first.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 16)
            self.assertEqual(
                {row["geometry_id"] for row in rows},
                {"lambda_1_4", "lambda_1_3", "lambda_2_3", "lambda_3_4"},
            )
            for row in rows:
                self.assertEqual(int(row["samples"]), 2)
                self.assertGreater(int(row["vertices"]), 0)
                self.assertGreater(int(row["edges"]), 0)
                self.assertLessEqual(
                    sum(int(row[f"count_{name}"]) for name in ("1234", "12_34", "14_23")),
                    int(row["samples"]),
                )

    def test_frozen_smoke_score_recomputes(self) -> None:
        directory = ROOT / "results" / "local-20260829" / "P263-boundary-qscore-smoke"
        actual = score(read_rows([directory / "level1_20k.batches.csv"]))
        frozen = json.loads((directory / "score.json").read_text(encoding="utf-8"))
        self.assertEqual(actual, frozen)
        self.assertEqual(actual["batch_count"], 20)
        self.assertEqual(actual["joint_gls"]["degrees_of_freedom"], 3)
        self.assertEqual(
            [row["events_14_23"] for row in actual["estimates"]],
            [19, 19, 133, 210],
        )


if __name__ == "__main__":
    unittest.main()
