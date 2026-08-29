from __future__ import annotations

import csv
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p263_local_stopped_score import render  # noqa: E402
from score_p263_local_stopped_qscore import read_rows, score  # noqa: E402


GEOMETRIES = (
    ("lambda_1_4", 1, 4, 15, 10950),
    ("lambda_1_3", 1, 3, 14, 9548),
    ("lambda_2_3", 2, 3, 15, 10950),
    ("lambda_3_4", 3, 4, 14, 9548),
)


class P263LocalStoppedQScoreTests(unittest.TestCase):
    def test_committed_exact_oracle_recomputes(self) -> None:
        path = (
            ROOT / "results" / "post-reveal-20260829"
            / "P263-local-stopped-qscore" / "exact.json"
        )
        self.assertEqual(render(), json.loads(path.read_text(encoding="utf-8")))

    def test_exact_coupled_identity_and_far_field_variance_cancellation(self) -> None:
        rows = render()["tiny_graph"]["rows"]
        self.assertEqual({row["target_covariance"] for row in rows}, {"-1/256"})
        for spectators, row in enumerate(rows):
            self.assertEqual(
                Fraction(row["global_centered_variance"]),
                Fraction(15 + 256 * spectators, 65536),
            )
            self.assertEqual(Fraction(row["ideal_stopped_variance"]), Fraction(15, 65536))
            self.assertEqual(
                Fraction(row["completion_noise_at_one_inner_draw"]),
                Fraction(47, 4096),
            )
        self.assertGreater(rows[0]["coupled_inner_replicates"]["8"]["ratio_to_global"], 1)
        self.assertLess(rows[4]["coupled_inner_replicates"]["8"]["ratio_to_global"], 0.11)

    def test_runner_is_deterministic_and_emits_local_sufficient_statistics(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / "pilot"
            first = Path(directory) / "first.csv"
            second = Path(directory) / "second.csv"
            subprocess.run(
                [
                    "c++", "-O1", "-std=c++17",
                    str(ROOT / "src" / "p263_local_stopped_qscore_pilot.cpp"),
                    "-o", str(executable),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            command = [
                str(executable), "--level", "1", "--samples", "8",
                "--batches", "4", "--inner", "2", "--outer-seed", "2026102633",
                "--completion-seed", "2026102634", "--output",
            ]
            subprocess.run(command + [str(first)], check=True, capture_output=True, text=True)
            subprocess.run(command + [str(second)], check=True, capture_output=True, text=True)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            with first.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 16)
            for row in rows:
                samples = int(row["samples"])
                inner = int(row["inner_replicates"])
                self.assertEqual(samples, 2)
                self.assertNotEqual(row["outer_seed"], row["completion_seed"])
                self.assertLessEqual(
                    abs(int(row["sum_delta_J_14_23"])),
                    inner * int(row["sum_revealed_14_23"]),
                )
                self.assertLessEqual(
                    int(row["sum_delta_J2_individual_14_23"]),
                    inner * int(row["sum_revealed2"]),
                )

    def test_scorer_recomputes_local_ratio_and_completion_noise(self) -> None:
        rows: list[dict[str, int | str]] = []
        for batch in range(5):
            for index, (geometry, numerator, denominator, span, edges) in enumerate(GEOMETRIES):
                delta = 40 + 4 * index + 2 * batch + index * batch
                rows.append({
                    "geometry_id": geometry,
                    "lambda_num": numerator,
                    "lambda_den": denominator,
                    "level": 1,
                    "span_L": span,
                    "nx": 1,
                    "ny": 1,
                    "vertices": 1,
                    "edges": edges,
                    "batch": batch,
                    "outer_seed": 2026102633,
                    "completion_seed": 2026102634,
                    "inner_replicates": 2,
                    "sample_begin": 100 * batch,
                    "samples": 100,
                    "count_14_23": 10,
                    "sum_delta_J_14_23": delta,
                    "sum_delta_J_inner_square_14_23": 2000 + delta * delta,
                    "sum_delta_J2_individual_14_23": 3000 + delta * delta,
                    "sum_revealed": 5000,
                    "sum_revealed2": 300000,
                    "sum_revealed_14_23": 1000,
                    "max_revealed": 80,
                })
        payload = score(rows)
        first = payload["estimates"][0]
        self.assertAlmostEqual(first["measure_tangent_14_23"], 220 / 2000)
        self.assertAlmostEqual(first["d_log_probability"], 220 / 200)
        self.assertEqual(payload["batch_count"], 5)
        self.assertGreater(payload["covariance_trace"], 0)
        self.assertIsNotNone(first["estimated_completion_noise_in_outer_variance"])

    def test_committed_frozen_pilot_score_recomputes(self) -> None:
        directory = (
            ROOT / "results" / "post-reveal-20260829"
            / "P263-local-stopped-qscore" / "pilot"
        )
        actual = score(read_rows([directory / "level1_20k.batches.csv"]))
        frozen = json.loads(
            (directory / "level1_20k.score.json").read_text(encoding="utf-8")
        )
        self.assertEqual(actual, frozen)
        self.assertAlmostEqual(actual["joint_gls"]["chi_square"], 6.586542432499845)
        self.assertEqual(
            [row["events_14_23"] for row in actual["estimates"]],
            [12, 21, 113, 209],
        )


if __name__ == "__main__":
    unittest.main()
