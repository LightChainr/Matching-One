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

from analyze_p263_crossfit_control import (  # noqa: E402
    ACTIVE_INDICES,
    ANCHOR_INDEX,
    GEOMETRY_ORDER,
    conditional_category_noop,
    crossfit_bond_control,
    read_rows,
)


RESULT_ROOT = ROOT / "results" / "server-20260829" / "P263-boundary-qscore-pilot"


class P263CrossfitControlTests(unittest.TestCase):
    def test_revealed_category_conditioning_is_exact_noop(self) -> None:
        for level in ("level1_200k", "level2_500k"):
            rows = read_rows([RESULT_ROOT / "raw" / f"{level}.batches.csv"])
            result = conditional_category_noop(rows)
            self.assertLess(result["maximum_absolute_difference"], 5e-12)
            self.assertEqual(result["variance_ratio_to_primary"], 1.0)

    def test_runner_emits_minimal_exact_bond_control(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / "pilot"
            output = Path(directory) / "rows.csv"
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
            subprocess.run(
                [
                    str(executable), "--level", "1", "--samples", "8",
                    "--batches", "4", "--seed", "263", "--output", str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            with output.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertIn("sum_b", rows[0])
            for row in rows:
                self.assertGreaterEqual(int(row["sum_b"]), 0)
                self.assertLessEqual(
                    int(row["sum_b"]), int(row["samples"]) * int(row["edges"])
                )

    def test_exact_zero_crossfit_can_reduce_synthetic_batch_variance(self) -> None:
        rows: list[dict[str, int | str]] = []
        batch_count = 20
        for batch in range(batch_count):
            controls = [
                (batch % 5) - 2,
                ((2 * batch) % 5) - 2,
                ((3 * batch) % 5) - 2,
                ((4 * batch) % 5) - 2,
            ]
            noises = [
                ((batch + index) % 3) - 1 for index in range(len(GEOMETRY_ORDER))
            ]
            for index, geometry in enumerate(GEOMETRY_ORDER):
                contribution = 3 * controls[index] + noises[index]
                rows.append(
                    {
                        "geometry_id": geometry,
                        "batch": batch,
                        "samples": 100,
                        "edges": 100,
                        "sum_b": 5000 + 50 * controls[index],
                        "sum_J": 10000,
                        "count_14_23": 10,
                        "sum_J_14_23": 1000 + 20 * contribution,
                    }
                )
        dlog = []
        for geometry in GEOMETRY_ORDER:
            selected = [row for row in rows if row["geometry_id"] == geometry]
            dlog.append(
                sum(int(row["sum_J_14_23"]) for row in selected)
                / (2 * sum(int(row["count_14_23"]) for row in selected))
                - sum(int(row["sum_J"]) for row in selected)
                / (2 * sum(int(row["samples"]) for row in selected))
            )
        residual = [dlog[index] - dlog[ANCHOR_INDEX] for index in ACTIVE_INDICES]
        score = {
            "residual": residual,
            "residual_covariance": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "joint_gls": {"chi_square": sum(value * value for value in residual)},
        }
        result = crossfit_bond_control(rows, score)
        self.assertLess(
            result["variance_comparison"]["trace_ratio_control_over_raw"], 0.5
        )

    def test_committed_primary_scores_remain_the_reported_comparison(self) -> None:
        expected = {
            "level1_200k": 1.0977457461886801,
            "level2_500k": 5.872171777302203,
        }
        for level, chi_square in expected.items():
            payload = json.loads(
                (RESULT_ROOT / "analysis" / f"{level}.score.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(payload["joint_gls"]["chi_square"], chi_square)


if __name__ == "__main__":
    unittest.main()
