from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts"))

import run_finite_size_grid as grid  # noqa: E402
import summarize_finite_size_grid as summary  # noqa: E402


class GridRunnerTests(unittest.TestCase):
    def test_default_grid_has_all_54_unique_jobs(self) -> None:
        jobs = grid.build_jobs((60, 100, 160), range(5, 11), (2, 3, 4))
        self.assertEqual(len(jobs), 54)
        self.assertEqual(len({job.name for job in jobs}), 54)
        self.assertEqual({job.dps for job in jobs}, {60, 100, 160})
        self.assertEqual({job.min_train for job in jobs}, set(range(5, 11)))
        self.assertEqual({job.holdout for job in jobs}, {2, 3, 4})

    def test_resume_requires_matching_metadata_and_causal_folds(self) -> None:
        job = grid.Job(dps=100, min_train=8, holdout=2)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "audit.json"
            payload = {
                "dps": 100,
                "min_train": 8,
                "holdout": 2,
                "summaries": [{"model": "4"}],
                "folds": [{"train_max": 18, "test_min": 19}],
            }
            path.write_text(json.dumps(payload), encoding="utf-8")
            self.assertTrue(grid.valid_output(path, job))
            payload["folds"][0] = {"train_max": 19, "test_min": 19}
            path.write_text(json.dumps(payload), encoding="utf-8")
            self.assertFalse(grid.valid_output(path, job))


class TrainingOnlySummaryTests(unittest.TestCase):
    @staticmethod
    def fold(model: str, train_max: int, test_min: int, rmse: str, intercept: str) -> dict:
        return {
            "model": model,
            "n_min": 1,
            "train_max": train_max,
            "test_min": test_min,
            "test_max": test_min,
            "intercept": intercept,
            "rmse": rmse,
            "max_abs": rmse,
        }

    def test_final_tail_cannot_affect_selection(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            csv_path = base / "sequence.csv"
            csv_path.write_text(
                "n,value\n"
                + "".join(
                    f"{n},{0.5 + 2.0 / n**4:.17g}\n" for n in range(1, 11)
                ),
                encoding="utf-8",
            )
            raw_dir = base / "raw"
            raw_dir.mkdir()
            payload = {
                "dps": 160,
                "min_train": 1,
                "holdout": 1,
                # These deliberately misleading full-data summaries must be ignored.
                "summaries": [
                    {"model": "4", "score": "100"},
                    {"model": "4,6", "score": "-100"},
                ],
                "folds": [
                    self.fold("4", 4, 5, "1e-12", "0.500000000000"),
                    self.fold("4", 5, 6, "2e-12", "0.500000000001"),
                    self.fold("4", 8, 9, "1e3", "2.0"),
                    self.fold("4,6", 4, 5, "1e-5", "0.50"),
                    self.fold("4,6", 5, 6, "2e-5", "0.51"),
                    self.fold("4,6", 8, 9, "1e-100", "0.500000"),
                ],
            }
            (raw_dir / "run.json").write_text(json.dumps(payload), encoding="utf-8")

            result = summary.summarize(
                csv_path=csv_path,
                raw_dir=raw_dir,
                final_tail=2,
                selection_dps=None,
                min_validation_folds=2,
            )

            self.assertEqual(result["knowledge_cutoff"], 8)
            self.assertEqual(result["selected"]["model"], "4")
            self.assertLessEqual(result["selected"]["validation_test_max"], 8)
            self.assertEqual(result["final_tail_score"]["training_n_max"], 8)
            self.assertEqual(result["withheld_widths"], [9, 10])


if __name__ == "__main__":
    unittest.main()
