from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import mpmath as mp


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


class HighPrecisionCsvParseTests(unittest.TestCase):
    """CSV decimals must be converted only after grid dps is applied."""

    # More than 30 non-binary-friendly decimal digits (pi). Default mpmath
    # dps=15 cannot represent this exactly, so a premature mpf() loses digits.
    HIGH_PRECISION = "0.31415926535897932384626433832795028841971693993751"
    GRID_DPS = 80

    def setUp(self) -> None:
        self._previous_dps = int(mp.mp.dps)

    def tearDown(self) -> None:
        mp.mp.dps = self._previous_dps

    def test_csv_mpf_matches_value_parsed_after_high_dps(self) -> None:
        mp.mp.dps = 15
        truncated = mp.mpf(self.HIGH_PRECISION)
        mp.mp.dps = self.GRID_DPS
        expected = mp.mpf(self.HIGH_PRECISION)
        self.assertNotEqual(truncated, expected)

        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            csv_path = base / "sequence.csv"
            lines = ["n,value\n"]
            for n in range(1, 11):
                if n >= 9:
                    lines.append(f"{n},{self.HIGH_PRECISION}\n")
                else:
                    lines.append(f"{n},{0.5 + 2.0 / n**4:.17g}\n")
            csv_path.write_text("".join(lines), encoding="utf-8")
            raw_dir = base / "raw"
            raw_dir.mkdir()
            payload = {
                "dps": self.GRID_DPS,
                "min_train": 1,
                "holdout": 1,
                "summaries": [{"model": "4", "score": "100"}],
                "folds": [
                    TrainingOnlySummaryTests.fold("4", 4, 5, "1e-12", "0.500000000000"),
                    TrainingOnlySummaryTests.fold("4", 5, 6, "2e-12", "0.500000000001"),
                    TrainingOnlySummaryTests.fold("4", 8, 9, "1e3", "2.0"),
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

        self.assertEqual(result["selection_dps"], self.GRID_DPS)
        self.assertEqual(int(mp.mp.dps), self.GRID_DPS)
        actual = result["final_tail_score"]["predictions"][-1]["actual"]
        self.assertEqual(result["final_tail_score"]["predictions"][-1]["n"], 10)

        # Reconstruct the reference mpf only after the high working precision
        # is in force, as required by the merge-blocker regression.
        mp.mp.dps = self.GRID_DPS
        expected_after = mp.mpf(self.HIGH_PRECISION)
        self.assertEqual(expected_after, expected)
        self.assertEqual(actual, mp.nstr(expected_after, 35, strip_zeros=False))
        self.assertNotEqual(actual, mp.nstr(truncated, 35, strip_zeros=False))
        self.assertEqual(
            summary.observations_from_decimals([(10, self.HIGH_PRECISION)])[0].value,
            expected_after,
        )


if __name__ == "__main__":
    unittest.main()
