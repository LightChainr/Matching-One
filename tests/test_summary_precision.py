from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import mpmath as mp


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import summarize_finite_size_grid as summary  # noqa: E402


class SummaryPrecisionTests(unittest.TestCase):
    def test_selected_precision_is_applied_before_csv_parsing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            exact_tail = "0.123456789012345678901234567890123456789"
            csv_path = base / "sequence.csv"
            csv_path.write_text(
                "n,value\n"
                "1,0.500000000000000000000000000000000000001\n"
                "2,0.400000000000000000000000000000000000002\n"
                "3,0.300000000000000000000000000000000000003\n"
                f"4,{exact_tail}\n",
                encoding="utf-8",
            )

            raw_dir = base / "raw"
            raw_dir.mkdir()
            payload = {
                "dps": 80,
                "min_train": 1,
                "holdout": 1,
                "folds": [
                    {
                        "model": "4",
                        "n_min": 1,
                        "train_max": 1,
                        "test_min": 2,
                        "test_max": 2,
                        "intercept": "0.5",
                        "rmse": "1e-6",
                        "max_abs": "1e-6",
                    },
                    {
                        "model": "4",
                        "n_min": 1,
                        "train_max": 2,
                        "test_min": 3,
                        "test_max": 3,
                        "intercept": "0.5000000001",
                        "rmse": "2e-6",
                        "max_abs": "2e-6",
                    },
                ],
            }
            (raw_dir / "run.json").write_text(
                json.dumps(payload), encoding="utf-8"
            )

            previous_dps = mp.mp.dps
            mp.mp.dps = 15
            try:
                result = summary.summarize(
                    csv_path=csv_path,
                    raw_dir=raw_dir,
                    final_tail=1,
                    selection_dps=80,
                    min_validation_folds=2,
                )
                with mp.workdps(80):
                    expected = mp.nstr(
                        mp.mpf(exact_tail), 35, strip_zeros=False
                    )
                actual = result["final_tail_score"]["predictions"][0]["actual"]
                self.assertEqual(actual, expected)
            finally:
                mp.mp.dps = previous_dps


if __name__ == "__main__":
    unittest.main()
