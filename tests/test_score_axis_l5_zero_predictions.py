#!/usr/bin/env python3
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from score_axis_l5_zero_predictions import score  # noqa: E402


class ScoreAxisL5ZeroPredictionsTests(unittest.TestCase):
    def test_rejects_non_l5_artifact(self) -> None:
        with self.assertRaisesRegex(ValueError, "axis L=5"):
            score({"geometry": "axis", "L": 4}, {}, 70)

    def test_rejects_changed_target(self) -> None:
        pilot = {"predictions": {"axis": {"prospective_next_size": [
            {"metric": name, "target_L": 6, "target_N": 36, "prediction": "1",
             "model": "frozen", "training_L": [3, 4], "training_N": [9, 16]}
            for name in ("physical_root_0_1", "imaginary_rms", "nonreal_fraction")
        ]}}}
        with self.assertRaisesRegex(ValueError, "not axis L=5"):
            score({"geometry": "axis", "L": 5}, pilot, 70)


if __name__ == "__main__":
    unittest.main()
