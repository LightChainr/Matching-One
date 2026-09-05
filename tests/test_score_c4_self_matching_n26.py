
from __future__ import annotations
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_c4_self_matching_n26 import score  # noqa: E402


RESULT = ROOT / "results" / "local-20260828" / "P115-n26-self-matching-exact"


class ScoreC4SelfMatchingN26Tests(unittest.TestCase):
    def test_scores_frozen_targets_before_structure(self) -> None:
        output = score(
            ROOT / "predictions" / "c4_self_matching_n26_beta_targets_20260828.json",
            RESULT / "raw" / "n26_threads10.json",
            RESULT / "raw" / "n26_threads1.json",
        )
        self.assertTrue(output["independent_reproduction_identical"])
        self.assertEqual(
            output["protocol_conclusion"],
            "BOTH_FROZEN_LAWS_FAILED_STOP_NO_GENERALIZED_BETA_FIT",
        )
        self.assertFalse(output["generalized_beta_fit_performed"])
        scores = output["hypothesis_scores"]
        self.assertEqual(
            [item["name"] for item in scores],
            ["geometry_shortest_support", "antipodal_orbit_majority"],
        )
        self.assertEqual(scores[0]["first_difference"], {
            "occupation_k": 5,
            "observed": -65624,
            "target": -65528,
            "observed_minus_target": -96,
        })
        self.assertEqual(scores[1]["first_difference"], {
            "occupation_k": 5,
            "observed": -65624,
            "target": -65780,
            "observed_minus_target": 156,
        })

    def test_reproduction_mismatch_is_a_hard_error(self) -> None:
        reproduction = json.loads((RESULT / "raw" / "n26_threads1.json").read_text())
        reproduction["channels"]["either"]["M_bernstein_integer_coefficients"][5] += 1
        with tempfile.TemporaryDirectory() as directory:
            altered = Path(directory) / "altered.json"
            altered.write_text(json.dumps(reproduction), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "independent enumeration"):
                score(
                    ROOT / "predictions" / "c4_self_matching_n26_beta_targets_20260828.json",
                    RESULT / "raw" / "n26_threads10.json",
                    altered,
                )


if __name__ == "__main__":
    unittest.main()
