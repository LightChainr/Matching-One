import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_certificate_boundary_ranker import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    frozen_candidates,
    rank_candidates,
    validate_result,
)


class ExactCertificateBoundaryRankerTests(unittest.TestCase):
    def test_frozen_ranking_is_exact_and_cost_aware(self) -> None:
        result = build_result()
        self.assertEqual(result["verification"]["top_candidate"], "morphism-sensitive-row")
        self.assertEqual(result["verification"]["disjoint_count"], 2)
        self.assertEqual(result["ranking"][0]["margin_per_cost"], "1/4")

    def test_overlapping_intervals_have_zero_margin(self) -> None:
        ranking = rank_candidates([{"id": "x", "feasible_interval": ["0", "1"], "forecast_interval": ["1", "2"], "cost": "1"}])
        self.assertFalse(ranking[0]["disjoint"])
        self.assertEqual(ranking[0]["separation_margin"], "0")

    def test_nonpositive_cost_and_duplicate_ids_fail_closed(self) -> None:
        candidates = frozen_candidates()
        candidates[0]["cost"] = "0"
        with self.assertRaisesRegex(ValueError, "positive"):
            rank_candidates(candidates)
        duplicate = [copy.deepcopy(frozen_candidates()[0]), copy.deepcopy(frozen_candidates()[0])]
        with self.assertRaisesRegex(ValueError, "unique"):
            rank_candidates(duplicate)

    def test_checked_in_result_reproduces(self) -> None:
        result = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        summary = validate_result(result)
        self.assertEqual(summary["status"], "valid_exact_certificate_boundary_ranking")
        self.assertEqual(summary["candidate_count"], 3)


if __name__ == "__main__":
    unittest.main()
