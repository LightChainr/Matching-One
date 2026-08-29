from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p252_relative_source_observable_closure import build_oracle  # noqa: E402


class RelativeSourceObserverClosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_oracle()
        cls.observers = cls.result["tiny_exact_oracle"]["observers"]

    def test_every_observer_sector_decomposition_is_exact(self) -> None:
        self.assertTrue(
            all(row["sector_reconstruction_exact"] for row in self.observers)
        )

    def test_every_observer_third_cumulant_closes(self) -> None:
        self.assertTrue(
            all(
                row["mixed_connected_cumulants"]["closure_exact"]
                for row in self.observers
            )
        )

    def test_all_higher_raw_source_rows_repeat(self) -> None:
        for row in self.observers:
            self.assertTrue(row["odd_raw_rows_equal"])
            self.assertTrue(row["positive_even_raw_rows_equal"])

    def test_source_response_rank_is_exactly_three(self) -> None:
        tiny = self.result["tiny_exact_oracle"]
        self.assertEqual(tiny["observer_by_source_power_matrix_rank"], 3)
        self.assertEqual(tiny["rank_bound"], 3)

    def test_claim_does_not_forbid_an_independent_rank3_observer(self) -> None:
        consequence = self.result["rank3_consequence"]
        self.assertIn("may couple", consequence["allowed"])
        self.assertIn("cannot be identified merely", consequence["issue_252"])

    def test_committed_artifact_is_reproducible(self) -> None:
        committed = json.loads(
            (
                ROOT
                / "results"
                / "exact-relative-source-observer-closure"
                / "latest.json"
            ).read_text()
        )
        self.assertEqual(committed, self.result)


if __name__ == "__main__":
    unittest.main()
