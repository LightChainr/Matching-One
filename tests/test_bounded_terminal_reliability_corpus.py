import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from bounded_terminal_reliability_corpus import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    selected_canonical_graphs,
    validate_result,
)


class BoundedTerminalReliabilityCorpusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_selection_matches_frozen_census(self):
        self.assertEqual(len(selected_canonical_graphs()), 27)
        self.assertEqual(self.result["selection"]["candidate_orbits"], 27)
        self.assertEqual(
            self.result["selection"]["edge_count_histogram"],
            {"4": 2, "5": 4, "6": 7, "7": 7, "8": 4, "9": 2, "10": 1},
        )

    def test_every_candidate_normalizes_by_open_edge_count(self):
        for candidate in self.result["candidates"]:
            edge_count = candidate["edge_count"]
            self.assertEqual(len(candidate["normalization_counts"]), edge_count + 1)
            self.assertEqual(sum(candidate["normalization_counts"]), candidate["configurations"])
            self.assertGreaterEqual(candidate["internal_degree"], 3)

    def test_star_fixture_and_claim_boundaries(self):
        self.assertEqual(self.result["enumeration"]["star_fixture_matches"], 1)
        self.assertFalse(self.result["conclusion"]["planarity_certified"])
        self.assertFalse(self.result["conclusion"]["self_duality_tested"])
        self.assertFalse(self.result["conclusion"]["new_percolation_bound"])

    def test_checked_in_corpus_reproduces_exactly(self):
        checked_in = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked_in, self.result)
        summary = validate_result(checked_in)
        self.assertEqual(summary["candidate_orbits"], 27)
        self.assertEqual(summary["total_configurations"], 4576)

    def test_digest_and_filter_tampering_fail_closed(self):
        digest = copy.deepcopy(self.result)
        digest["enumeration"]["corpus_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(digest)
        selection = copy.deepcopy(self.result)
        selection["selection"]["candidate_orbits"] = 26
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(selection)


if __name__ == "__main__":
    unittest.main()
