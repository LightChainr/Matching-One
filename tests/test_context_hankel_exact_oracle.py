from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from context_hankel_exact_oracle import build_artifact, determinant, rank  # noqa: E402


class ContextHankelExactOracleTests(unittest.TestCase):
    def test_exact_rank_helpers(self) -> None:
        matrix = [[Fraction(1), Fraction(1)], [Fraction(1), Fraction(2)]]
        self.assertEqual(rank(matrix), 2)
        self.assertEqual(determinant(matrix), 1)
        self.assertEqual(rank([[Fraction(1), Fraction(2)], [Fraction(2), Fraction(4)]]), 1)

    def test_same_endpoint_words_are_blind_in_unmarked_block(self) -> None:
        artifact = build_artifact()
        words = artifact["exact_gaussian_words"]
        endpoint = artifact["endpoint_product_sector"]
        self.assertIs(words["same_endpoint"], True)
        self.assertEqual(words["endpoint_smith"]["aa"], [2, 2])
        self.assertEqual(words["endpoint_smith"]["direct_2i"], [2, 2])
        self.assertEqual(endpoint["aa_equals_direct_column"], [True, True])
        self.assertEqual(endpoint["rank"], 2)

    def test_allowed_charged_context_raises_declared_rank_to_three(self) -> None:
        artifact = build_artifact()
        enriched = artifact["sector_resolved_enrichment"]
        self.assertEqual(enriched["rank"], 3)
        self.assertEqual(enriched["determinant"], "1")
        self.assertIn("direct-sum", enriched["interpretation"])

    def test_morphism_defect_is_rank_one_and_endpoint_invisible(self) -> None:
        artifact = build_artifact()
        witness = artifact["morphism_sensitive_exact_witness"]
        self.assertEqual(witness["defect_rank"], 1)
        self.assertEqual(
            witness["A_direct_minus_A_a_squared"],
            [["0", "0", "0"], ["0", "0", "0"], ["0", "0", "-1"]],
        )

    def test_archive_gate_forbids_completion_by_fit(self) -> None:
        artifact = build_artifact()
        gate = artifact["archive_gate"]
        self.assertIs(gate["rank_scored_from_current_archive"], False)
        self.assertIn("Do not fill", gate["do_not_do"])
        self.assertEqual(artifact["frozen_physical_prediction"]["combined_unstructured_rank"], 3)


if __name__ == "__main__":
    unittest.main()
