from __future__ import annotations

import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "global_matching_spin4_selection.py"
SPEC = importlib.util.spec_from_file_location("global_matching_selection", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class GlobalMatchingSpin4SelectionTests(unittest.TestCase):
    def test_invariant_covector_annihilates_non_singlets(self) -> None:
        for q in (4, 5, 6, 7):
            projectors = MODULE.P262.unordered_pair_projectors(q)
            invariant = [Fraction(1)] * (q * (q - 1) // 2)
            self.assertTrue(all(x == 0 for x in MODULE.row_times_matrix(invariant, projectors["standard"])))
            self.assertTrue(all(x == 0 for x in MODULE.row_times_matrix(invariant, projectors["two_row_2"])))
            self.assertEqual(MODULE.row_times_matrix(invariant, projectors["singlet"]), invariant)

    def test_charged_vector_is_invisible_linearly_but_not_quadratically(self) -> None:
        row = MODULE.integer_selection_check(5)["charged_positive_control"]
        self.assertEqual(row["invariant_one_point"]["text"], "0")
        self.assertGreater(Fraction(row["two_point_norm"]["text"]), 0)

    def test_s4_transposition_is_one_versus_zero_discriminator(self) -> None:
        rows = MODULE.s4_transposition_oracle()
        self.assertEqual(rows["singlet"]["twist_to_identity_trace_ratio"]["text"], "1")
        self.assertEqual(rows["two_row_2"]["twist_to_identity_trace_ratio"]["text"], "0")
        self.assertEqual(rows["standard"]["transposition_character"]["text"], "1")

    def test_claim_is_explicitly_not_an_arbitrary_trace_exclusion(self) -> None:
        payload = MODULE.analyze()
        self.assertEqual(payload["claim_boundary"]["categorical_trace_at_Q1"]["text"], "-1")
        self.assertIn("linear matrix-element null", payload["claim_boundary"]["minimal_counterexample"])

    def test_frozen_artifact_is_reproducible(self) -> None:
        expected = json.loads((ROOT / "predictions" / "p257_global_matching_spin4_selection_20260829.json").read_text())
        self.assertEqual(MODULE.analyze(), expected)


if __name__ == "__main__":
    unittest.main()
