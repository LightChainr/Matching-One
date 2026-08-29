from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p144_typed_incidence_spine import build_oracle  # noqa: E402


class TypedIncidenceSpineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_oracle()
        cls.exact = cls.result["honest_torus_exact_oracle"]

    def test_all_local_face_partitions_match(self) -> None:
        local = self.result["local_16_pattern_certificate"]
        self.assertEqual(len(local["patterns"]), 16)
        self.assertTrue(local["all_16_pass"])

    def test_axis_L3_spine_matches_reference(self) -> None:
        self.assertEqual(self.exact["configurations"], 512)
        self.assertEqual(self.exact["spine_failure_count"], 0)
        self.assertEqual(self.exact["boundary_descriptor_failure_count"], 0)

    def test_state_sum_specializes_to_matching_coefficients(self) -> None:
        self.assertEqual(
            self.exact["matching_Bernstein_coefficients_direct"],
            self.exact["matching_Bernstein_coefficients_from_state_sum"],
        )

    def test_complement_exchanges_local_binary_symbols(self) -> None:
        self.assertEqual(self.exact["complement_local_symbol_failure_count"], 0)
        states = self.result["fixed_embedded_object"]
        self.assertIn("exchange B and W", states["complement"])

    def test_pairwise_smoothings_need_junction_type(self) -> None:
        obstruction = self.result["naive_two_smoothing_obstruction"]
        self.assertFalse(obstruction["either_smoothing_equals_target"])
        self.assertIn("J4", obstruction["minimal_missing_partition_type"])

    def test_committed_artifact_is_reproducible(self) -> None:
        committed = json.loads(
            (ROOT / "results" / "exact-typed-incidence-spine" / "latest.json").read_text()
        )
        self.assertEqual(committed, self.result)


if __name__ == "__main__":
    unittest.main()
