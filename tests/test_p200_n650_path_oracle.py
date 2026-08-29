from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p200_n650_path_oracle import render  # noqa: E402


class P200N650PathOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = render()

    def test_both_lineages_have_exact_same_endpoint_by_path(self) -> None:
        endpoints = []
        for lineage in self.payload["n650_lineages"]:
            self.assertTrue(lineage["same_integer_period_graph_not_merely_isomorphic"])
            self.assertEqual(lineage["final_smith_invariants"], [1, 650])
            endpoints.append(lineage["common_final_gaussian"])
        self.assertEqual(endpoints, [[23, 11], [17, 19]])

    def test_crt_and_configurationwise_character_commutator(self) -> None:
        oracle = self.payload["single_parent_fiber_exhaustive_oracle"]
        self.assertTrue(oracle["crt_is_bijection"])
        self.assertEqual(len(set(map(tuple, oracle["crt_pairs_in_full_representative_order"]))), 10)
        self.assertEqual(oracle["sequential_character_transform_nonzero_commutators"], 0)
        direct = oracle["functorial_boolean_direct_image"]
        self.assertEqual(direct["OR_nonzero_order_commutators"], 0)
        self.assertEqual(direct["AND_nonzero_order_commutators"], 0)
        mixed = oracle["symmetric_mixed_partition_witness"]
        self.assertTrue(mixed["join_orders_equal"])
        self.assertEqual(mixed["ranks"]["h_Pi_join_R2"], 5)
        self.assertEqual(mixed["ranks"]["h_Pi_join_R5"], 8)
        self.assertEqual(mixed["Delta25_h"], -4)

    def test_generic_binary_mask_does_not_descend(self) -> None:
        counts = self.payload["single_parent_fiber_exhaustive_oracle"]["descent_counts"]
        self.assertEqual(counts["to_Z5_intermediate"], 32)
        self.assertEqual(counts["to_Z2_intermediate"], 4)
        self.assertEqual(counts["to_both"], 2)
        self.assertIn(
            "undefined as a Bernoulli intermediate mask",
            self.payload["exact_conclusions"]["intermediate_homology_flag_without_pushdown_rule"],
        )

    def test_frozen_artifact_matches(self) -> None:
        artifact = json.loads(
            (ROOT / "results" / "exact-p200-n650-path" / "latest.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(artifact, self.payload)


if __name__ == "__main__":
    unittest.main()
