from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p144_surface_polynomial_classification as classification  # noqa: E402


class P144SurfacePolynomialClassificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = classification.build_artifact()

    def test_one_edge_and_transition_states_miss_J4(self) -> None:
        no_go = self.artifact["local_partition_no_go"]
        self.assertFalse(no_go["ordinary_single_edge_contains_J4"])
        self.assertFalse(no_go["transition_pairings_contain_J4"])
        self.assertEqual(no_go["minimum_independent_edges_needed_for_J4"], 3)
        self.assertEqual(no_go["independent_subsets_of_minimum_gadget"], 8)
        self.assertEqual(
            no_go["unwanted_mixed_subsets_if_only_all_off_and_all_on_are_kept"], 6
        )

    def test_named_families_fail_state_index_not_only_variable_choice(self) -> None:
        rows = self.artifact["named_family_classification"]
        self.assertEqual(len(rows), 5)
        self.assertTrue(all(not row["direct_site_specialization"] for row in rows))
        krushkal = next(row for row in rows if row["family"] == "Krushkal")
        self.assertIn("sufficient for q", krushkal["topological_readout"])
        self.assertIn("state", krushkal["reason"])

    def test_rank_relative_and_genus_derivatives_agree_exactly(self) -> None:
        quotient = self.artifact["exact_site_rank_image_quotient"]
        expected = [-1, -9, -36, -78, -90, -36, 36, 36, 9, 1]
        self.assertEqual(quotient["configurations"], 512)
        self.assertEqual(quotient["matching_Bernstein_coefficients"], expected)
        self.assertEqual(quotient["matching_from_black_rank_derivative"], expected)
        self.assertEqual(quotient["matching_from_relative_rank_derivative"], expected)
        self.assertEqual(quotient["matching_from_Krushkal_genus_derivative"], expected)
        self.assertTrue(quotient["all_derivative_routes_agree"])

    def test_unrestricted_rank_sum_reduces_terminal_topology_to_one_source(self) -> None:
        quotient = self.artifact["exact_site_rank_image_quotient"]
        self.assertEqual(quotient["rank_pairs"], [[0, 2], [1, 1], [2, 0]])
        self.assertEqual(quotient["carrier_genus_pairs"], [[0, 0], [0, 1], [1, 0]])
        sources = quotient["independent_output_sources"]
        self.assertEqual(sources["extra_topology_sources_beyond_occupation"], 1)
        self.assertEqual(sources["fixed_N_reduced_form"], ["t=a/b", "Q_relative_topology"])

    def test_order_two_is_the_minimal_topology_source_witness(self) -> None:
        witness = self.artifact["minimal_topology_source_witness"]
        self.assertEqual(witness["first_distinguishing_order"], 2)
        row = witness["witness"]
        self.assertEqual(row["shared_occupation_only_counts"], [1, 2, 1])
        self.assertNotEqual(row["matching_A"], row["matching_B"])
        self.assertEqual(witness["audits"][0]["distinct_matching_polynomials"], 1)
        self.assertGreater(witness["audits"][1]["distinct_matching_polynomials"], 1)

    def test_checked_artifact_reproduces(self) -> None:
        checked = json.loads(
            (
                ROOT
                / "results"
                / "p144-surface-polynomial-classification"
                / "latest.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(checked, self.artifact)
        self.assertTrue(self.artifact["all_machine_gates_pass"])


if __name__ == "__main__":
    unittest.main()
