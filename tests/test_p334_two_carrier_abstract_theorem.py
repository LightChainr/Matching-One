import json
from pathlib import Path
import unittest

from p334_two_carrier_abstract_theorem import build_result, dual_system, upsets


ROOT = Path(__file__).resolve().parents[1]


class P334TwoCarrierAbstractTheoremTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result(ROOT)

    def test_boolean_upset_and_system_scan_is_exhaustive(self):
        self.assertEqual([len(upsets(n)) for n in (2, 3, 4)], [6, 20, 168])
        summary = self.result["symbolic_boolean_scan"]["summary"]
        self.assertEqual(
            [row["admissible_systems"] for row in summary], [5, 111, 7076]
        )
        self.assertTrue(
            all(
                row["boundary_association_failures"] == 0
                and row["transport_moment_failures"] == 0
                and row["ulc_failures"] == 0
                for row in summary[:2]
            )
        )
        self.assertEqual(summary[2]["boundary_association_failures"], 746)
        self.assertEqual(summary[2]["transport_moment_failures"], 36)
        self.assertEqual(summary[2]["ulc_failures"], 48)

    def test_boundary_association_axiom_is_independent(self):
        witness = self.result["symbolic_boolean_scan"][
            "association_axiom_independence_witness"
        ]
        self.assertEqual(witness["N"], 4)
        self.assertEqual(witness["rank_one_sector_masks"], [1, 6, 9, 10, 12, 14])
        row = witness["witness_row"]
        self.assertEqual((row["transport_moment_left"], row["transport_moment_right"]), (6, 4))
        self.assertEqual(
            (row["boundary_association_left"], row["boundary_association_right"]),
            (32, 35),
        )
        self.assertEqual(row["xi_delta"], "-1/24")

    def test_transport_moment_axiom_is_independent(self):
        witness = self.result["symbolic_boolean_scan"][
            "transport_axiom_independence_witness"
        ]
        self.assertEqual(witness["rank_one_sector_masks"], [1, 2, 6, 10, 14])
        row = witness["witness_row"]
        self.assertEqual(
            (row["boundary_association_left"], row["boundary_association_right"]),
            (4, 4),
        )
        self.assertEqual(
            (row["transport_moment_left"], row["transport_moment_right"]),
            (12, 16),
        )
        self.assertEqual(row["xi_delta"], "-1/6")

    def test_basic_axioms_have_a_minimal_n4_ulc_counterexample(self):
        witness = self.result["symbolic_boolean_scan"]["first_ULC_failure"]
        self.assertEqual(witness["N"], 4)
        self.assertEqual(witness["rank_one_sector_masks"], [1, 5, 9, 13, 14])
        self.assertEqual(witness["failing_rows"][0]["margin"], "-1/72")

    def test_dual_system_is_exact_complement_reversal(self):
        dual_one, dual_two, dual_sector = dual_system(
            4,
            frozenset([1, 3, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15]),
            frozenset([3, 5, 7, 11, 13, 15]),
        )
        full = 15
        sector = frozenset([1, 6, 9, 10, 12, 14])
        self.assertEqual(dual_sector, frozenset(full ^ mask for mask in sector))
        self.assertTrue(dual_two <= dual_one)

    def test_checked_in_artifact_reproduces(self):
        checked = json.loads(
            (ROOT / "results/p334-two-carrier-abstract-theorem/latest.json").read_text()
        )
        self.assertEqual(checked, self.result)


if __name__ == "__main__":
    unittest.main()
