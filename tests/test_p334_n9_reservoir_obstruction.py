from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p334_n9_reservoir_obstruction import build_result, selected_rows  # noqa: E402


class P334N9ReservoirObstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = build_result()

    def test_strict_descriptor_selects_exactly_rows_one_and_three(self) -> None:
        self.assertEqual([index for index, _row in selected_rows()], [1, 3])
        self.assertTrue(all(self.result["same_descriptor_isomorphism"]["face_family_bijections"].values()))
        self.assertTrue(self.result["same_descriptor_isomorphism"]["translation_group_conjugated_to_itself"])

    def test_exact_all_site_minimum_cut(self) -> None:
        obstruction = self.result["obstruction"]
        self.assertEqual(
            (obstruction["source_demand"], obstruction["reachable_targets"], obstruction["deficiency"]),
            (6912, 4752, 2160),
        )
        self.assertEqual(obstruction["deficiency_fraction"], "5/16")
        self.assertEqual(obstruction["minimum_cut"]["class_count"], 768)
        self.assertEqual(obstruction["minimum_cut"]["replica_histogram"], {0: 192, 1: 192, 2: 192, 3: 192})

    def test_two_mark_fixed_base_repair_is_ordinary_and_saturates(self) -> None:
        repair = self.result["minimal_legal_repair"]
        self.assertEqual(repair["coarse_matching"], [6912, 6912])
        self.assertEqual(repair["reachable_MM_orbits"], 20736)
        self.assertEqual(repair["new_MM_orbits"], 15984)
        self.assertEqual(repair["raw_new_MM_tokens"], 143856)
        self.assertEqual(repair["regular_degree"], 216)
        for row in self.result["rows"]:
            fixed = row["two_mark_fixed_base_repair"]
            self.assertTrue(fixed["ordinary_untagged_targets_only"])
            self.assertTrue(fixed["saturates"])
            self.assertEqual(fixed["reachable_MM_orbits"], fixed["total_MM_orbits"])


if __name__ == "__main__":
    unittest.main()
