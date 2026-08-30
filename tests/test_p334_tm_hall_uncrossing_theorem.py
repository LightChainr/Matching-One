import json
from pathlib import Path
import unittest

from p334_tm_hall_uncrossing_theorem import build_result


ROOT = Path(__file__).resolve().parents[1]


class P334TMHallUncrossingTheoremTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_translation_makes_both_pair_graphs_regular(self):
        audit = self.result["bounded_audit"]
        self.assertEqual(audit["line_layer_rows"], 984)
        self.assertEqual(audit["demand_regular_rows"], 984)
        self.assertEqual(audit["supply_regular_rows"], 984)

    def test_regular_cut_decomposition_closes_every_Hall_cut(self):
        audit = self.result["bounded_audit"]
        self.assertEqual(audit["nonempty_site_cuts_checked"], 2_470_440)
        self.assertEqual(audit["regular_cut_decomposition_failures"], 0)
        self.assertEqual(audit["negative_Hall_cuts"], 0)
        self.assertEqual(audit["ratio_bound_failures"], 0)

    def test_all_site_is_unique_ratio_maximizer(self):
        audit = self.result["bounded_audit"]
        self.assertEqual(audit["positive_demand_rows"], 688)
        self.assertEqual(audit["zero_demand_rows"], 296)
        self.assertEqual(audit["unique_all_site_ratio_max_rows"], 688)
        self.assertEqual(audit["nonunique_or_proper_ratio_max_rows"], 0)

    def test_tight_and_near_tight_cut_classification(self):
        cuts = self.result["tight_cut_classification"]
        self.assertEqual(cuts["largest_ratio"], "24/25")
        self.assertEqual(cuts["maximizer_count"], 6)
        self.assertEqual(cuts["largest_proper_ratio"], "576/715")
        self.assertTrue(
            all(
                row["class"] == "all_site"
                for row in cuts["near_tight_ratio_at_least_9_over_10"]
            )
        )

    def test_uncrossing_and_checked_artifact(self):
        self.assertEqual(
            self.result["uncrossing_lemma"]["truth_table_failures"], 0
        )
        checked = json.loads(
            (
                ROOT / "results/p334-tm-hall-uncrossing-theorem/latest.json"
            ).read_text()
        )
        self.assertEqual(checked, self.result)


if __name__ == "__main__":
    unittest.main()
