import json
from pathlib import Path
import unittest

from p334_linewise_marked_path_oracle import build_result


ROOT = Path(__file__).resolve().parents[1]


class P334LinewiseMarkedPathOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_all_linewise_path_conditions_pass(self):
        audit = self.result["bounded_audit"]
        self.assertEqual(audit["line_layer_rows"], 984)
        self.assertEqual(audit["BA_pass"], 984)
        self.assertEqual(audit["TM_pass"], 984)
        self.assertEqual(audit["square_switch_failures"], 0)
        self.assertEqual(audit["complement_rank_one_states"], 59922)
        self.assertEqual(audit["complement_state_failures"], 0)
        self.assertEqual(audit["complement_degree_swap_failures"], 0)

    def test_BA_has_nontrivial_concordance_equalities(self):
        ba = self.result["BA_concordance_extreme"]
        self.assertEqual(ba["maximum_discordance_over_concordance"], "1")
        self.assertEqual(ba["nonzero_equality_count"], 8)
        self.assertTrue(
            all(
                row["BA"]["concordance_mass"]
                == row["BA"]["discordance_mass"]
                == 1600
                for row in ba["nonzero_equalities"]
            )
        )

    def test_TM_near_extreme_requires_independent_reservoir(self):
        tm = self.result["TM_path_extreme"]
        self.assertEqual(tm["maximum_demand_over_supply"], "24/25")
        self.assertEqual(tm["maximizer_count"], 6)
        for row in tm["maximizers"]:
            self.assertEqual(row["N"], 12)
            self.assertEqual(row["TM"]["new_exit_triples"], 0)
            self.assertEqual(row["TM"]["supply"], 43200)
            self.assertEqual(row["TM"]["demand"], 41472)
        self.assertEqual(
            self.result["bounded_audit"][
                "TM_synergy_absent_but_demand_positive"
            ],
            74,
        )

    def test_marked_path_algebra_is_checked_in_extremes(self):
        rows = self.result["TM_path_extreme"]["maximizers"]
        for row in rows:
            tm = row["TM"]
            self.assertEqual(
                tm["supply"], tm["supply_synergy"] + tm["supply_independent"]
            )
            self.assertEqual(tm["margin"], tm["supply"] - tm["demand"])

    def test_checked_in_artifact_reproduces(self):
        checked = json.loads(
            (ROOT / "results/p334-linewise-marked-path-oracle/latest.json").read_text()
        )
        self.assertEqual(checked, self.result)


if __name__ == "__main__":
    unittest.main()
