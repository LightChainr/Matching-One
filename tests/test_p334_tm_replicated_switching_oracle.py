import json
from pathlib import Path
import unittest

from p334_tm_replicated_switching_oracle import build_result


ROOT = Path(__file__).resolve().parents[1]


class P334TMReplicatedSwitchingOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_one_common_site_flow_is_complete(self):
        audit = self.result["bounded_audit"]
        self.assertEqual(audit["line_layer_rows"], 984)
        self.assertEqual(audit["one_common_site_flow_pass"], 984)
        self.assertEqual(audit["one_common_site_flow_fail"], 0)
        self.assertEqual(audit["rows_using_one_endpoint_switch"], 658)

    def test_exact_pair_locking_has_minimal_bounded_failure(self):
        audit = self.result["bounded_audit"]
        self.assertEqual(audit["exact_pair_fail"], 658)
        first = self.result["exact_pair_locking"]["first_failure"]
        self.assertEqual(first["N"], 6)
        self.assertEqual(first["matrix"], [[2, 0], [0, 3]])
        self.assertEqual(first["carrier"], "matching")
        self.assertTrue(
            all(row["demand"] > row["same_pair_supply"] for row in first["deficits"])
        )

    def test_every_induced_Hall_cut_passes(self):
        audit = self.result["bounded_audit"]
        self.assertEqual(audit["Hall_families_checked"], 1_176_258)
        self.assertEqual(audit["Hall_failures"], 0)
        tight = self.result["tightest_Hall_cut"]
        self.assertEqual(tight["required_over_available"], "24/25")
        self.assertEqual(tight["row"]["active_vertices"], list(range(12)))

    def test_flow_is_an_integer_partition_of_demand(self):
        flow = self.result["flow_digest"]
        self.assertEqual(
            flow["demand_tokens"],
            flow["exact_pair_tokens"] + flow["one_endpoint_switched_tokens"],
        )
        self.assertEqual(
            flow["demand_tokens"],
            flow["synergy_tokens_used"] + flow["reservoir_tokens_used"],
        )
        self.assertEqual(
            self.result["largest_canonical_switch_fraction"]["fraction"],
            "88/91",
        )

    def test_checked_artifact_reproduces(self):
        checked = json.loads(
            (
                ROOT
                / "results/p334-tm-replicated-switching-oracle/latest.json"
            ).read_text()
        )
        self.assertEqual(checked, self.result)


if __name__ == "__main__":
    unittest.main()
