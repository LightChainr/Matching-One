import json
from pathlib import Path
import unittest

from p334_hazard_transport_bound import build_result


ROOT = Path(__file__).resolve().parents[1]


class P334HazardTransportBoundTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_variance_association_gate_covers_all_existing_pairs(self):
        counts = self.result["bounded_counts"]
        self.assertEqual(counts["adjacent_pairs"], 984)
        self.assertEqual(counts["negative_degree_bias"], 484)
        self.assertEqual(counts["variance_domination_pass"], 984)
        self.assertEqual(counts["association_nonnegative"], 984)

    def test_negative_bias_ratio_has_exact_non_saturated_extreme(self):
        extreme = self.result["negative_bias_ratio_extreme"]
        self.assertEqual(extreme["maximum"], "5/14")
        self.assertEqual(extreme["achiever_count"], 8)
        self.assertTrue(
            all(
                row["N"] == 12
                and row["carrier"] == "matching"
                and row["lower_layer"] == 6
                and row["edge_slack"] == "2/15"
                and row["degree_bias"] == "-1/21"
                for row in extreme["achievers"]
            )
        )
        equalities = self.result["uniform_hazard_equalities"]
        self.assertEqual(equalities["count"], 78)
        self.assertTrue(equalities["all_trivial_zero_slack_zero_bias"])

    def test_first_order_transport_has_exact_n12_counterexample(self):
        row = self.result["stronger_route_counterexamples"][
            "first_order_stochastic_transport"
        ]
        self.assertEqual(row["N"], 12)
        self.assertEqual(row["matrix"], [[3, 0], [0, 4]])
        self.assertEqual(row["carrier"], "matching")
        threshold = next(
            item for item in row["tail_table"] if item["threshold"] == "3/4"
        )
        self.assertEqual(threshold["uniform_lower_tail"], "4/19")
        self.assertEqual(threshold["edge_upper_tail"], "1/5")
        self.assertEqual(threshold["difference"], "-1/95")

    def test_pointwise_comonotonicity_and_cauchy_are_too_strong(self):
        routes = self.result["stronger_route_counterexamples"]
        row = routes["pointwise_birth_exit_comonotonicity"]
        self.assertEqual(row["N"], 8)
        self.assertEqual(row["left_birth_exit_degrees"], [0, 2])
        self.assertEqual(row["right_birth_exit_degrees"], [2, 1])
        self.assertEqual(row["layer_birth_exit_covariance"], "1/90")
        self.assertEqual(routes["cauchy_worst_case_failure_count"], 4)
        first = routes["cauchy_worst_case_failures"][0]
        self.assertEqual(first["cauchy_left_square"], "4/103041")
        self.assertEqual(
            first["cauchy_variance_product"], "50545/1179716409"
        )

    def test_checked_in_artifact_reproduces(self):
        checked = json.loads(
            (ROOT / "results/p334-hazard-transport-bound/latest.json").read_text()
        )
        self.assertEqual(checked, self.result)


if __name__ == "__main__":
    unittest.main()
