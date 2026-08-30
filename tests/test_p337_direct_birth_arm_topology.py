import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p337_direct_birth_arm_topology import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    carrier_descriptor,
)


class DirectBirthArmTopologyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_every_direct_edge_has_one_of_two_carriers(self):
        allowed = {"one_carrier_theta", "two_carrier_figure_eight"}
        for geometry in self.result["geometries"]:
            self.assertTrue(set(geometry["topology_counts"]) <= allowed)
            self.assertEqual(sum(geometry["topology_counts"].values()), geometry["direct_edges"])

    def test_reference_counts_and_arm_lower_bounds(self):
        rows = {row["N"]: row for row in self.result["geometries"]}
        self.assertEqual(rows[9]["topology_counts"], {
            "one_carrier_theta": 36,
            "two_carrier_figure_eight": 9,
        })
        self.assertEqual(rows[16]["topology_counts"]["two_carrier_figure_eight"], 336)
        self.assertEqual(rows[17]["direct_edges"], 8823)
        self.assertEqual(self.result["theorem"]["theta_occupied_arm_lower_bound"], 3)
        self.assertEqual(self.result["theorem"]["figure_eight_occupied_arm_lower_bound"], 4)

    def test_minimal_n9_figure_eight_certificate(self):
        witness = carrier_descriptor(3, 0, 30, 0)
        self.assertEqual(witness["topology"], "two_carrier_figure_eight")
        directions = []
        for group in witness["groups"]:
            self.assertEqual(group["affine_rank"], 1)
            self.assertEqual(group["distinct_deck_addresses"], 2)
            p, q = (item["deck_address"] for item in group["items"])
            directions.append((q[0] - p[0], q[1] - p[1]))
        self.assertNotEqual(directions[0][0] * directions[1][1] - directions[0][1] * directions[1][0], 0)

    def test_checked_in_result_reproduces(self):
        self.assertEqual(json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8")), self.result)


if __name__ == "__main__":
    unittest.main()
