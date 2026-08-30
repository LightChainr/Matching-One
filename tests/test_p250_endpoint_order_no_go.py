import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p250_endpoint_order_no_go import (  # noqa: E402
    DEFAULT_OUTPUT,
    build_result,
    endpoint_no_go,
    endpoint_symbol,
    p333_rectangle_gate,
    runner_embedding_audit,
)


class EndpointOrderNoGoTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_xy_and_yx_alias_at_every_radius_and_moment(self):
        for radius in range(1, 9):
            for order in range(1, 13):
                self.assertEqual(endpoint_symbol("xy", radius, order), endpoint_symbol("yx", radius, order))
        self.assertEqual(endpoint_no_go()["violations"], [])

    def test_actual_parent_translation_endpoints_commute(self):
        gate = self.result["actual_parent_translation_gate"]
        self.assertEqual(gate["parent_vertices_checked"], 65)
        self.assertEqual(gate["unit_xy_yx_endpoint_failures"], [])

    def test_p333_control_retains_both_orders(self):
        gate = p333_rectangle_gate()
        self.assertEqual(gate["order_commutator_rank"], 2)
        self.assertEqual(gate["connected_rectangle_rank"], 1)
        self.assertTrue(gate["connected_rectangle_square_zero"])
        self.assertTrue(gate["first_jet_gram_self_adjoint"])

    def test_current_runner_fails_semantic_embedding_gate(self):
        audit = runner_embedding_audit()
        self.assertEqual(audit["embedding_decision"], "NOT_IMPLEMENTED_NO_PHYSICAL_SEMANTICS")
        self.assertEqual(audit["current_ProjectiveLegIndex_public_methods"], ["chi4", "scalar"])
        self.assertNotIn("morphism", audit["charged_rows_parameters"])
        self.assertTrue(any("intermediate" in item for item in audit["minimum_new_intermediate_fields"]))

    def test_checked_in_result_reproduces(self):
        self.assertEqual(json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8")), self.result)


if __name__ == "__main__":
    unittest.main()
