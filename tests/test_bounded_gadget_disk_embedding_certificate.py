import importlib.util
import json
from itertools import combinations
from pathlib import Path
import unittest

from scripts.bounded_gadget_planarity_certificate import minimum_orientable_genus


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "bounded_gadget_disk_embedding_certificate.py"
ARTIFACT = ROOT / "analysis" / "bounded_gadget_disk_embedding_certificate.json"
SPEC = importlib.util.spec_from_file_location("bounded_gadget_disk_embedding_certificate", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class BoundedGadgetDiskEmbeddingCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact_regenerates_exactly(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_four_cycle_accepts_only_its_boundary_order(self):
        cycle = ((0, 1), (1, 2), (2, 3), (0, 3))
        self.assertTrue(MODULE.disk_embedding(4, cycle, (0, 1, 2, 3))["disk_planar"])
        self.assertFalse(MODULE.disk_embedding(4, cycle, (0, 1, 3, 2))["disk_planar"])

    def test_four_spoke_star_accepts_all_order_classes(self):
        star = ((0, 4), (1, 4), (2, 4), (3, 4))
        self.assertTrue(all(MODULE.disk_embedding(5, star, order)["disk_planar"] for order in MODULE.canonical_cyclic_orders(4)))

    def test_k4_is_planar_but_not_disk_planar_with_all_vertices_terminal(self):
        k4 = tuple(combinations(range(4), 2))
        self.assertFalse(MODULE.disk_embedding(4, k4, (0, 1, 2, 3))["disk_planar"])

    def test_k5_is_not_disk_planar(self):
        k5 = tuple(combinations(range(5), 2))
        self.assertFalse(MODULE.disk_embedding(5, k5, (0, 1, 2, 3))["disk_planar"])

    def test_k5_minus_edge_is_planar_but_not_terminal_disk_planar(self):
        graph = tuple(edge for edge in combinations(range(5), 2) if edge != (0, 1))
        self.assertEqual(minimum_orientable_genus(5, graph)[0], 0)
        results = [MODULE.disk_embedding(5, graph, order)["disk_planar"] for order in MODULE.canonical_cyclic_orders(4)]
        self.assertEqual(results, [False, False, False])

    def test_crossing_partition_is_absent_for_disk_witnesses(self):
        for row in self.artifact["rows"]:
            for record in row["records"]:
                for result in record["orders"]:
                    if result["disk_planar"] and row["terminal_count"] == 4:
                        self.assertEqual(result["crossing_subgraphs"], 0)

    def test_crossing_partition_control_is_detected(self):
        graph = ((0, 2), (1, 3))
        count, checked = MODULE.crossing_subgraph_count(4, graph, (0, 1, 2, 3))
        self.assertEqual((count, checked), (1, 4))

    def test_invalid_orders_and_graphs_fail_closed(self):
        cycle = ((0, 1), (1, 2), (2, 0))
        with self.assertRaises(ValueError):
            MODULE.disk_embedding(3, cycle, (0, 1, 1))
        with self.assertRaises(ValueError):
            MODULE.disk_embedding(3, cycle, (0, 1))
        with self.assertRaises(ValueError):
            MODULE.disk_embedding(4, cycle, (0, 1, 2))
        with self.assertRaises(TypeError):
            MODULE.disk_embedding(3, cycle, (0, 1, True))

    def test_tampered_artifact_fails_closed(self):
        payload = json.loads(json.dumps(self.artifact))
        payload["rows"][1]["disk_planar_orbit_order_pairs"] += 1
        with self.assertRaises(ValueError):
            MODULE.validate_artifact(payload)


if __name__ == "__main__":
    unittest.main()
