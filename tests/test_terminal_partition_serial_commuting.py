import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts" / "terminal_partition_serial_commuting.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_serial_commuting_certificate.json"
SPEC = importlib.util.spec_from_file_location("serial_commuting", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SerialCommutingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_center_is_unique(self):
        self.assertEqual(self.artifact["center_indices"], [6])

    def test_commutative_submonoid_census(self):
        self.assertEqual(self.artifact["commutative_submonoid_count"], 29)
        self.assertEqual(len(self.artifact["maximal_commutative_submonoids"]), 11)

    def test_every_reported_maximal_sector_commutes(self):
        table = MODULE.multiplication_table()
        for sector in self.artifact["maximal_commutative_submonoids"]:
            self.assertTrue(all(table[a][b] == table[b][a] for a in sector for b in sector))

    def test_center_removed_graph_invariants(self):
        graph = self.artifact["commuting_graph_without_center"]
        self.assertEqual(len(graph["components"]), 1)
        self.assertEqual(graph["clique_number"], 3)


if __name__ == "__main__":
    unittest.main()
