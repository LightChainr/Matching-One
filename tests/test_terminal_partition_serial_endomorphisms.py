import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts" / "terminal_partition_serial_endomorphisms.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_serial_endomorphisms_certificate.json"
SPEC = importlib.util.spec_from_file_location("serial_endomorphisms", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SerialEndomorphismTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = MODULE.multiplication_table()
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_every_map_preserves_products(self):
        for record in self.artifact["endomorphisms"]:
            mapping = record["mapping"]
            for a in range(15):
                for b in range(15):
                    self.assertEqual(mapping[self.table[a][b]], self.table[mapping[a]][mapping[b]])

    def test_idempotent_indices_are_retractions(self):
        for index in self.artifact["idempotent_retraction_indices"]:
            mapping = self.artifact["endomorphisms"][index]["mapping"]
            self.assertEqual(MODULE.compose_maps(mapping, mapping), tuple(mapping))

    def test_orbits_partition_census(self):
        self.assertEqual(sorted(i for orbit in self.artifact["automorphism_conjugacy_orbits"] for i in orbit), list(range(self.artifact["endomorphism_count"])))

    def test_map_width_mismatch_fails_closed(self):
        with self.assertRaises(ValueError):
            MODULE.compose_maps([0], [0, 1])

    def test_tampering_fails_closed(self):
        value = json.loads(json.dumps(self.artifact))
        value["endomorphism_count"] += 1
        with self.assertRaises(ValueError):
            MODULE.validate_artifact(value)


if __name__ == "__main__":
    unittest.main()
