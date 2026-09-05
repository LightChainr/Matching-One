import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts" / "terminal_partition_serial_automorphisms.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_serial_automorphisms_certificate.json"
SPEC = importlib.util.spec_from_file_location("serial_automorphisms", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SerialAutomorphismTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = MODULE.multiplication_table()
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_complete_candidate_count(self):
        self.assertEqual(self.artifact["ordered_generator_image_candidates"], 48)

    def test_automorphism_group_is_c2(self):
        maps = [tuple(value) for value in self.artifact["automorphisms"]]
        self.assertEqual(len(maps), 2)
        self.assertTrue(all(MODULE.compose_maps(value, value) == tuple(range(15)) for value in maps))

    def test_all_anti_maps_reverse_products(self):
        for value in self.artifact["anti_automorphisms"]:
            for a in range(15):
                for b in range(15):
                    self.assertEqual(value[self.table[a][b]], self.table[value[b]][value[a]])


if __name__ == "__main__":
    unittest.main()
