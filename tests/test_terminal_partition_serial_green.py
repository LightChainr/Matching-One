import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts" / "terminal_partition_serial_green.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_serial_green_certificate.json"
SPEC = importlib.util.spec_from_file_location("serial_green", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SerialGreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_d_equals_j(self):
        classes = self.artifact["green_classes"]
        self.assertEqual(classes["D"], classes["J"])
        self.assertEqual(sorted(map(len, classes["J"])), [2, 4, 9])

    def test_all_elements_are_regular(self):
        self.assertEqual(self.artifact["regular_element_indices"], list(range(15)))

    def test_ideal_lattice_is_exact_chain(self):
        ideals = self.artifact["two_sided_ideals"]
        self.assertEqual(list(map(len, ideals)), [0, 4, 13, 15])
        self.assertTrue(all(set(ideals[i]) < set(ideals[i + 1]) for i in range(3)))

    def test_unique_nontrivial_h_class(self):
        self.assertEqual([c for c in self.artifact["green_classes"]["H"] if len(c) > 1], [[6, 8]])

    def test_tampering_fails_closed(self):
        value = json.loads(json.dumps(self.artifact))
        value["two_sided_ideals"][1].append(0)
        with self.assertRaises(ValueError):
            MODULE.validate_artifact(value)


if __name__ == "__main__":
    unittest.main()
