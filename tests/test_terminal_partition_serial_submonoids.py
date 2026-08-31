import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts" / "terminal_partition_serial_submonoids.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_serial_submonoids_certificate.json"
SPEC = importlib.util.spec_from_file_location("serial_submonoids", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SerialSubmonoidTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_exact_census_and_histogram_mass(self):
        self.assertEqual(self.artifact["submonoid_count"], 228)
        self.assertEqual(sum(self.artifact["submonoid_size_histogram"].values()), 228)

    def test_maximal_proper_submonoids(self):
        self.assertEqual(len(self.artifact["maximal_proper_submonoids"]), 5)
        self.assertEqual(sorted(map(len, self.artifact["maximal_proper_submonoids"])), [9, 9, 12, 12, 14])

    def test_symmetry_orbit_profile(self):
        symmetry = self.artifact["symmetry"]
        self.assertEqual(symmetry["orbit_count"], 95)
        self.assertEqual(symmetry["orbit_size_histogram"], {"1": 22, "2": 43, "4": 30})
        self.assertEqual(symmetry["reversal_stable_count"], 32)
        self.assertEqual(symmetry["lane_swap_stable_count"], 84)

    def test_invalid_symmetry_map_fails_closed(self):
        with self.assertRaises(ValueError):
            MODULE.image(frozenset({6}), list(range(14)))

    def test_tampering_fails_closed(self):
        value = json.loads(json.dumps(self.artifact))
        value["submonoid_count"] = 227
        with self.assertRaises(ValueError):
            MODULE.validate_artifact(value)


if __name__ == "__main__":
    unittest.main()
