import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts" / "terminal_partition_serial_generators.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_serial_generators_certificate.json"
SPEC = importlib.util.spec_from_file_location("serial_generators", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SerialGeneratorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = MODULE.multiplication_table()
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_rank_and_all_minimum_sets(self):
        self.assertEqual(self.artifact["monoid_rank"], 3)
        self.assertEqual(self.artifact["minimal_generating_set_count"], 8)

    def test_reported_sets_generate_all_states(self):
        for seed in self.artifact["minimal_generating_sets"]:
            self.assertEqual(MODULE.generated_closure(seed, self.table), frozenset(range(15)))

    def test_word_metrics_are_exact(self):
        diameters = []
        for profile in self.artifact["word_metric_profiles"]:
            lengths = MODULE.shortest_word_lengths(profile["generators"], self.table)
            self.assertEqual(list(lengths), profile["shortest_word_lengths"])
            diameters.append(max(lengths))
        self.assertEqual(sorted(diameters), [4, 4, 4, 4, 5, 5, 5, 5])

    def test_invalid_generator_fails_closed(self):
        with self.assertRaises(ValueError):
            MODULE.generated_closure([15], self.table)

    def test_nongenerating_word_metric_fails_closed(self):
        with self.assertRaises(ValueError):
            MODULE.shortest_word_lengths([6], self.table)

    def test_tampering_fails_closed(self):
        value = json.loads(json.dumps(self.artifact))
        value["monoid_rank"] = 2
        with self.assertRaises(ValueError):
            MODULE.validate_artifact(value)


if __name__ == "__main__":
    unittest.main()
