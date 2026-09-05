import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts" / "terminal_partition_serial_shortlex.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_serial_shortlex_certificate.json"
SPEC = importlib.util.spec_from_file_location("serial_shortlex", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SerialShortlexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = MODULE.multiplication_table()
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_every_word_and_trace(self):
        for profile in self.artifact["profiles"]:
            for record in profile["normal_forms"]:
                value, trace = MODULE.evaluate_word(record["word"], self.table)
                self.assertEqual(value, record["target"])
                self.assertEqual(trace, record["evaluation_trace"])

    def test_all_minimum_generating_sets_are_covered(self):
        self.assertEqual(self.artifact["profile_count"], 8)
        self.assertEqual(self.artifact["diameter_histogram"], {"4": 4, "5": 4})


if __name__ == "__main__":
    unittest.main()
