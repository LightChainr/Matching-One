import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts" / "terminal_partition_serial_congruences.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_serial_congruences_certificate.json"
SPEC = importlib.util.spec_from_file_location("serial_congruences", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SerialCongruenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = MODULE.multiplication_table()
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_complete_census_is_compatible(self):
        for record in self.artifact["congruences"]:
            self.assertIsNone(MODULE.congruence_failure(record["labels"], self.table))

    def test_binary_join_closure(self):
        values = [tuple(record["labels"]) for record in self.artifact["congruences"]]
        for left in values:
            for right in values:
                self.assertIn(MODULE.join(left, right, self.table), values)

    def test_raw_merge_witness_is_real(self):
        witness = self.artifact["raw_merge_failure_witness"]
        labels = list(range(15))
        labels[witness["unsupported_merge"][1]] = labels[witness["unsupported_merge"][0]]
        self.assertEqual(MODULE.congruence_failure(labels, self.table), {key: value for key, value in witness.items() if key != "unsupported_merge"})

    def test_invalid_width_fails_closed(self):
        with self.assertRaises(ValueError):
            MODULE.congruence_failure([0], self.table)

    def test_tampering_fails_closed(self):
        value = json.loads(json.dumps(self.artifact))
        value["congruence_count"] += 1
        with self.assertRaises(ValueError):
            MODULE.validate_artifact(value)


if __name__ == "__main__":
    unittest.main()
