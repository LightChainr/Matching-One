import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "terminal_partition_serial_category.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_serial_category_certificate.json"
SPEC = importlib.util.spec_from_file_location("terminal_partition_serial_category", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TerminalPartitionSerialCategoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.partitions = MODULE.enumerate_rgs(4)
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact_regenerates_exactly(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_wire_partition_is_two_sided_identity(self):
        identity = (0, 1, 0, 1)
        for state in self.partitions:
            self.assertEqual(MODULE.serial_compose(identity, state), state)
            self.assertEqual(MODULE.serial_compose(state, identity), state)

    def test_all_3375_triples_are_associative(self):
        for left in self.partitions:
            for middle in self.partitions:
                for right in self.partitions:
                    self.assertEqual(
                        MODULE.serial_compose(MODULE.serial_compose(left, middle), right),
                        MODULE.serial_compose(left, MODULE.serial_compose(middle, right)),
                    )

    def test_all_3375_triples_match_explicit_graph_oracle(self):
        for left in self.partitions:
            for middle in self.partitions:
                for right in self.partitions:
                    self.assertEqual(
                        MODULE.serial_compose(MODULE.serial_compose(left, middle), right),
                        MODULE.triple_graph_output(left, middle, right),
                    )

    def test_serial_product_is_order_sensitive(self):
        failures = self.artifact["monoid"]["ordered_commutativity_failures"]
        self.assertGreater(failures, 0)

    def test_exact_measures_are_associative(self):
        table = MODULE.interface_table(self.partitions)
        denominator = sum(range(1, 16))
        a = tuple(Fraction(index + 1, denominator) for index in range(15))
        b = tuple(reversed(a))
        c = tuple(a[(index + 4) % 15] for index in range(15))
        self.assertEqual(
            MODULE.compose_measures(MODULE.compose_measures(a, b, table), c, table),
            MODULE.compose_measures(a, MODULE.compose_measures(b, c, table), table),
        )

    def test_invalid_width_fails_closed(self):
        with self.assertRaises(ValueError):
            MODULE.serial_compose((0, 1, 2), (0, 1, 2, 3))
        with self.assertRaises(ValueError):
            MODULE.triple_graph_output((0, 1, 2, 3), (0, 1), (0, 1, 2, 3))

    def test_tampered_artifact_fails_closed(self):
        payload = json.loads(json.dumps(self.artifact))
        payload["monoid"]["identity_index"] = 0
        with self.assertRaises(ValueError):
            MODULE.validate_artifact(payload)


if __name__ == "__main__":
    unittest.main()
