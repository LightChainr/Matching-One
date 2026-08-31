import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "terminal_partition_serial_reversal.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_serial_reversal_certificate.json"
SPEC = importlib.util.spec_from_file_location("terminal_partition_serial_reversal", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TerminalPartitionSerialReversalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.partitions = MODULE.enumerate_rgs(4)
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact_regenerates_exactly(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_reversal_is_an_involution_on_all_states(self):
        for state in self.partitions:
            self.assertEqual(MODULE.reverse_ports(MODULE.reverse_ports(state)), state)

    def test_all_225_products_reverse_order(self):
        for left in self.partitions:
            for right in self.partitions:
                self.assertEqual(
                    MODULE.reverse_ports(MODULE.serial_compose(left, right)),
                    MODULE.serial_compose(MODULE.reverse_ports(right), MODULE.reverse_ports(left)),
                )

    def test_wire_identity_is_fixed(self):
        self.assertEqual(MODULE.reverse_ports((0, 1, 0, 1)), (0, 1, 0, 1))

    def test_vector_transport_is_involutive_for_any_catalog_order(self):
        catalog = tuple(reversed(self.partitions))
        values = tuple(range(15))
        self.assertEqual(
            MODULE.reverse_index_vector(MODULE.reverse_index_vector(values, catalog), catalog),
            values,
        )

    def test_invalid_state_width_fails_closed(self):
        with self.assertRaises(ValueError):
            MODULE.reverse_ports((0, 1, 2))

    def test_incomplete_or_duplicate_catalog_fails_closed(self):
        values = tuple(range(15))
        with self.assertRaises(ValueError):
            MODULE.reverse_index_vector(values, self.partitions[:-1])
        duplicate = self.partitions[:-1] + (self.partitions[0],)
        with self.assertRaises(ValueError):
            MODULE.reverse_index_vector(values, duplicate)

    def test_mismatched_vector_length_fails_closed(self):
        with self.assertRaises(ValueError):
            MODULE.reverse_index_vector(tuple(range(14)), self.partitions)

    def test_tampered_artifact_fails_closed(self):
        payload = json.loads(json.dumps(self.artifact))
        payload["reversal_index_map"][0] = 1
        with self.assertRaises(ValueError):
            MODULE.validate_artifact(payload)


if __name__ == "__main__":
    unittest.main()
