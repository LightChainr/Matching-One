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


if __name__ == "__main__":
    unittest.main()
