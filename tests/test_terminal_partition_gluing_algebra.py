import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "terminal_partition_gluing_algebra.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_gluing_algebra_certificate.json"
SPEC = importlib.util.spec_from_file_location("terminal_partition_gluing_algebra", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TerminalPartitionGluingAlgebraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.partitions = MODULE.enumerate_rgs(4)
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact_regenerates_exactly(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_four_terminal_catalog_has_bell_number_15(self):
        self.assertEqual(len(self.partitions), 15)

    def test_join_laws_exhaustively(self):
        discrete = (0, 1, 2, 3)
        connected = (0, 0, 0, 0)
        for left in self.partitions:
            self.assertEqual(MODULE.partition_join(left, discrete), left)
            self.assertEqual(MODULE.partition_join(left, connected), connected)
            self.assertEqual(MODULE.partition_join(left, left), left)
            for right in self.partitions:
                self.assertEqual(MODULE.partition_join(left, right), MODULE.partition_join(right, left))

    def test_join_associativity_for_all_3375_triples(self):
        for left in self.partitions:
            for middle in self.partitions:
                for right in self.partitions:
                    self.assertEqual(
                        MODULE.partition_join(MODULE.partition_join(left, middle), right),
                        MODULE.partition_join(left, MODULE.partition_join(middle, right)),
                    )

    def test_interface_table_is_deterministic_15_by_15(self):
        table = self.artifact["two_port_interface"]["deterministic_output_index_table"]
        self.assertEqual(len(table), 15)
        self.assertTrue(all(len(row) == 15 for row in table))
        self.assertTrue(all(0 <= output < 15 for row in table for output in row))

    def test_all_225_interface_outputs_match_tiny_graphs(self):
        for left in self.partitions:
            for right in self.partitions:
                self.assertEqual(MODULE.interface_glue(left, right), MODULE.graph_composition_output(left, right))

    def test_known_interface_extremes(self):
        discrete = (0, 1, 2, 3)
        connected = (0, 0, 0, 0)
        self.assertEqual(MODULE.interface_glue(discrete, discrete), discrete)
        self.assertEqual(MODULE.interface_glue(connected, connected), connected)

    def test_bilinear_composition_preserves_exact_mass(self):
        table = MODULE.interface_table(self.partitions)
        denominator = sum(range(1, 16))
        left = tuple(Fraction(index + 1, denominator) for index in range(15))
        right = tuple(reversed(left))
        output = MODULE.bilinear_compose(left, right, table)
        self.assertEqual(sum(output), Fraction(1))
        self.assertTrue(all(value >= 0 for value in output))

    def test_invalid_inputs_fail_closed(self):
        with self.assertRaises(ValueError):
            MODULE.partition_join((0, 1), (0, 1, 2))
        with self.assertRaises(ValueError):
            MODULE.interface_glue((0, 1, 2), (0, 1, 2, 3))
        with self.assertRaises(TypeError):
            MODULE.bilinear_compose([0.0] * 15, [0] * 15, MODULE.interface_table(self.partitions))
        with self.assertRaises(ValueError):
            MODULE.bilinear_compose([0] * 14, [0] * 15, MODULE.interface_table(self.partitions))

    def test_tampered_artifact_fails_closed(self):
        payload = json.loads(json.dumps(self.artifact))
        payload["join_cayley_table"][0][0] = 14
        with self.assertRaises(ValueError):
            MODULE.validate_artifact(payload)


if __name__ == "__main__":
    unittest.main()
