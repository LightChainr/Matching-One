import importlib.util
import json
from fractions import Fraction
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "terminal_partition_d4_orbit_compression.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_d4_orbit_compression_certificate.json"
SPEC = importlib.util.spec_from_file_location("terminal_partition_d4_orbit_compression", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TerminalPartitionD4OrbitCompressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact_regenerates_exactly(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_d4_compression_has_seven_orbits(self):
        self.assertEqual(len(self.artifact["d4_orbits"]), 7)
        self.assertEqual(sum(map(len, self.artifact["d4_orbits"])), 15)

    def test_deterministic_quotient_has_26_ambiguous_pairs(self):
        boundary = self.artifact["deterministic_quotient"]
        self.assertEqual(boundary["orbit_pairs"], 49)
        self.assertEqual(boundary["ambiguous_pairs"], 26)
        self.assertEqual(boundary["unambiguous_pairs"], 23)

    def test_smallest_counterexample_has_two_outputs(self):
        witness = self.artifact["deterministic_quotient"]["smallest_counterexample"]
        self.assertEqual(witness["labelled_pair_count"], 2)
        self.assertEqual(len(witness["output_counts"]), 2)
        self.assertEqual(len(witness["labelled_witnesses"]), 2)

    def test_every_uniform_kernel_cell_is_normalized(self):
        for row in self.artifact["uniform_orbit_pair_kernel"]:
            for cell in row:
                self.assertEqual(sum(Fraction(value) for value in cell), 1)

    def test_uniform_orbit_average_is_not_associative(self):
        self.assertEqual(self.artifact["averaging_boundary"]["basis_triples_checked"], 343)
        self.assertEqual(self.artifact["averaging_boundary"]["associativity_failures"], 166)
        failure = self.artifact["averaging_boundary"]["first_failure"]
        self.assertNotEqual(failure["left_grouped"], failure["right_grouped"])


if __name__ == "__main__":
    unittest.main()
