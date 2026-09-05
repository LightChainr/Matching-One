import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "terminal_partition_incidence_algebra.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_incidence_algebra_certificate.json"
SPEC = importlib.util.spec_from_file_location("terminal_partition_incidence_algebra", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TerminalPartitionIncidenceAlgebraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.partitions = MODULE.enumerate_rgs(4)
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact_regenerates_exactly(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_partial_order_is_reflexive_antisymmetric_transitive(self):
        for left in self.partitions:
            self.assertTrue(MODULE.refines(left, left))
            for right in self.partitions:
                if MODULE.refines(left, right) and MODULE.refines(right, left):
                    self.assertEqual(left, right)
                for third in self.partitions:
                    if MODULE.refines(left, right) and MODULE.refines(right, third):
                        self.assertTrue(MODULE.refines(left, third))

    def test_closed_form_matches_independent_recurrence(self):
        zeta, mobius = MODULE.matrices(self.partitions)
        self.assertEqual(mobius, MODULE.recurrence_mobius(zeta))

    def test_zeta_and_mobius_are_two_sided_inverses(self):
        zeta, mobius = MODULE.matrices(self.partitions)
        identity = [[int(i == j) for j in range(15)] for i in range(15)]
        self.assertEqual(MODULE.matmul(zeta, mobius), identity)
        self.assertEqual(MODULE.matmul(mobius, zeta), identity)

    def test_all_integer_vectors_round_trip(self):
        zeta, mobius = MODULE.matrices(self.partitions)
        for shift in range(-3, 4):
            vector = tuple((index + 1) * shift - 2 for index in range(15))
            self.assertEqual(MODULE.transform(MODULE.transform(vector, zeta), mobius), vector)

    def test_extreme_interval_values(self):
        discrete = (0, 1, 2, 3)
        connected = (0, 0, 0, 0)
        self.assertEqual(MODULE.mobius_closed_form(discrete, connected), -6)
        self.assertEqual(MODULE.mobius_closed_form(connected, discrete), 0)

    def test_exact_determinants_are_units(self):
        zeta, mobius = MODULE.matrices(self.partitions)
        self.assertEqual(abs(MODULE.determinant_bareiss(zeta)), 1)
        self.assertEqual(abs(MODULE.determinant_bareiss(mobius)), 1)


if __name__ == "__main__":
    unittest.main()
