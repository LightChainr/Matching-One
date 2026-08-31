import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "terminal_partition_d4_equivariance.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_d4_equivariance_certificate.json"
SPEC = importlib.util.spec_from_file_location("terminal_partition_d4_equivariance", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TerminalPartitionD4EquivarianceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.partitions = MODULE.enumerate_rgs(4)
        cls.group = MODULE.d4_group()
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact_regenerates_exactly(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_d4_is_valid_order_eight_group(self):
        self.assertEqual(len(self.group), 8)
        self.assertEqual(MODULE.validate_group(self.group, 4), self.group)

    def test_d4_orbits_partition_all_states(self):
        orbits = MODULE.d4_orbits(self.partitions, self.group)
        self.assertEqual(sum(len(orbit) for orbit in orbits), 15)
        self.assertEqual(len({item for orbit in orbits for item in orbit}), 15)

    def test_join_is_d4_equivariant_for_all_1800_cases(self):
        for permutation in self.group:
            for left in self.partitions:
                for right in self.partitions:
                    self.assertEqual(
                        MODULE.partition_join(
                            MODULE.apply_permutation(left, permutation),
                            MODULE.apply_permutation(right, permutation),
                        ),
                        MODULE.apply_permutation(MODULE.partition_join(left, right), permutation),
                    )

    def test_declared_interface_is_coordinate_covariant(self):
        standard = {
            "left_outer": (0, 1),
            "left_interface": (2, 3),
            "right_interface": (0, 1),
            "right_outer": (2, 3),
        }
        for permutation in self.group:
            moved = {key: tuple(permutation[item] for item in value) for key, value in standard.items()}
            for left in self.partitions:
                for right in self.partitions:
                    self.assertEqual(
                        MODULE.declared_interface_glue(
                            MODULE.apply_permutation(left, permutation),
                            MODULE.apply_permutation(right, permutation),
                            **moved,
                        ),
                        MODULE.declared_interface_glue(left, right, **standard),
                    )

    def test_known_orbit_stabilizers(self):
        for row in self.artifact["d4_orbits"]:
            self.assertEqual(row["orbit_size"] * row["stabilizer_size"], 8)

    def test_invalid_declarations_fail_closed(self):
        discrete = (0, 1, 2, 3)
        with self.assertRaises(ValueError):
            MODULE.declared_interface_glue(discrete, discrete, left_outer=(0, 1), left_interface=(1, 2), right_interface=(0, 1), right_outer=(2, 3))
        with self.assertRaises(ValueError):
            MODULE.declared_interface_glue(discrete, discrete, left_outer=(0,), left_interface=(1, 2), right_interface=(0, 1), right_outer=(2, 3))

    def test_tampered_artifact_fails_closed(self):
        payload = json.loads(json.dumps(self.artifact))
        payload["counts"]["group_order"] = 7
        with self.assertRaises(ValueError):
            MODULE.validate_artifact(payload)


if __name__ == "__main__":
    unittest.main()
