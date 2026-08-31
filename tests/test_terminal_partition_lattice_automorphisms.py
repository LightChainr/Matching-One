import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "terminal_partition_lattice_automorphisms.py"
ARTIFACT = ROOT / "analysis" / "terminal_partition_lattice_automorphisms_certificate.json"
SPEC = importlib.util.spec_from_file_location("terminal_partition_lattice_automorphisms", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class TerminalPartitionLatticeAutomorphismsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.partitions = MODULE.enumerate_rgs(4)
        cls.atoms = MODULE.atom_catalog(cls.partitions)
        cls.maps = MODULE.enumerate_lattice_automorphisms(cls.partitions)
        cls.artifact = MODULE.build_artifact()

    def test_committed_artifact_regenerates_exactly(self):
        self.assertEqual(json.loads(ARTIFACT.read_text()), self.artifact)

    def test_six_atoms_and_720_candidates(self):
        self.assertEqual(len(self.atoms), 6)
        self.assertEqual(self.artifact["enumeration"]["atom_permutations_checked"], 720)

    def test_exactly_24_lattice_automorphisms(self):
        self.assertEqual(len(self.maps), 24)
        self.assertEqual(len(set(self.maps)), 24)

    def test_every_automorphism_preserves_join(self):
        index = {partition: position for position, partition in enumerate(self.partitions)}
        for mapping in self.maps:
            for left_index, left in enumerate(self.partitions):
                for right_index, right in enumerate(self.partitions):
                    self.assertEqual(
                        mapping[index[MODULE.partition_join(left, right)]],
                        index[MODULE.partition_join(self.partitions[mapping[left_index]], self.partitions[mapping[right_index]])],
                    )

    def test_automorphisms_equal_terminal_relabelings(self):
        self.assertEqual(self.maps, MODULE.terminal_relabeling_maps(self.partitions))

    def test_partition_type_orbits(self):
        rows = self.artifact["partition_type_orbits"]
        self.assertEqual([row["orbit_size"] for row in rows], [1, 4, 3, 6, 1])
        self.assertEqual([row["stabilizer_size"] for row in rows], [24, 6, 8, 4, 24])

    def test_non_graph_atom_permutation_is_rejected(self):
        accepted = {
            tuple(mapping) for mapping in self.artifact["automorphism_maps"]
        }
        rejected = 0
        from itertools import permutations
        for candidate in permutations(range(6)):
            mapping = MODULE.induced_map_from_atom_permutation(candidate, self.partitions, self.atoms)
            if mapping is None:
                rejected += 1
            else:
                self.assertIn(mapping, accepted)
        self.assertEqual(rejected, 696)

    def test_invalid_atom_permutation_fails_closed(self):
        with self.assertRaises(ValueError):
            MODULE.induced_map_from_atom_permutation((0, 0, 1, 2, 3, 4), self.partitions, self.atoms)
        with self.assertRaises(ValueError):
            MODULE.join_atoms((6,), self.atoms)

    def test_tampered_artifact_fails_closed(self):
        payload = json.loads(json.dumps(self.artifact))
        payload["enumeration"]["lattice_automorphisms"] = 25
        with self.assertRaises(ValueError):
            MODULE.validate_artifact(payload)


if __name__ == "__main__":
    unittest.main()
