import copy
import json
import unittest

from scripts.gadget_graph_canonical import (
    DEFAULT_MANIFEST,
    canonical_graph,
    decode_graph,
    encode_graph,
    graph_orbit_catalog,
    relabel_graph,
    validate_graph,
    validate_manifest,
)
from scripts.terminal_partition_canonical import full_symmetric_group


class GadgetGraphCanonicalTests(unittest.TestCase):
    def test_encode_decode_round_trip(self):
        edges = ((0, 3), (1, 3), (2, 3))
        encoded = encode_graph(4, 3, edges)
        self.assertEqual(decode_graph(encoded), (3, 4, edges))

    def test_terminal_and_internal_relabeling_invariance(self):
        group = full_symmetric_group(3)
        graph = ((0, 3), (1, 3), (1, 4), (2, 4), (3, 4))
        key = canonical_graph(5, 3, graph, group)
        for mapping in ((1, 2, 0, 3, 4), (0, 1, 2, 4, 3), (2, 0, 1, 4, 3)):
            moved = relabel_graph(5, 3, graph, mapping)
            self.assertEqual(canonical_graph(5, 3, moved, group), key)

    def test_explicit_group_does_not_assume_full_symmetry(self):
        identity = ((0, 1, 2),)
        graph = ((0, 3),)
        moved = ((1, 3),)
        self.assertNotEqual(canonical_graph(4, 3, graph, identity), canonical_graph(4, 3, moved, identity))
        full = full_symmetric_group(3)
        self.assertEqual(canonical_graph(4, 3, graph, full), canonical_graph(4, 3, moved, full))

    def test_exhaustive_reference_orbits(self):
        self.assertEqual(len(graph_orbit_catalog(4, 3, full_symmetric_group(3))), 20)
        catalog = graph_orbit_catalog(5, 4, full_symmetric_group(4))
        self.assertEqual(len(catalog), 90)
        self.assertEqual(sum(catalog.values()), 1024)

    def test_checked_in_manifest_and_boundary(self):
        manifest = json.loads(DEFAULT_MANIFEST.read_text(encoding="utf-8"))
        summary = validate_manifest(manifest)
        self.assertEqual(summary["audited"]["3"]["canonical_orbits"], 20)
        self.assertEqual(summary["audited"]["4"]["canonical_orbits"], 90)
        tampered = copy.deepcopy(manifest)
        tampered["exhaustive_checks"][1]["canonical_orbits"] += 1
        with self.assertRaisesRegex(ValueError, "orbit count drift"):
            validate_manifest(tampered)


if __name__ == "__main__":
    unittest.main()
