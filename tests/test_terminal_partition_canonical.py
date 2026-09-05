#!/usr/bin/env python3
"""Tests for the Issue 13 terminal-partition canonicalization primitive."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import terminal_partition_canonical as canonical  # noqa: E402


MANIFEST_PATH = ROOT / "analysis" / "terminal_partition_canonical_manifest.json"


class TerminalPartitionCanonicalTests(unittest.TestCase):
    def setUp(self) -> None:
        self.manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))

    def test_rgs_enumeration_and_block_round_trip(self) -> None:
        for n, expected in ((3, 5), (4, 15)):
            with self.subTest(n=n):
                values = canonical.enumerate_rgs(n)
                self.assertEqual(len(values), expected)
                self.assertEqual(len(values), len(set(values)))
                for value in values:
                    self.assertEqual(canonical.blocks_to_rgs(canonical.rgs_to_blocks(value), n), value)
        self.assertEqual(canonical.blocks_to_rgs(((2, 0), (3, 1)), 4), (0, 1, 0, 1))

    def test_full_symmetric_orbits_match_integer_partition_types(self) -> None:
        expected = {
            3: ((0, 0, 0), (0, 0, 1), (0, 1, 2)),
            4: ((0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 1, 1), (0, 0, 1, 2), (0, 1, 2, 3)),
        }
        for n in (3, 4):
            with self.subTest(n=n):
                catalog = canonical.orbit_catalog(n, canonical.full_symmetric_group(n))
                self.assertEqual(tuple(catalog), expected[n])
                self.assertEqual(sum(map(len, catalog.values())), len(canonical.enumerate_rgs(n)))

    def test_explicit_group_action_is_invariant_and_not_assumed_symmetric(self) -> None:
        cyclic_four = (
            (0, 1, 2, 3),
            (1, 2, 3, 0),
            (2, 3, 0, 1),
            (3, 0, 1, 2),
        )
        partition = (0, 0, 1, 2)
        representative = canonical.canonical_orbit(partition, cyclic_four)
        for permutation in cyclic_four:
            moved = canonical.apply_permutation(partition, permutation)
            self.assertEqual(canonical.canonical_orbit(moved, cyclic_four), representative)
        self.assertGreater(
            len(canonical.orbit_catalog(4, cyclic_four)),
            len(canonical.orbit_catalog(4, canonical.full_symmetric_group(4))),
        )

    def test_invalid_encodings_partitions_and_groups_fail_closed(self) -> None:
        for invalid in ((1, 0), (0, 2), (0, -1), (0, True)):
            with self.subTest(rgs=invalid), self.assertRaises(ValueError):
                canonical.validate_rgs(invalid)
        with self.assertRaisesRegex(ValueError, "multiple blocks"):
            canonical.blocks_to_rgs(((0, 1), (1, 2)), 3)
        with self.assertRaisesRegex(ValueError, "every terminal"):
            canonical.blocks_to_rgs(((0,), (2,)), 3)
        with self.assertRaisesRegex(ValueError, "composition"):
            canonical.validate_group(((0, 1, 2), (1, 0, 2), (0, 2, 1)), 3)
        with self.assertRaisesRegex(ValueError, "bijection"):
            canonical.validate_group(((0, 1, 2), (0, 0, 2)), 3)


if __name__ == "__main__":
    unittest.main()
