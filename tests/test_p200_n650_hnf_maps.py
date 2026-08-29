from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p200_n650_hnf_maps import render  # noqa: E402


class P200N650HNFMapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = render()

    def test_real_hnf_quotient_maps_commute(self) -> None:
        finals = []
        for lineage in self.payload["lineages"]:
            self.assertTrue(lineage["composition_to_N65_equal_by_both_paths"])
            self.assertTrue(lineage["N130_N325_label_pair_is_bijective_on_each_N65_fiber"])
            self.assertEqual(lineage["map_cardinalities"]["N650_to_N130"]["fiber_size"], 5)
            self.assertEqual(lineage["map_cardinalities"]["N650_to_N325"]["fiber_size"], 2)
            finals.append(lineage["lattices"]["final_N650"]["gaussian"])
        self.assertEqual(finals, [[23, 11], [17, 19]])

    def test_hnf_is_runner_convention(self) -> None:
        hnfs = [lineage["lattices"]["final_N650"]["column_HNF"] for lineage in self.payload["lineages"]]
        self.assertEqual(hnfs, [[[650, 593], [0, 1]], [[650, 343], [0, 1]]])

    def test_frozen_artifact_matches(self) -> None:
        artifact = json.loads(
            (ROOT / "results" / "exact-cover-character-oracles" / "n650_hnf_maps.json").read_text(encoding="utf-8")
        )
        self.assertEqual(artifact, self.payload)


if __name__ == "__main__":
    unittest.main()
