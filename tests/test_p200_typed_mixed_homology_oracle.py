from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p200_typed_mixed_homology_oracle import render  # noqa: E402


class P200TypedMixedHomologyOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = render()

    def test_real_tiny_gaussian_hnf_lineage(self) -> None:
        geometry = self.payload["geometry"]
        self.assertEqual(geometry["final_column_HNF"], [[10, 3], [0, 1]])
        self.assertEqual(geometry["N2_column_HNF"], [[2, 1], [0, 1]])
        self.assertEqual(geometry["N5_column_HNF"], [[5, 3], [0, 1]])

    def test_exact_symmetry_and_nonzero_witnesses(self) -> None:
        checks = self.payload["exhaustive_checks"]
        self.assertEqual(checks["configurations"], 1024)
        self.assertTrue(checks["join_path_symmetry_every_configuration_and_colour"])
        self.assertTrue(checks["typed_matching_layer_swap_negates_odd_rows"])
        self.assertTrue(self.payload["witnesses"]["nonzero_ambient_colour_odd"])
        self.assertTrue(self.payload["witnesses"]["nonzero_partition_residual_colour_odd"])

    def test_local_incidence_normalization(self) -> None:
        normalization = self.payload["exact_local_normalization"]
        self.assertEqual(normalization["E_J_local_black_at_p_half"]["exact"], "499/1024")
        self.assertEqual(normalization["E_JB_minus_JW"]["exact"], "0")
        self.assertEqual(normalization["Var_JB_minus_JW"]["exact"], "681/512")

    def test_frozen_artifact_matches(self) -> None:
        artifact = json.loads(
            (ROOT / "results" / "exact-cover-character-oracles" / "p200_typed_mixed_homology.json").read_text(encoding="utf-8")
        )
        self.assertEqual(artifact, self.payload)


if __name__ == "__main__":
    unittest.main()
