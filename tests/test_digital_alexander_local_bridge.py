from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import digital_alexander_local_bridge as bridge  # noqa: E402


class DigitalAlexanderLocalBridgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        config = json.loads(
            (ROOT / "analysis" / "digital_alexander_local_bridge_manifest.json").read_text()
        )
        cls.result = bridge.analyze(config)
        cls.patterns = cls.result["face_certificate"]["patterns"]

    def test_all_sixteen_face_patterns_are_certified(self) -> None:
        self.assertEqual(len(self.patterns), 16)
        self.assertTrue(self.result["face_certificate"]["all_local_cases_pass"])
        for row in self.patterns:
            with self.subTest(mask=row["black_mask"]):
                self.assertTrue(row["connectivity_preserved"])
                self.assertTrue(row["embedded_diagonal_gate"])
                self.assertFalse(row["replacement_failures"])

    def test_only_opposite_white_pairs_retain_a_diagonal(self) -> None:
        self.assertEqual(
            self.result["face_certificate"]["retained_diagonal_masks"],
            [5, 10],
        )
        for mask in (5, 10):
            row = self.patterns[mask]
            self.assertEqual(len(row["white_corners"]), 2)
            self.assertEqual(len(row["retained_diagonals"]), 1)

    def test_every_removed_diagonal_has_a_kept_white_path(self) -> None:
        self.assertEqual(
            self.result["face_certificate"]["removed_diagonal_replacement_count"],
            6,
        )
        for row in self.patterns:
            kept = {
                tuple(edge)
                for edge in row["white_boundary_edges"] + row["retained_diagonals"]
            }
            for replacement in row["removed_diagonal_replacements"]:
                path = replacement["boundary_path"]
                self.assertEqual(path[0], replacement["diagonal"][0])
                self.assertEqual(path[-1], replacement["diagonal"][1])
                for first, second in zip(path, path[1:]):
                    self.assertIn(bridge.canonical_edge(first, second), kept)

    def test_surface_proof_chain_names_every_required_isomorphism(self) -> None:
        chain = " ".join(self.result["surface_duality_theorem"]["proof_chain"])
        self.assertIn("long exact sequence", chain)
        self.assertIn("excision", chain)
        self.assertIn("Poincare-Lefschetz", chain)
        self.assertIn("naturality", chain)

    def test_strong_rank_sum_implies_weak_identity(self) -> None:
        rows = self.result["rank_consequences"]
        self.assertEqual(
            [(row["rank_black"], row["rank_white"], row["q"]) for row in rows],
            [(0, 2, -1), (1, 1, 0), (2, 0, 1)],
        )
        self.assertTrue(all(row["strong_residual"] == 0 for row in rows))
        self.assertTrue(all(row["weak_residual"] == 0 for row in rows))

    def test_claim_scope_keeps_degenerate_quotients_separate(self) -> None:
        scope = self.result["lattice_bridge"]["scope"]
        self.assertIn("four distinct corners", scope)
        boundary = " ".join(self.result["scientific_boundary"])
        self.assertIn("separate finite oracle", boundary)


if __name__ == "__main__":
    unittest.main()
