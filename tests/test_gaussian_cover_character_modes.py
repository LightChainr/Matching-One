from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gaussian_cover_character_modes import render  # noqa: E402


class GaussianCoverCharacterModesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = render()
        cls.artifact = json.loads(
            (
                ROOT
                / "results"
                / "exact-cover-character-oracles"
                / "gaussian_cover_characters.json"
            ).read_text(encoding="utf-8")
        )

    def test_required_smith_groups_are_derived(self) -> None:
        expected = {
            "1+i": ([1, 2], "Z/2"),
            "2+i": ([1, 5], "Z/5"),
            "2-i": ([1, 5], "Z/5"),
            "2i": ([2, 2], "Z/2 x Z/2"),
            "3+i": ([1, 10], "Z/10"),
        }
        for label, (smith, structure) in expected.items():
            group = self.payload["groups"][label]
            self.assertEqual(group["smith_invariants"], smith)
            self.assertEqual(group["additive_group"], structure)
            self.assertTrue(group["exact_checks"]["character_homomorphism"])
            self.assertTrue(group["exact_checks"]["character_orthogonality"])

    def test_norm4_is_noncyclic_and_exponent_two(self) -> None:
        group = self.payload["groups"]["2i"]
        self.assertEqual(group["group_exponent"], 2)
        addition = group["addition_table_indices"]
        for index in range(4):
            self.assertEqual(addition[index][index], 0)
        conclusion = self.payload["norm4_conclusion"]
        self.assertFalse(conclusion["cyclic_Z4"])
        self.assertEqual(conclusion["derived_smith_invariants_for_2i"], [2, 2])

    def test_degree2_squared_character_composition(self) -> None:
        composition = self.payload["degree2_to_degree4_composition"]
        self.assertEqual(composition["identity"], "(1+i)^2=2i")
        self.assertEqual(composition["kernel_element_indices"], [0, 3])
        self.assertEqual(
            composition["hadamard_exponent_table_mod_2"],
            [
                [0, 0, 0, 0],
                [0, 1, 0, 1],
                [0, 0, 1, 1],
                [0, 1, 1, 0],
            ],
        )
        meanings = {
            row["meaning"]: row["stage_frequencies"]
            for row in composition["character_composition"]
        }
        self.assertEqual(meanings["coarse_pullback"], [1, 0])
        self.assertEqual(meanings["new_detail"], [0, 1])
        self.assertEqual(meanings["coarse_times_detail"], [1, 1])

    def test_d4_and_conjugation_actions(self) -> None:
        actions = self.payload["d4_and_conjugation_actions"]
        conjugate = {
            label: next(
                row
                for row in rows
                if row["operation"] == "reflection_after_rotation_0"
            )
            for label, rows in actions.items()
        }
        self.assertEqual(conjugate["2+i"]["target_multiplier"], "2-i")
        self.assertEqual(conjugate["2-i"]["target_multiplier"], "2+i")
        self.assertEqual(conjugate["3+i"]["target_multiplier"], "3-i")
        self.assertEqual(conjugate["1+i"]["target_multiplier"], "1+i")
        self.assertEqual(conjugate["2i"]["element_image_indices"], [0, 1, 2, 3])
        rotation = next(
            row for row in actions["2i"] if row["operation"] == "rotation_90"
        )
        self.assertEqual(rotation["element_image_indices"], [0, 2, 1, 3])
        self.assertEqual(rotation["character_pushforward_indices"], [0, 2, 1, 3])

    def test_frozen_artifact_matches_oracle(self) -> None:
        self.assertEqual(self.artifact, self.payload)


if __name__ == "__main__":
    unittest.main()
