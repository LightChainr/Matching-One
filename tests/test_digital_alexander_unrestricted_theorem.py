from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import digital_alexander_unrestricted_theorem as theorem  # noqa: E402


class DigitalAlexanderUnrestrictedTheoremTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        config = json.loads(
            (ROOT / "analysis/digital_alexander_unrestricted_manifest.json")
            .read_text(encoding="utf-8")
        )
        cls.artifact = theorem.build_artifact(config)

    def test_all_universal_face_chain_replacements_are_exact(self) -> None:
        face = self.artifact["machine_certificates"]["universal_face_chains"]
        self.assertEqual(face["pattern_count"], 16)
        self.assertEqual(face["retained_diagonal_masks"], [5, 10])
        self.assertEqual(face["removed_diagonal_count"], 6)
        self.assertTrue(face["all_patterns_pass"])
        for row in face["patterns"]:
            for replacement in row["replacement_certificates"]:
                self.assertTrue(replacement["same_relative_boundary"])
                self.assertEqual(replacement["difference_lift_displacement"], [0, 0])

    def test_projection_certificate_exercises_degenerate_incidence(self) -> None:
        projection = self.artifact["machine_certificates"]["degenerate_quotient_projection"]
        self.assertEqual(projection["representatives"], 86)
        self.assertEqual(projection["honest_face_representatives"], 51)
        self.assertEqual(projection["self_identifying_face_representatives"], 35)
        self.assertGreater(projection["representatives_with_primal_loops"], 0)
        self.assertGreater(projection["representatives_with_matching_loops"], 0)
        self.assertGreater(projection["representatives_with_repeated_matching_endpoints"], 0)
        self.assertEqual(projection["failure_count"], 0)
        self.assertIsNone(projection["first_failure"])

    def test_twice_period_cover_is_four_sheeted_and_honest(self) -> None:
        cover = self.artifact["machine_certificates"]["finite_honest_cover"]
        self.assertEqual(cover["representatives"], 86)
        self.assertEqual(cover["cover_degree"], 4)
        self.assertEqual(cover["honest_cover_representatives"], 86)
        self.assertEqual(cover["failure_count"], 0)
        self.assertIsNone(cover["first_failure"])
        for row in cover["rows"]:
            self.assertEqual(row["cover_degree"], 4)
            self.assertEqual(row["H1_map_in_period_bases"], [[2, 0], [0, 2]])
            self.assertEqual(row["intersection_scale"], 4)

    def test_cover_lemma_is_basis_and_orientation_independent(self) -> None:
        matrices = (
            ((1, 0), (0, 1)),
            ((0, 1), (-1, 0)),
            ((2, 3), (5, 7)),
            ((-3, 2), (1, 4)),
        )
        for matrix in matrices:
            with self.subTest(matrix=matrix):
                row = theorem.finite_regular_cover_row(matrix)
                self.assertEqual(row["cover_degree"], 4)
                self.assertTrue(row["cover_has_four_distinct_face_corners"])
                self.assertEqual(row["H1_map_in_period_bases"], [[2, 0], [0, 2]])
                self.assertEqual(row["intersection_scale"], 4)

    def test_cached_subset_regression_includes_all_index_ten_states(self) -> None:
        states = self.artifact["machine_certificates"]["cached_subset_regression"]
        self.assertEqual(states["representatives"], 86)
        self.assertEqual(states["states"], 31068)
        self.assertEqual(states["rank_one_states"], 17248)
        for gate in (
            "rank_sum_failures", "rank_mark_failures",
            "primitive_line_failures", "nonsaturated_rank_one_states",
        ):
            self.assertEqual(states[gate], 0)
        self.assertIsNone(states["first_failure"])

    def test_symplectic_and_filtration_consequences_are_symbolic(self) -> None:
        machine = self.artifact["machine_certificates"]
        self.assertGreater(machine["symplectic_line"]["orthogonal_nonzero_pairs_checked"], 0)
        self.assertEqual(machine["symplectic_line"]["failure_count"], 0)
        self.assertEqual(machine["filtration_algebra"]["threshold_pairs_checked"], 364)
        self.assertEqual(machine["filtration_algebra"]["failure_count"], 0)

    def test_theorem_has_no_honest_quotient_cell_hypothesis(self) -> None:
        result = self.artifact["theorem"]
        self.assertEqual(
            result["period_lattice"],
            "L=P Z^2 for every nonsingular integer 2x2 matrix P",
        )
        self.assertIsNone(result["extra_honest_cell_hypothesis"])
        self.assertTrue(self.artifact["machine_certificates"]["all_pass"])

    def test_checked_in_artifacts_reproduce(self) -> None:
        expected_json = json.loads(
            (ROOT / "results/digital-alexander-unrestricted/latest.json")
            .read_text(encoding="utf-8")
        )
        expected_markdown = (
            ROOT / "results/digital-alexander-unrestricted/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(expected_json, self.artifact)
        self.assertEqual(expected_markdown, theorem.render_markdown(self.artifact))


if __name__ == "__main__":
    unittest.main()
