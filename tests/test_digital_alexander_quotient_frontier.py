from __future__ import annotations

import json
import itertools
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import digital_alexander_quotient_frontier as frontier  # noqa: E402


class DigitalAlexanderQuotientFrontierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = frontier.build_artifact()

    def test_all_hnf_representatives_and_permutations_are_exhausted(self) -> None:
        self.assertEqual(self.artifact["status"], "no_counterexample_through_index_8")
        self.assertEqual(self.artifact["HNF_representatives"], 55)
        self.assertEqual(self.artifact["filtration_paths"], 654678)
        self.assertEqual(len(self.artifact["geometries"]), 55)

    def test_face_degeneracy_partition_is_complete(self) -> None:
        self.assertEqual(
            self.artifact["honest_face_representatives"]
            + self.artifact["self_identifying_face_representatives"],
            55,
        )
        self.assertEqual(self.artifact["honest_face_representatives"], 28)
        self.assertEqual(self.artifact["self_identifying_face_representatives"], 27)

    def test_index_seven_extension_is_complete(self) -> None:
        rows = [row for row in self.artifact["geometries"] if row["order"] == 7]
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(row["permutations"] == 5040 for row in rows))
        self.assertEqual(sum(row["four_distinct_face_corners"] for row in rows), 4)

    def test_index_eight_extension_is_complete(self) -> None:
        rows = [row for row in self.artifact["geometries"] if row["order"] == 8]
        self.assertEqual(len(rows), 15)
        self.assertTrue(all(row["permutations"] == 40320 for row in rows))
        self.assertEqual(sum(row["four_distinct_face_corners"] for row in rows), 11)
        self.assertEqual(sum(row["rank_one_plateau_steps"] for row in rows), 1817856)

    def test_cached_and_uncached_paths_are_identical(self) -> None:
        geometry = frontier.integer_torus_geometry(((2, 0), (0, 2)), name="cache-equivalence")
        state_table = frontier.build_state_table(geometry)
        self.assertEqual(len(state_table), 1 << geometry.n)
        for permutation in itertools.permutations(range(geometry.n)):
            self.assertEqual(
                frontier.analyze_permutation(geometry, permutation),
                frontier.analyze_permutation(
                    geometry,
                    permutation,
                    state_table=state_table,
                ),
            )

    def test_all_five_exact_gates_have_zero_failures(self) -> None:
        expected = {
            "birth": 0,
            "line": 0,
            "rank_sum": 0,
            "reconstruction": 0,
            "reflection": 0,
        }
        self.assertEqual(self.artifact["total_failure_counts"], expected)
        self.assertEqual(
            self.artifact["first_counterexamples"],
            {gate: None for gate in expected},
        )

    def test_index_and_plateau_diagnostics_are_frozen(self) -> None:
        geometries = self.artifact["geometries"]
        self.assertEqual(max(row["maximum_saturation_index"] for row in geometries), 1)
        self.assertEqual(sum(row["permutations_with_index_evolution"] for row in geometries), 0)
        self.assertEqual(sum(row["rank_one_plateau_steps"] for row in geometries), 1967984)
        self.assertEqual(sum(row["cached_subsets"] for row in geometries), 5980)
        self.assertTrue(all(row["cached_subsets"] == 1 << row["order"] for row in geometries))

    def test_failure_counts_match_stored_examples(self) -> None:
        for gate, count in self.artifact["total_failure_counts"].items():
            example = self.artifact["first_counterexamples"][gate]
            self.assertEqual(example is None, count == 0)

    def test_each_matrix_has_the_declared_permutation_count(self) -> None:
        for row in self.artifact["geometries"]:
            with self.subTest(matrix=row["matrix"]):
                expected = 1
                for value in range(2, row["order"] + 1):
                    expected *= value
                self.assertEqual(row["permutations"], expected)

    def test_checked_in_results_reproduce(self) -> None:
        checked_json = json.loads(
            (ROOT / "results/digital-alexander-quotient-frontier/latest.json").read_text(encoding="utf-8")
        )
        checked_markdown = (
            ROOT / "results/digital-alexander-quotient-frontier/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(checked_json, self.artifact)
        self.assertEqual(checked_markdown, frontier.render_markdown(self.artifact))


if __name__ == "__main__":
    unittest.main()
