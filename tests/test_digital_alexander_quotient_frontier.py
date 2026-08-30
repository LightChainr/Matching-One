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
        self.assertEqual(self.artifact["status"], "no_counterexample_through_index_12")
        self.assertEqual(self.artifact["HNF_representatives"], 126)
        self.assertEqual(self.artifact["filtration_paths"], 13961736918)
        self.assertEqual(len(self.artifact["geometries"]), 126)

    def test_face_degeneracy_partition_is_complete(self) -> None:
        self.assertEqual(
            self.artifact["honest_face_representatives"]
            + self.artifact["self_identifying_face_representatives"],
            126,
        )
        self.assertEqual(self.artifact["honest_face_representatives"], 83)
        self.assertEqual(self.artifact["self_identifying_face_representatives"], 43)

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

    def test_index_nine_extension_is_complete(self) -> None:
        rows = [row for row in self.artifact["geometries"] if row["order"] == 9]
        self.assertEqual(len(rows), 13)
        self.assertTrue(all(row["permutations"] == 362880 for row in rows))
        self.assertEqual(sum(row["four_distinct_face_corners"] for row in rows), 9)
        self.assertEqual(sum(row["rank_one_plateau_steps"] for row in rows), 15759360)

    def test_index_ten_extension_is_complete(self) -> None:
        rows = [row for row in self.artifact["geometries"] if row["order"] == 10]
        self.assertEqual(len(rows), 18)
        self.assertTrue(all(row["permutations"] == 3628800 for row in rows))
        self.assertEqual(sum(row["four_distinct_face_corners"] for row in rows), 14)
        self.assertEqual(sum(row["rank_one_plateau_steps"] for row in rows), 245174400)

    def test_index_eleven_extension_is_complete(self) -> None:
        rows = [row for row in self.artifact["geometries"] if row["order"] == 11]
        self.assertEqual(len(rows), 12)
        self.assertTrue(all(row["permutations"] == 39916800 for row in rows))
        self.assertEqual(sum(row["four_distinct_face_corners"] for row in rows), 8)
        self.assertEqual(sum(row["rank_one_plateau_steps"] for row in rows), 2167672320)

    def test_index_twelve_extension_is_complete(self) -> None:
        rows = [row for row in self.artifact["geometries"] if row["order"] == 12]
        self.assertEqual(len(rows), 28)
        self.assertTrue(all(row["permutations"] == 479001600 for row in rows))
        self.assertEqual(sum(row["four_distinct_face_corners"] for row in rows), 24)
        self.assertEqual(sum(row["rank_one_plateau_steps"] for row in rows), 51318696960)

    def test_subset_dp_matches_permutation_oracle_through_index_six(self) -> None:
        matrices = frontier.hnf_matrices(maximum_order=6)
        self.assertEqual(len(matrices), 32)
        for matrix in matrices:
            with self.subTest(matrix=matrix):
                self.assertEqual(
                    frontier.summarize_matrix(matrix),
                    frontier.summarize_matrix_by_permutation(matrix),
                )

    def test_subset_dp_failure_and_first_mark_semantics(self) -> None:
        self.assertEqual(frontier._count_paths_avoiding_nodes(3, {1}), 4)
        self.assertEqual(frontier._first_bad_node_permutation(3, {1}), (0, 1, 2))

        marks = [None, (1, 0), (0, 1), (1, 1)]
        failures, first_marks = frontier._marked_path_summary(2, marks, {1})
        self.assertEqual(failures, 2)
        self.assertEqual(first_marks, {"0,1": 1, "1,0": 1})
        self.assertEqual(
            frontier._first_mark_failure_permutation(2, marks, {1}),
            (0, 1),
        )

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
        self.assertEqual(self.artifact["maximum_saturation_index"], 1)
        self.assertEqual(self.artifact["permutations_with_index_evolution"], 0)
        self.assertEqual(self.artifact["rank_one_plateau_steps"], 53749271024)
        self.assertEqual(self.artifact["cached_subsets"], 170332)
        self.assertEqual(max(row["maximum_saturation_index"] for row in geometries), 1)
        self.assertEqual(sum(row["permutations_with_index_evolution"] for row in geometries), 0)
        self.assertEqual(sum(row["rank_one_plateau_steps"] for row in geometries), 53749271024)
        self.assertEqual(sum(row["cached_subsets"] for row in geometries), 170332)
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
