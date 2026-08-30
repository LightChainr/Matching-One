from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import digital_alexander_filtration_oracle as oracle  # noqa: E402


class DigitalAlexanderFiltrationOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = oracle.build_artifact()

    def test_all_declared_permutations_are_exhausted(self) -> None:
        counts = {row["id"]: row["permutations"] for row in self.artifact["geometries"]}
        self.assertEqual(counts, {"axis-L2": 24, "gaussian-2-1": 120})
        self.assertEqual(self.artifact["totals"]["permutations"], 144)

    def test_historical_endpoints_are_direct_rank_births(self) -> None:
        for row in self.artifact["geometries"]:
            self.assertFalse(row["birth_counterexamples"])
            self.assertFalse(row["reconstruction_counterexamples"])
            self.assertGreater(row["positive_rank_one_plateaus"], 0)

    def test_both_endpoint_reflections_hold(self) -> None:
        for row in self.artifact["geometries"]:
            self.assertFalse(row["reflection_counterexamples"])

    def test_rank_sum_and_projective_line_hold_pathwise(self) -> None:
        for row in self.artifact["geometries"]:
            self.assertFalse(row["rank_sum_counterexamples"])
            self.assertFalse(row["line_counterexamples"])
            self.assertGreater(row["rank_one_plateau_steps"], 0)
            self.assertGreaterEqual(row["maximum_saturation_index"], 1)

    def test_one_permutation_reconstructs_every_rank(self) -> None:
        geometry = oracle.gaussian_integer_torus(2, 1)
        permutation = (0, 1, 2, 3, 4)
        result = oracle.analyze_permutation(geometry, permutation)
        self.assertEqual(result["k_first_direct"], result["k_minus_historical"])
        self.assertEqual(result["k_second_direct"], result["k_plus_historical"])
        self.assertEqual(result["endpoint_reflection_residuals"], [0, 0])
        self.assertFalse(result["reconstruction_failures"])
        self.assertFalse(result["line_failures"])

    def test_checked_in_results_reproduce(self) -> None:
        checked_json = json.loads(
            (ROOT / "results/digital-alexander-filtration/latest.json").read_text(encoding="utf-8")
        )
        checked_markdown = (
            ROOT / "results/digital-alexander-filtration/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(checked_json, self.artifact)
        self.assertEqual(checked_markdown, oracle.render_markdown(self.artifact))


if __name__ == "__main__":
    unittest.main()
