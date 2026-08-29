from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_p200_n650_context_morphism import analyze  # noqa: E402


SOURCE = ROOT / "results" / "local-20260829" / "P200-n650-mixed-join-smoke"
RESULT = ROOT / "results" / "post-reveal-20260829" / "P200-n650-context-morphism"


class P200N650ContextMorphismTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.actual = analyze(
            SOURCE / "n650_20k.batches.csv",
            SOURCE / "n650_20k.metadata.json",
            ROOT / "predictions" / "p200_n650_mixed_join_phaseB_20260829.json",
            SOURCE / "score.json",
        )

    def test_committed_result_recomputes_from_existing_20k_only(self) -> None:
        frozen = json.loads((RESULT / "analysis.json").read_text(encoding="utf-8"))
        self.assertEqual(self.actual, frozen)
        self.assertEqual(frozen["input"]["samples"], 20000)
        self.assertEqual(frozen["status"], "post_reveal_opportunity_map_no_new_production")

    def test_static_interaction_survives_but_geometry_difference_does_not(self) -> None:
        primary = self.actual["primary_partition_residual"]
        self.assertGreater(primary["common_geometry_S"]["chi_square"], 900000)
        self.assertAlmostEqual(
            primary["geometry_difference_D"]["p_value"], 0.26494695308493754
        )
        self.assertEqual(primary["leave_one_sign_stability"], {"ES": 1.0, "OS": 1.0})

    def test_local_baseline_reconstruction_and_density_freeze(self) -> None:
        rows = self.actual["local_incidence_subtraction"]["context_reconstruction"]
        for row in rows:
            self.assertAlmostEqual(
                row["exact_local_incidence_baseline_mean"]
                + row["connected_residual_mean"],
                row["implied_full_mixed_join_mean"],
            )
        freeze = self.actual["morphism_parameter_freeze"]
        self.assertAlmostEqual(
            freeze["common_connected_residual_density_per_source_fiber"]["black"]["estimate"],
            -0.6771680769230769,
        )
        self.assertEqual(
            freeze["conditional_N1300_prediction"]["mean"],
            [-107.41575, 0.0, -68.64794999999998, 0.0],
        )

    def test_unrun_N1300_geometry_pair_has_exact_factor_chain(self) -> None:
        geometries = self.actual["morphism_parameter_freeze"]["conditional_N1300_prediction"]["unrun_geometry_pair"]

        def multiply(first: list[int], second: tuple[int, int]) -> list[int]:
            a, b = first
            c, d = second
            return [a * c - b * d, a * d + b * c]

        for geometry in geometries:
            source = geometry["source_gaussian"]
            self.assertEqual(source[0] ** 2 + source[1] ** 2, 130)
            self.assertEqual(multiply(source, (1, 1)), geometry["N260_after_1_plus_i"])
            self.assertEqual(multiply(source, (2, -1)), geometry["N650_after_2_minus_i"])
            final = geometry["N1300_after_3_plus_i"]
            self.assertEqual(multiply(source, (3, 1)), final)
            self.assertEqual(final[0] ** 2 + final[1] ** 2, 1300)

    def test_typed_H1_is_common_even_and_not_path_evidence(self) -> None:
        ambient = self.actual["typed_ambient_H1"]
        self.assertAlmostEqual(ambient["mean"][0], -1.99905)
        self.assertGreater(ambient["geometry_difference_D"]["p_value"], 0.7)
        boundary = self.actual["identifiability_boundary"]
        self.assertIn("not identified", boundary["path_or_state_memory"])
        self.assertIn("not chronology", ambient["interpretation"])

    def test_scientific_card_is_exactly_five_lines(self) -> None:
        lines = (RESULT / "SCIENTIFIC_CARD.md").read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(lines), 5)
        self.assertEqual(lines, self.actual["scientific_card"])


if __name__ == "__main__":
    unittest.main()
