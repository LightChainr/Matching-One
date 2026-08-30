from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import digital_alexander_integral_saturation as theorem  # noqa: E402


class DigitalAlexanderIntegralSaturationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.certificate = theorem.build_certificate()

    def test_integer_white_replacements_preserve_boundary_and_displacement(self) -> None:
        audit = self.certificate["machine_certificates"]["integer_white_face_chains"]
        self.assertEqual(audit["pattern_count"], 16)
        self.assertEqual(audit["replacement_count"], 6)
        self.assertTrue(audit["all_patterns_pass"])
        for row in audit["patterns"]:
            for replacement in row["replacements"]:
                self.assertTrue(replacement["same_integral_boundary"])
                self.assertTrue(replacement["same_integral_lift_displacement"])

    def test_coprime_cover_preserves_every_nontrivial_smith_defect(self) -> None:
        audit = self.certificate["machine_certificates"]["coprime_smith_descent"]
        self.assertTrue(audit["all_rows_pass"])
        self.assertGreater(audit["rank_one_rows"], 0)
        self.assertGreater(audit["rank_two_rows"], 0)
        for row in audit["rows"]:
            self.assertTrue(row["q_coprime_to_d"])
            self.assertTrue(row["H_intersect_qL_equals_qH"])
            self.assertEqual(
                row["index_of_H_intersect_qL_in_qS"],
                row["saturation_defect_d"],
            )
            self.assertGreater(row["saturation_defect_d"], 1)

    def test_rank_one_and_rank_two_examples(self) -> None:
        rank_one = theorem.smith_intersection_row((12,))
        rank_two = theorem.smith_intersection_row((2, 6))
        self.assertEqual(rank_one["rank"], 1)
        self.assertEqual(rank_one["index_of_H_intersect_qL_in_qS"], 12)
        self.assertEqual(rank_two["rank"], 2)
        self.assertEqual(rank_two["index_of_H_intersect_qL_in_qS"], 12)

    def test_theorem_closes_all_rank_boundaries(self) -> None:
        boundaries = self.certificate["degenerate_quotient_descent"]["rank_boundaries"]
        self.assertEqual(set(boundaries), {"rank_zero", "rank_one", "rank_two"})
        self.assertEqual(
            self.certificate["status"],
            "unrestricted_integral_saturation_theorem",
        )
        self.assertTrue(self.certificate["machine_certificates"]["all_pass"])
        self.assertIn(
            "downstairs integral ambient image is exactly H",
            self.certificate["degenerate_quotient_descent"]["component_stabilizer"],
        )

    def test_checked_in_artifacts_reproduce(self) -> None:
        expected_json = json.loads(
            (ROOT / "results/digital-alexander-integral/latest.json")
            .read_text(encoding="utf-8")
        )
        expected_markdown = (
            ROOT / "results/digital-alexander-integral/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(expected_json, self.certificate)
        self.assertEqual(expected_markdown, theorem.render_markdown(self.certificate) + "\n")


if __name__ == "__main__":
    unittest.main()
