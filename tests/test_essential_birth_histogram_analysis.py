
from __future__ import annotations
import json
import sys
import unittest
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import essential_birth_histogram_analysis as analysis  # noqa: E402


class EssentialBirthHistogramAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = analysis.build_artifact()

    def test_joint_reconstructs_marginals_and_moments(self) -> None:
        audit = self.artifact["archive_validation"]
        self.assertTrue(audit["joint_reconstructs_both_marginals"])
        self.assertTrue(audit["moments_match_metadata"])
        self.assertFalse(audit["K_minus_above_K_plus_support_violations"])

    def test_rank_probability_and_mixture_identities_at_half(self) -> None:
        row = self.artifact["exact_at_p_half"]
        f1 = Fraction(row["F_first"])
        f2 = Fraction(row["F_second"])
        self.assertEqual(Fraction(row["P_rank_0"]), 1 - f1)
        self.assertEqual(Fraction(row["P_rank_1"]), f1 - f2)
        self.assertEqual(Fraction(row["P_rank_2"]), f2)
        self.assertEqual(Fraction(row["M"]), f1 + f2 - 1)
        self.assertEqual(Fraction(row["essential_birth_mixture_cdf"]), (f1 + f2) / 2)

    def test_archived_declared_evaluation_reproduces(self) -> None:
        row = self.artifact["archived_declared_point_reproduction"]
        self.assertLess(abs(Decimal(row["M_residual"])), Decimal("1e-45"))
        self.assertLess(abs(Decimal(row["M_prime_residual"])), Decimal("1e-45"))

    def test_archived_root_is_birth_mixture_median(self) -> None:
        root = self.artifact["archived_root_evaluation"]
        self.assertLess(abs(Decimal(root["M"])), Decimal("1e-45"))
        self.assertLess(
            abs(Decimal(root["essential_birth_mixture_cdf"]) - Decimal("0.5")),
            Decimal("1e-45"),
        )
        self.assertLess(abs(Decimal(root["Phi_log_P2_over_P0"])), Decimal("1e-44"))

    def test_neutral_area_is_expected_priority_lifetime(self) -> None:
        moments = self.artifact["priority_center_lifetime_moments"]
        self.assertEqual(moments["E_lifetime"], moments["neutral_area_integral"])
        self.assertEqual(Fraction(moments["E_lifetime"]), Fraction(110129, 1300000))

    def test_missing_marks_are_not_claimed_recoverable(self) -> None:
        missing = self.artifact["recoverability_boundary"]["missing_from_archive"]
        self.assertIn("projective winding line ell", missing)
        self.assertIn("integral saturation index", missing)
        self.assertIn("first/second birth-site local marks", missing)

    def test_checked_in_results_reproduce(self) -> None:
        checked_json = json.loads(
            (ROOT / "results/essential-birth-histogram/latest.json").read_text(encoding="utf-8")
        )
        checked_markdown = (
            ROOT / "results/essential-birth-histogram/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(checked_json, self.artifact)
        self.assertEqual(checked_markdown, analysis.render_markdown(self.artifact))


if __name__ == "__main__":
    unittest.main()
