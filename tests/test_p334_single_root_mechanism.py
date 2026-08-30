from __future__ import annotations

from fractions import Fraction
import json
from math import comb
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p334_single_root_mechanism as mechanism  # noqa: E402


class P334SingleRootMechanismTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = mechanism.build_certificate(ROOT)

    def test_focused_atlas_and_gaussians_are_strict_ulc(self) -> None:
        focused = self.payload["focused_atlas_n13_n17"]
        self.assertEqual(focused["row_count"], 16)
        self.assertTrue(focused["all_strict_single_peak_ulc_ratio"])
        for row in focused["rows"]:
            audit = row["audit"]
            self.assertTrue(audit["strict_single_peak"])
            self.assertTrue(audit["strict_log_concave_on_positive_support"])
            self.assertTrue(audit["adjacent_ratio_strictly_decreasing"])
            self.assertTrue(audit["derivative_coefficient_identity_pass"])

    def test_broad_scan_separates_weak_from_strict_peak(self) -> None:
        broad = self.payload["bounded_scan"]
        self.assertEqual(broad["matrix_count"], 83)
        self.assertEqual(broad["line_sequences"], 240)
        self.assertEqual(broad["orbit_sequences"], 217)
        self.assertEqual(broad["line_pass_counts"]["weak_single_peak"], 240)
        self.assertEqual(broad["orbit_pass_counts"]["weak_single_peak"], 217)
        self.assertEqual(broad["line_pass_counts"]["strict_single_peak"], 228)
        self.assertEqual(broad["orbit_pass_counts"]["strict_single_peak"], 205)
        for label, total in (("line_pass_counts", 240), ("orbit_pass_counts", 217)):
            self.assertEqual(
                broad[label]["strict_log_concave_on_positive_support"], total
            )
            self.assertEqual(
                broad[label]["adjacent_ratio_strictly_decreasing"], total
            )
            self.assertEqual(
                broad[label]["single_nonzero_derivative_sign_change"], total
            )

    def test_minimal_plateau_still_has_one_sign_change(self) -> None:
        row = self.payload["bounded_scan"]["minimal_non_strict_orbit_peak"]
        self.assertEqual(row["matrix"], [[2, 0], [0, 3]])
        self.assertEqual(row["N"], 6)
        self.assertEqual(row["group"], [[1, 0]])
        self.assertEqual(row["audit"]["q_on_support"], ["1/5", "3/5", "3/5"])
        self.assertEqual(row["audit"]["modes"], [3, 4])
        self.assertFalse(row["audit"]["strict_single_peak"])
        self.assertTrue(row["audit"]["weak_single_peak"])
        self.assertTrue(row["audit"]["adjacent_ratio_strictly_decreasing"])

    def test_derivative_coefficient_identity_is_exact(self) -> None:
        counts = [0, 0, 3, 12, 9, 0, 0]
        audit = mechanism.audit_counts(counts, 6)
        self.assertTrue(audit["derivative_coefficient_identity_pass"])
        q = [Fraction(counts[k], comb(6, k)) for k in range(7)]
        for k in range(6):
            left = (k + 1) * counts[k + 1] - (6 - k) * counts[k]
            right = (6 - k) * comb(6, k) * (q[k + 1] - q[k])
            self.assertEqual(left, right)

    def test_matroid_and_nested_upset_obstacles_are_exact(self) -> None:
        witness = self.payload["bounded_scan"]["minimal_direct_rank2_witness"]
        self.assertEqual(witness["matrix"], [[2, 0], [0, 2]])
        self.assertEqual(witness["rank_jump"], 2)
        # U\\V={{1},{2,3,4}} on four elements: q1=q3=1/4, q2=0.
        q = [Fraction(0), Fraction(1, 4), Fraction(0), Fraction(1, 4), Fraction(0)]
        self.assertLess(q[2] * q[2], q[1] * q[3])

    def test_checked_in_artifacts_reproduce(self) -> None:
        expected_json = json.loads(
            (ROOT / "results/p334-single-root-mechanism/latest.json").read_text(
                encoding="utf-8"
            )
        )
        expected_markdown = (
            ROOT / "results/p334-single-root-mechanism/latest.md"
        ).read_text(encoding="utf-8")
        self.assertEqual(expected_json, self.payload)
        self.assertEqual(expected_markdown, mechanism.render_markdown(self.payload) + "\n")


if __name__ == "__main__":
    unittest.main()
