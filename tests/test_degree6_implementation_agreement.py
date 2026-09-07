#!/usr/bin/env python3
"""Locks the two-implementation cross-check of the C(1..6, 3) and C(1..6, 4) censuses."""

from __future__ import annotations

import copy
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import degree6_implementation_agreement as agreement  # noqa: E402
import degree6_independent_replication as replication  # noqa: E402

ARTIFACT = ROOT / "results" / "pslq-degree6-implementation-agreement" / "latest.json"
INTERVALS = (
    "jacobsen-2015-eigenvalue",
    "mertens-2022-p-med",
    "mertens-2022-p-cell",
    "yang-zhou-2024-corrected",
)
# The narrowest interval, and the one carrying Result A.  One interval is
# recensused here rather than four: the census is ~12 s per interval of pure
# Python rational arithmetic, and the other three are covered by the agreement
# check against the primary implementation, which CI rebuilds in full.
RECENSUSED = "jacobsen-2015-eigenvalue"
# The global minimiser over the whole class differs between the two heights.
CLOSEST = {3: [0, -2, 2, 2, -1, 2, 1], 4: [1, -3, -1, 2, 2, 4, 4]}
CLASS_TOTAL = {3: 409_584, 4: 2_351_328}


def _rows_at(rows, height):
    return [row for row in rows if row["coefficient_height_max"] == height]


class ImplementationAgreementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.committed = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        cls.built = agreement.build_result()

    def test_committed_artifact_reproduces(self) -> None:
        self.assertEqual(self.committed, self.built)
        report = agreement.validate_result(self.committed)
        self.assertEqual(report["status"], "valid")
        self.assertTrue(report["implementations_agree"])

    def test_every_cell_agrees_at_both_heights(self) -> None:
        self.assertEqual(self.built["heights_compared"], [3, 4])
        self.assertEqual(self.built["cells_compared"], 48)
        self.assertEqual(self.built["cells_in_agreement"], 48)
        # 4 intervals x 2 heights, each interval by height with all six degrees.
        self.assertEqual(len(self.built["intervals"]), 8)
        for height in (3, 4):
            rows = _rows_at(self.built["intervals"], height)
            self.assertEqual([row["interval_id"] for row in rows], list(INTERVALS))
            for row in rows:
                self.assertTrue(row["exclusion_verdict"], row["interval_id"])
                self.assertTrue(row["exclusion_verdicts_agree"], row["interval_id"])
                for cell in row["by_degree"]:
                    self.assertEqual(cell["disagreements"], [], (height, row["interval_id"], cell["degree"]))
                    self.assertEqual(cell["values"]["screen_survivors"], 0)
                    self.assertEqual(cell["values"]["root_containing_polynomials"], 0)
            total = sum(
                cell["values"]["polynomials_in_class"]
                for cell in rows[0]["by_degree"]
            )
            self.assertEqual(total, CLASS_TOTAL[height])

    def test_closest_member_matches_and_respects_the_mean_value_bound(self) -> None:
        """Equal counts could be coincidence; equal residuals cannot.

        The two implementations evaluate the minimiser at different points, so
        their residuals must differ -- but by no more than ``D*(u-l)/2``.
        """
        for row in self.built["intervals"]:
            height = row["coefficient_height_max"]
            closest = row["closest_member"]
            self.assertTrue(closest["coefficients_agree"], row["interval_id"])
            self.assertEqual(
                closest["coefficients_ascending"], CLOSEST[height], (height, row["interval_id"])
            )
            gap = Fraction(closest["residual_gap_text"])
            allowance = Fraction(closest["mean_value_allowance_text"])
            self.assertGreater(gap, 0, row["interval_id"])
            self.assertLessEqual(gap, allowance, row["interval_id"])
            self.assertTrue(closest["within_mean_value_bound"], row["interval_id"])

    def test_derivative_bound_is_the_polynomial_specific_one(self) -> None:
        # x^6+2x^5-x^4+2x^3+2x^2-2x: 1*2+2*2+3*2+4*1+5*2+6*1 = 32, well under the
        # class-wide bound of 3*6*7/2 = 63 that the censuses screen against.
        self.assertEqual(agreement.derivative_bound([0, -2, 2, 2, -1, 2, 1]), 32)
        self.assertLess(32, 3 * 6 * 7 // 2)
        self.assertEqual(agreement.derivative_bound([1]), 0)
        for row in self.built["intervals"]:
            closest = row["closest_member"]
            self.assertEqual(
                closest["polynomial_derivative_bound"],
                agreement.derivative_bound(closest["coefficients_ascending"]),
            )

    def test_the_shared_sturm_path_is_declared_and_did_not_run(self) -> None:
        """The replication is partial and the artifact has to say so."""
        shared = self.built["shared_code"]
        self.assertEqual(shared["path"], "scripts/exact_polynomial_root_certificate.py")
        self.assertIn("never runs", shared["contribution_to_this_result"])
        self.assertIn("degree-4 height-100 census", self.built["claim_boundary"]["excluded"])
        paths = {row["path"] for row in self.built["implementations"]}
        self.assertEqual(
            paths,
            {
                "scripts/degree6_low_height_exclusion.py",
                "scripts/degree6_independent_replication.py",
            },
        )
        screens = {row["screen"] for row in self.built["implementations"]}
        self.assertEqual(len(screens), 2)

    def test_replication_artifact_regenerates(self) -> None:
        """One interval is recensused so the second implementation is not a fossil."""
        for height in (3, 4):
            prefix = "low-height" if height == 3 else "height4"
            committed = json.loads(
                (
                    ROOT / "results" / f"pslq-degree6-{prefix}-replication-{RECENSUSED}" / "latest.json"
                ).read_text(encoding="utf-8")
            )
            self.assertEqual(replication.build_result(RECENSUSED, height=height), committed)
            self.assertTrue(committed["interval_result"]["excluded"])
            self.assertEqual(committed["interval_result"]["class_size_total"], CLASS_TOTAL[height])

    def test_a_planted_disagreement_is_reported(self) -> None:
        """Guards the comparison itself: agreement must be earned, not default."""
        original = agreement._load

        def patched(directory: str) -> dict:
            payload = copy.deepcopy(original(directory))
            if directory.startswith("pslq-degree6-low-height-replication-"):
                payload["interval_result"]["by_degree"][3]["screen_survivors"] += 1
            return payload

        agreement._load = patched
        try:
            tampered = agreement.build_result()
        finally:
            agreement._load = original
        self.assertFalse(tampered["implementations_agree"])
        # 48 cells, four height-3 intervals each losing one degree-4 cell.
        self.assertEqual(tampered["cells_in_agreement"], 44)
        for row in _rows_at(tampered["intervals"], 3):
            degree_four = next(c for c in row["by_degree"] if c["degree"] == 4)
            self.assertEqual(degree_four["disagreements"], ["screen_survivors"])
        for row in _rows_at(tampered["intervals"], 4):
            for cell in row["by_degree"]:
                self.assertEqual(cell["disagreements"], [])

    def test_a_mismatched_contract_is_refused(self) -> None:
        original = agreement._load

        def patched(directory: str) -> dict:
            payload = copy.deepcopy(original(directory))
            payload["contract_sha256"] = "0" * 64
            return payload

        agreement._load = patched
        try:
            with self.assertRaises(ValueError):
                agreement.build_result()
        finally:
            agreement._load = original


if __name__ == "__main__":
    unittest.main()
