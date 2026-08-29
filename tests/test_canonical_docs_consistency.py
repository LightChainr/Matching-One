#!/usr/bin/env python3
"""Guard claim boundaries without forcing canonical documents to duplicate prose."""

from __future__ import annotations

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
STATUS = ROOT / "docs" / "STATUS.md"
ROADMAP = ROOT / "docs" / "ROADMAP.md"
SYNTHESIS = ROOT / "notes" / "SYNTHESIS-20260828.md"


class CanonicalDocsConsistencyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.status = STATUS.read_text(encoding="utf-8")
        cls.roadmap = ROADMAP.read_text(encoding="utf-8")
        cls.synthesis = SYNTHESIS.read_text(encoding="utf-8")
        cls.documents = {
            "STATUS": cls.status,
            "ROADMAP": cls.roadmap,
            "SYNTHESIS": cls.synthesis,
        }

    def test_issue43_channel_correction_remains_in_claim_ledger(self) -> None:
        self.assertIn("DeltaS_cross = -DeltaS_either", self.status)
        self.assertIn("0.5700315436/2", self.status.replace(" ", ""))

    def test_superseded_even_sector_claims_do_not_reappear(self) -> None:
        stale_claims = (
            "falsified the old positive `P4[S] ~ N^-1` sign assignment",
            "the simple even companion law fails in sign",
            "more N=185/265 replicas for the failed conjunction",
        )
        for name, document in self.documents.items():
            for claim in stale_claims:
                with self.subTest(document=name, claim=claim):
                    self.assertNotIn(claim, document)

    def test_active_execution_points_to_n290(self) -> None:
        self.assertIn("#50", self.status)
        self.assertIn("#50", self.roadmap)
        self.assertIn("#50", self.synthesis)
        self.assertIn("145", self.roadmap)
        self.assertIn("290", self.roadmap)

    def test_norm5_completion_is_visible_without_requiring_issue_number_everywhere(self) -> None:
        for name, document in self.documents.items():
            with self.subTest(document=name):
                self.assertTrue("P57" in document or "#57" in document or "norm-5" in document)
                self.assertIn("H12", document)

    def test_n26_beta_result_remains_visible_in_claim_ledger(self) -> None:
        self.assertIn("Beta(5,5)", self.status)
        self.assertIn("Beta(7,7)", self.status)
        self.assertIn("N=26", self.status)

    def test_russo_pivotal_progress_is_visible_in_claim_ledger(self) -> None:
        self.assertIn("Russo", self.status)
        self.assertIn("pivotal", self.status.lower())

    def test_roadmap_is_priority_not_claim_duplication(self) -> None:
        self.assertIn("information", self.roadmap.lower())
        self.assertIn("not a permission", self.roadmap.lower())
        self.assertNotIn("remain gated behind", self.roadmap)
        self.assertNotIn("do not start", self.roadmap.lower())

    def test_scalar_width_is_not_promoted_after_failed_diagnostics(self) -> None:
        self.assertIn("scalar", self.status.lower())
        self.assertIn("low-rank", self.synthesis.lower())
        self.assertIn("low-rank", self.roadmap.lower())

    def test_norm5_typed_entrypoints_exist(self) -> None:
        required = (
            "scripts/score_norm5_harmonic_primary_typed.py",
            "scripts/score_intrinsic_functional_cocycle_typed.py",
            "predictions/norm5_harmonic_semantic_gate_20260829.yaml",
            "predictions/intrinsic_functional_cocycle_semantic_gate_20260829.yaml",
        )
        for relative in required:
            with self.subTest(path=relative):
                self.assertTrue((ROOT / relative).is_file(), relative)


if __name__ == "__main__":
    unittest.main()
