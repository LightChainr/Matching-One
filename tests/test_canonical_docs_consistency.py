#!/usr/bin/env python3
"""Guard scientific claim boundaries without coupling roadmap prose to them."""

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
        cls.documents = {
            "STATUS": STATUS.read_text(encoding="utf-8"),
            "ROADMAP": ROADMAP.read_text(encoding="utf-8"),
            "SYNTHESIS": SYNTHESIS.read_text(encoding="utf-8"),
        }
        cls.claim_documents = {
            "STATUS": cls.documents["STATUS"],
            "SYNTHESIS": cls.documents["SYNTHESIS"],
        }

    def test_issue43_channel_correction_is_locked_in_claim_documents(self) -> None:
        for name, document in self.claim_documents.items():
            with self.subTest(document=name):
                self.assertIn("DeltaS_cross = -DeltaS_either", document)
                self.assertIn("0.5700315436 / 2", document)

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

    def test_near_term_priorities_are_named_everywhere(self) -> None:
        for name, document in self.documents.items():
            with self.subTest(document=name):
                self.assertIn("#57", document)
                self.assertIn("#50", document)

    def test_n26_beta_result_remains_visible_in_claim_documents(self) -> None:
        for name, document in self.claim_documents.items():
            with self.subTest(document=name):
                self.assertIn("Beta(5,5)", document)
                self.assertIn("Beta(7,7)", document)
                self.assertIn("N=26", document)

        stale_future_phrases = (
            "falsify/extend the exact N=10 `Beta(3,3)` threshold law on N=26 (#115)",
            "Pre-frozen N=26 exact falsification #115",
            "N=10 self-matching `Beta(3,3)` extends to a finite exact family",
        )
        for name, document in self.documents.items():
            for phrase in stale_future_phrases:
                with self.subTest(document=name, phrase=phrase):
                    self.assertNotIn(phrase, document)

    def test_russo_pivotal_progress_is_visible_in_claim_documents(self) -> None:
        for name, document in self.claim_documents.items():
            with self.subTest(document=name):
                self.assertIn("Russo", document)
                self.assertIn("pivotal", document.lower())

    def test_roadmap_is_priority_not_claim_duplication(self) -> None:
        roadmap = self.documents["ROADMAP"]
        self.assertIn("information", roadmap.lower())
        self.assertIn("not a permission", roadmap.lower())
        self.assertNotIn("remain gated behind", roadmap)

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
