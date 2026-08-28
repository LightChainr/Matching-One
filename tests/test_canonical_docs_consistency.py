#!/usr/bin/env python3
"""Guard the canonical claim documents against known protocol regressions."""

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

    def test_issue43_channel_correction_is_canonical_everywhere(self) -> None:
        for name, document in self.documents.items():
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


if __name__ == "__main__":
    unittest.main()
