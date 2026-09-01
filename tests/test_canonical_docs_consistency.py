#!/usr/bin/env python3
"""Guard claim boundaries and current execution direction without duplicating prose."""

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
        cls.documents = {"STATUS": cls.status, "ROADMAP": cls.roadmap, "SYNTHESIS": cls.synthesis}

    def test_issue43_channel_correction_remains_in_claim_ledger(self) -> None:
        self.assertIn("DeltaS_cross = -DeltaS_either", self.status)
        self.assertIn("0.5700315436/2", self.status.replace(" ", ""))

    def test_superseded_even_sector_claims_do_not_reappear(self) -> None:
        stale = (
            "falsified the old positive `P4[S] ~ N^-1` sign assignment",
            "the simple even companion law fails in sign",
            "more N=185/265 replicas for the failed conjunction",
        )
        for name, document in self.documents.items():
            for claim in stale:
                with self.subTest(document=name, claim=claim):
                    self.assertNotIn(claim, document)

    def test_n290_is_completed_not_active_compute(self) -> None:
        self.assertIn("N145->290", self.status)
        self.assertIn("#50 N145->290 full curve", self.roadmap)
        self.assertIn("complete", self.roadmap.lower())
        active_prefix = self.roadmap.split("## Ready", 1)[0]
        self.assertNotIn("### 1. N=145", active_prefix)

    def test_p50_nullspace_amendment_is_visible_without_new_queue(self) -> None:
        self.assertIn("cutoff-sensitive", self.status)
        self.assertIn("#543", self.status)
        self.assertIn("#543", self.roadmap)
        self.assertIn("complete and closed", self.roadmap)

    def test_current_discriminators_are_visible(self) -> None:
        for needle in ("coalescence", "modulus", "Norm-4", "local pivotal"):
            with self.subTest(needle=needle):
                self.assertIn(needle.lower(), self.roadmap.lower())
        self.assertIn("x≈4", self.synthesis)
        self.assertIn("x=21/4", self.synthesis)

    def test_norm5_and_independent_signal_boundary_remain_visible(self) -> None:
        self.assertIn("H12", self.status)
        self.assertIn("global zero", self.status.lower())
        self.assertIn("child block alone", self.status.lower())

    def test_n26_beta_result_remains_visible(self) -> None:
        self.assertIn("Beta(5,5)", self.status)
        self.assertIn("Beta(7,7)", self.status)
        self.assertIn("N=26", self.status)

    def test_russo_pivotal_progress_is_visible(self) -> None:
        self.assertIn("Russo", self.status)
        self.assertIn("pivotal", self.status.lower())

    def test_roadmap_is_priority_not_permission_system(self) -> None:
        self.assertIn("information", self.roadmap.lower())
        self.assertIn("not a permission", self.roadmap.lower())
        self.assertNotIn("remain gated behind", self.roadmap)
        self.assertNotIn("do not start", self.roadmap.lower())

    def test_failed_scalar_shortcuts_stay_retired(self) -> None:
        self.assertIn("scalar width", self.status.lower())
        self.assertIn("low-rank", self.synthesis.lower())
        self.assertIn("free exponent", self.roadmap.lower())

    def test_same_readout_bruteforce_is_not_recommended(self) -> None:
        self.assertIn("N130/N170", self.roadmap)
        self.assertIn("stop adding replicas", self.roadmap.lower())
        self.assertIn("third primitive norm-2", self.roadmap.lower())

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
