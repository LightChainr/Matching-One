#!/usr/bin/env python3
"""Guard claim boundaries and current execution direction without duplicating prose."""

from __future__ import annotations

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
STATUS = ROOT / "docs" / "STATUS.md"
ROADMAP = ROOT / "docs" / "ROADMAP.md"
RESEARCH_MAP = ROOT / "docs" / "RESEARCH-MAP.md"
NEXT_TARGETS = ROOT / "docs" / "NEXT-TARGETS.md"
SYNTHESIS = ROOT / "notes" / "SYNTHESIS-20260828.md"


class CanonicalDocsConsistencyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.readme = README.read_text(encoding="utf-8")
        cls.status = STATUS.read_text(encoding="utf-8")
        cls.roadmap = ROADMAP.read_text(encoding="utf-8")
        cls.research_map = RESEARCH_MAP.read_text(encoding="utf-8")
        cls.next_targets = NEXT_TARGETS.read_text(encoding="utf-8")
        cls.synthesis = SYNTHESIS.read_text(encoding="utf-8")
        cls.documents = {
            "README": cls.readme,
            "STATUS": cls.status,
            "ROADMAP": cls.roadmap,
            "RESEARCH_MAP": cls.research_map,
            "NEXT_TARGETS": cls.next_targets,
            "SYNTHESIS": cls.synthesis,
        }

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

    def test_five_scientific_coordinates_are_visible(self) -> None:
        combined = (self.readme + self.research_map + self.next_targets).lower()
        for coordinate in ("state", "source", "observer", "geometry", "acquisition"):
            with self.subTest(coordinate=coordinate):
                self.assertIn(coordinate, combined)

    def test_frontier_lifecycle_is_distinct_from_main(self) -> None:
        for lifecycle in ("main_integrated", "open_pr", "branch_only", "hypothesis"):
            with self.subTest(lifecycle=lifecycle):
                self.assertIn(lifecycle, self.status)
                self.assertIn(lifecycle, self.next_targets)

    def test_two_activation_and_observer_sector_tension_are_visible(self) -> None:
        combined = self.readme + self.status + self.roadmap + self.next_targets
        self.assertIn("K1=K_minus", combined.replace(" ", ""))
        self.assertIn("K2=K_plus", combined.replace(" ", ""))
        self.assertIn("71/21/8", combined.replace(" ", ""))
        self.assertIn("charged", combined.lower())
        self.assertIn("global", combined.lower())
        self.assertIn("ten", combined.lower())
        self.assertIn("K1", combined)
        self.assertIn("N=265/325/425", combined)

    def test_new_attention_order_is_decision_output_driven(self) -> None:
        headings = [line for line in self.roadmap.splitlines() if line.startswith("### ")]
        numbered = [line for line in headings if line[4:5].isdigit()]
        self.assertGreaterEqual(len(numbered), 7)
        self.assertIn("K1/K2", self.roadmap)
        self.assertIn("connectivity/defect radical", self.roadmap)
        self.assertIn("Next decision output", self.roadmap)


if __name__ == "__main__":
    unittest.main()
