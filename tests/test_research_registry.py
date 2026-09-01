#!/usr/bin/env python3
"""Validate lightweight research navigation without recreating process gates."""

from __future__ import annotations

from pathlib import Path
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
LEDGER_PATH = ROOT / "analysis" / "research_ledger.yaml"
REGISTRY_PATH = ROOT / "analysis" / "artifact_registry.yaml"


def load_yaml(path: Path) -> dict:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError(f"{path} must contain a YAML mapping")
    return payload


class ResearchRegistryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.ledger = load_yaml(LEDGER_PATH)
        cls.registry = load_yaml(REGISTRY_PATH)

    def test_canonical_views_exist(self) -> None:
        for label, relative in self.ledger["canonical_views"].items():
            with self.subTest(label=label, path=relative):
                self.assertTrue((ROOT / relative).exists(), relative)
        for row in self.registry["canonical_documents"]:
            with self.subTest(path=row["path"]):
                self.assertIn(row["status"], {"canonical_current", "canonical_audit"})
                self.assertTrue((ROOT / row["path"]).exists(), row["path"])

    def test_only_three_hard_constraints_are_registered(self) -> None:
        expected = {
            "preserve_frozen_predictions_and_committed_result_history",
            "require_identical_observable_semantics_or_a_registered_exact_map_for_claim_bearing_scores",
            "do_not_add_correlated_views_of_one_raw_random_block_as_independent_primary_evidence",
        }
        self.assertEqual(set(self.ledger["hard_constraints"]), expected)
        self.assertFalse(self.registry["integration_policy"]["registry_conflicts_block_science"])

    def test_research_track_ids_and_evidence_references(self) -> None:
        tracks = self.ledger["research_tracks"]
        ids = [row["id"] for row in tracks]
        self.assertEqual(ids, list("ABCDEFGHI"))
        all_blocks = (
            self.ledger["evidence_blocks"]
            + self.ledger["decision_evidence_blocks"]
            + self.ledger["unmerged_candidate_support"]
        )
        evidence = {row["id"]: row for row in all_blocks}
        self.assertEqual(len(evidence), len(all_blocks))
        required = {
            "E_CHANNEL_SEMANTICS",
            "E_RUSSO_PIVOTAL",
            "E_PREQUENTIAL",
            "E_KRAWTCHOUK_MODES",
            "E_RECTANGULAR_Q4",
            "E_N26_BETA_NEGATIVE",
        }
        self.assertTrue(required.issubset(evidence), required - set(evidence))
        for track in tracks:
            for evidence_id in track.get("evidence", []):
                with self.subTest(track=track["id"], evidence=evidence_id):
                    self.assertIn(evidence_id, evidence)
            for candidate_id in track.get("unmerged_candidate_support", []):
                with self.subTest(track=track["id"], candidate=candidate_id):
                    self.assertIn(candidate_id, evidence)
        for block in self.ledger["evidence_blocks"]:
            for relative in block.get("paths", []):
                with self.subTest(evidence=block["id"], path=relative):
                    self.assertTrue((ROOT / relative).exists(), relative)

        for block in self.ledger["decision_evidence_blocks"]:
            with self.subTest(decision=block["id"]):
                self.assertRegex(block["sha"], r"^[0-9a-f]{40}$")
                self.assertIn("independent_primary", block)

        decision_ids = {row["id"] for row in self.ledger["decision_evidence_blocks"]}
        candidate_prs = []
        for block in self.ledger["unmerged_candidate_support"]:
            with self.subTest(candidate=block["id"]):
                self.assertNotIn(block["id"], decision_ids)
                self.assertRegex(block["sha"], r"^[0-9a-f]{40}$")
                self.assertFalse(block["independent_primary"])
                candidate_prs.append(int(block["pr"]))
        self.assertEqual(len(candidate_prs), len(set(candidate_prs)))

    def test_frozen_prediction_paths_are_unique_and_real(self) -> None:
        rows = self.registry["frozen_predictions"]
        paths = [row["path"] for row in rows]
        self.assertEqual(len(paths), len(set(paths)))
        for row in rows:
            with self.subTest(path=row["path"]):
                self.assertEqual(row["status"], "frozen_prediction")
                self.assertTrue((ROOT / row["path"]).is_file(), row["path"])

    def test_evidence_and_topic_paths_exist(self) -> None:
        for section in ("evidence_archives", "historical_protocols"):
            for row in self.registry[section]:
                with self.subTest(section=section, path=row["path"]):
                    self.assertTrue((ROOT / row["path"]).exists(), row["path"])

    def test_channel_semantics_minimum_contract(self) -> None:
        semantics = self.registry["channel_semantics"]
        for key in ("implementation", "audit", "primary_norm5_entrypoint", "fullcurve_norm5_entrypoint"):
            with self.subTest(key=key):
                self.assertTrue((ROOT / semantics[key]).is_file(), semantics[key])
        self.assertTrue({"orientation_order", "normalization", "quantity"}.issubset(semantics["required_fields"]))

    def test_unmerged_assets_are_unique_and_not_integrated_history(self) -> None:
        unmerged = [int(row["pr"]) for row in self.registry["unmerged_assets"]]
        self.assertEqual(len(unmerged), len(set(unmerged)))
        absorbed = {int(row["pr"]) for row in self.registry["absorbed_docs_prs"]}
        merged = {int(row["pr"]) for row in self.registry["same_day_exact_pipeline_history"]}
        self.assertFalse(set(unmerged) & absorbed)
        self.assertFalse(set(unmerged) & merged)

    def test_completed_work_is_not_still_queued_as_new_compute(self) -> None:
        active_issues = {int(row["issue"]) for row in self.ledger["active_execution"] if row.get("issue") is not None}
        exact_completed = {int(row["issue"]) for row in self.ledger["completed_exact_tasks"] if row.get("issue") is not None}
        analysis_completed = {int(row["issue"]) for row in self.ledger.get("completed_analysis_tasks", []) if row.get("issue") is not None}
        self.assertFalse(active_issues & exact_completed)
        self.assertFalse(active_issues & analysis_completed)
        self.assertIn(95, analysis_completed)
        self.assertIn(115, exact_completed)

    def test_active_queue_is_explicit_and_has_no_gated_state(self) -> None:
        def walk(value):
            if isinstance(value, dict):
                for key, item in value.items():
                    if key == "state":
                        yield item
                    yield from walk(item)
            elif isinstance(value, list):
                for item in value:
                    yield from walk(item)

        self.assertNotIn("gated", set(walk(self.ledger)))
        self.assertEqual(
            self.ledger["active_execution"],
            [
                {
                    "issue": 537,
                    "parent_issue": 337,
                    "kind": "theory",
                    "target": "freeze_formal_ordinary_no_extra_landing_contract_then_complete_surviving_signed_rate",
                    "status": "active_no_sampling",
                    "random_sample_budget": 0,
                    "next_falsifier": "formal_contract_includes_or_excludes_provisional_clean_two_bridge_counterexample",
                    "stop_rule": "no_N_random_or_five_six_arm_work_before_one_semantic_inclusion_decision",
                }
            ],
        )


if __name__ == "__main__":
    unittest.main()
