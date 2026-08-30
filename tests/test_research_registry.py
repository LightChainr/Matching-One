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
                self.assertEqual(row["status"], "canonical_current")
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
        evidence = {row["id"]: row for row in self.ledger["evidence_blocks"]}
        self.assertEqual(len(evidence), len(self.ledger["evidence_blocks"]))
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
        for block in evidence.values():
            for relative in block.get("paths", []):
                with self.subTest(evidence=block["id"], path=relative):
                    self.assertTrue((ROOT / relative).exists(), relative)

    def test_frozen_prediction_paths_are_unique_and_real(self) -> None:
        rows = self.registry["frozen_predictions"]
        paths = [row["path"] for row in rows]
        self.assertEqual(len(paths), len(set(paths)))
        for row in rows:
            with self.subTest(path=row["path"]):
                self.assertEqual(row["status"], "frozen_prediction")
                self.assertTrue((ROOT / row["path"]).is_file(), row["path"])

    def test_evidence_and_topic_paths_exist(self) -> None:
        for section in ("evidence_archives", "topic_derivations", "historical_protocols"):
            for row in self.registry[section]:
                with self.subTest(section=section, path=row["path"]):
                    self.assertTrue((ROOT / row["path"]).exists(), row["path"])

    def test_channel_semantics_minimum_contract(self) -> None:
        semantics = self.registry["channel_semantics"]
        for key in ("implementation", "audit", "primary_norm5_entrypoint", "fullcurve_norm5_entrypoint"):
            with self.subTest(key=key):
                self.assertTrue((ROOT / semantics[key]).is_file(), semantics[key])
        self.assertTrue({"orientation_order", "normalization", "quantity"}.issubset(semantics["required_fields"]))

    def test_frontier_is_disjoint_from_integrated_and_closed_history(self) -> None:
        frontier = [int(row["pr"]) for row in self.registry["frontier_open_as_of_2026_08_29"]]
        self.assertEqual(len(frontier), len(set(frontier)))
        canonical = {int(row["pr"]) for row in self.registry["canonical_integration_history"]}
        manual = {int(row["pr"]) for row in self.registry.get("manual_integrations", [])}
        closed = {int(row["pr"]) for row in self.registry["superseded_active_paths_closed"]}
        self.assertFalse(set(frontier) & canonical)
        self.assertFalse(set(frontier) & manual)
        self.assertFalse(set(frontier) & closed)

    def test_completed_work_is_not_still_queued_as_new_compute(self) -> None:
        new_issues = {int(row["issue"]) for row in self.ledger["new_compute_queue"] if row.get("issue") is not None}
        exact_completed = {int(row["issue"]) for row in self.ledger["completed_exact_tasks"] if row.get("issue") is not None}
        analysis_completed = {int(row["issue"]) for row in self.ledger.get("completed_analysis_tasks", []) if row.get("issue") is not None}
        self.assertFalse(new_issues & exact_completed)
        self.assertFalse(new_issues & analysis_completed)
        self.assertIn(95, analysis_completed)
        self.assertIn(115, exact_completed)

    def test_no_gated_state_remains_in_machine_readable_work_queue(self) -> None:
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
        allowed_new_compute_states = {"active", "ready", "later"}
        self.assertTrue({row["state"] for row in self.ledger["new_compute_queue"]}.issubset(allowed_new_compute_states))


if __name__ == "__main__":
    unittest.main()
