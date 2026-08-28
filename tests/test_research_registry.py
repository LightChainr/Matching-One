#!/usr/bin/env python3
"""Validate the repository research-map and artifact-registry contracts."""

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
            relative = row["path"]
            with self.subTest(path=relative):
                self.assertEqual(row["status"], "canonical_current")
                self.assertTrue((ROOT / relative).exists(), relative)

    def test_research_track_ids_are_unique(self) -> None:
        tracks = self.ledger["research_tracks"]
        ids = [row["id"] for row in tracks]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(ids, list("ABCDEFGHI"))

    def test_evidence_references_are_declared_and_paths_exist(self) -> None:
        blocks = self.ledger["evidence_blocks"]
        evidence_ids = {row["id"] for row in blocks}
        self.assertEqual(len(evidence_ids), len(blocks))

        required = {
            "E_CHANNEL_SEMANTICS",
            "E_RELIABILITY_SIGNATURE",
            "E_RUSSO_PIVOTAL",
            "E_SELFMATCHING_TANGENT",
            "E_P48_NEW_GEOMETRY",
            "E_NORM5_SCORE_SEMANTICS",
            "E_N26_BETA_NEGATIVE",
        }
        self.assertTrue(required.issubset(evidence_ids), required - evidence_ids)

        for track in self.ledger["research_tracks"]:
            for evidence_id in track.get("evidence", []):
                with self.subTest(track=track["id"], evidence=evidence_id):
                    self.assertIn(evidence_id, evidence_ids)

        for block in blocks:
            for relative in block.get("paths", []):
                with self.subTest(evidence=block["id"], path=relative):
                    self.assertTrue((ROOT / relative).exists(), relative)

    def test_frozen_prediction_registry_points_to_real_files(self) -> None:
        rows = self.registry["frozen_predictions"]
        paths = [row["path"] for row in rows]
        self.assertEqual(len(paths), len(set(paths)))
        for row in rows:
            with self.subTest(path=row["path"]):
                self.assertEqual(row["status"], "frozen_prediction")
                self.assertTrue((ROOT / row["path"]).is_file(), row["path"])

    def test_evidence_archives_are_not_declared_canonical_documents(self) -> None:
        canonical_paths = {
            row["path"] for row in self.registry["canonical_documents"]
        }
        for row in self.registry["evidence_archives"]:
            with self.subTest(path=row["path"]):
                self.assertNotIn(row["path"], canonical_paths)
                self.assertTrue((ROOT / row["path"]).exists(), row["path"])

    def test_channel_semantics_contract_is_registered(self) -> None:
        semantics = self.registry["channel_semantics"]
        required_paths = (
            semantics["implementation"],
            semantics["audit"],
            semantics["primary_norm5_entrypoint"],
            semantics["fullcurve_norm5_entrypoint"],
        )
        for relative in required_paths:
            with self.subTest(path=relative):
                self.assertTrue((ROOT / relative).is_file(), relative)
        self.assertIn("orientation_order", semantics["required_fields"])
        self.assertIn("normalization", semantics["required_fields"])

    def test_completed_n26_beta_is_not_still_new_compute(self) -> None:
        new_issue_numbers = {
            int(row["issue"])
            for row in self.ledger["new_compute_queue"]
            if row.get("issue") is not None
        }
        self.assertNotIn(115, new_issue_numbers)
        completed = {
            int(row["issue"]): row
            for row in self.ledger["completed_exact_tasks"]
            if row.get("issue") is not None
        }
        self.assertIn(115, completed)
        self.assertEqual(
            completed[115]["outcome"], "both_frozen_Beta5_and_Beta7_laws_fail"
        )

    def test_frontier_excludes_closed_superseded_paths(self) -> None:
        rows = self.registry["frontier_open_as_of_2026_08_29"]
        numbers = [int(row["pr"]) for row in rows]
        self.assertEqual(len(numbers), len(set(numbers)))
        superseded = {
            168,
            128,
            135,
            137,
            110,
            163,
            157,
            169,
            148,
            167,
            152,
            162,
        }
        self.assertFalse(superseded.intersection(numbers))

        recorded_superseded = {
            int(row["pr"]) for row in self.registry["superseded_active_paths_closed"]
        }
        self.assertTrue(superseded.issubset(recorded_superseded))

    def test_canonical_integration_history_is_unique(self) -> None:
        rows = self.registry["canonical_integration_history"]
        numbers = [int(row["pr"]) for row in rows]
        self.assertEqual(len(numbers), len(set(numbers)))
        self.assertTrue({170, 171, 172, 173, 174, 175, 176, 177}.issubset(numbers))


if __name__ == "__main__":
    unittest.main()
