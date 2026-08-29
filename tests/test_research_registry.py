#!/usr/bin/env python3
"""Validate lightweight research navigation without recreating process gates."""

from __future__ import annotations

from pathlib import Path
from pathlib import PurePosixPath
import re
import subprocess
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
LEDGER_PATH = ROOT / "analysis" / "research_ledger.yaml"
REGISTRY_PATH = ROOT / "analysis" / "artifact_registry.yaml"
TWO_ACTIVATION_MANIFEST = ROOT / "analysis" / "two_activation_h4_manifest.yaml"
HEX40 = re.compile(r"^[0-9a-f]{40}$")


def load_yaml(path: Path) -> dict:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError(f"{path} must contain a YAML mapping")
    return payload


def assert_repo_relative(testcase: unittest.TestCase, relative: str) -> None:
    path = PurePosixPath(relative)
    testcase.assertFalse(path.is_absolute(), relative)
    testcase.assertNotIn("..", path.parts, relative)
    testcase.assertNotEqual(relative, "", relative)


def git_object_exists(commit: str, relative: str) -> bool:
    result = subprocess.run(
        ["git", "cat-file", "-e", f"{commit}:{relative}"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    return result.returncode == 0


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
        rows = self.registry["frontier_open_as_of_2026_08_29"]
        frontier = [int(row["pr"]) for row in rows]
        self.assertEqual(len(frontier), len(set(frontier)))
        self.assertTrue(
            {84, 196, 197, 228, 229, 230, 245, 246, 247, 267, 273, 277, 281}.issubset(frontier)
        )
        canonical = {int(row["pr"]) for row in self.registry["canonical_integration_history"]}
        manual = {int(row["pr"]) for row in self.registry.get("manual_integrations", [])}
        closed = {int(row["pr"]) for row in self.registry["superseded_active_paths_closed"]}
        self.assertFalse(set(frontier) & canonical)
        self.assertFalse(set(frontier) & manual)
        self.assertFalse(set(frontier) & closed)
        for row in rows:
            with self.subTest(pr=row["pr"]):
                self.assertEqual(row["state"], "open")
                self.assertRegex(row["commit"], HEX40)
                assert_repo_relative(self, row["representative_path"])
                self.assertTrue(
                    git_object_exists(row["commit"], row["representative_path"]),
                    f"PR #{row['pr']} {row['commit']}:{row['representative_path']}",
                )

    def test_branch_sources_are_immutable_and_paths_resolve(self) -> None:
        rows = self.ledger["branch_sources"]
        ids = [row["id"] for row in rows]
        self.assertEqual(len(ids), len(set(ids)))
        allowed = set(self.ledger["frontier_schema"]["integration_statuses"])
        self.assertEqual(allowed, {"main_integrated", "open_pr", "branch_only", "hypothesis"})
        for row in rows:
            with self.subTest(source=row["id"]):
                self.assertIn(row["integration_status"], allowed)
                self.assertRegex(row["commit"], HEX40)
                self.assertNotIn("origin/", row["ref"])
                if row["integration_status"] == "open_pr":
                    self.assertIsInstance(row.get("pr"), int)
                for relative in row["paths"]:
                    assert_repo_relative(self, relative)
                    self.assertTrue(
                        git_object_exists(row["commit"], relative),
                        f"{row['commit']}:{relative}",
                    )

    def test_scientific_nodes_use_declared_dimensions_and_references(self) -> None:
        schema = self.ledger["frontier_schema"]
        nodes = self.ledger["scientific_nodes"]
        node_ids = [row["id"] for row in nodes]
        self.assertEqual(len(node_ids), len(set(node_ids)))
        self.assertGreaterEqual(len(nodes), 12)
        source_ids = {row["id"] for row in self.ledger["branch_sources"]}
        observer_ids = {row["id"] for row in self.ledger["observer_sectors"]}
        experiment_ids = {row["id"] for row in self.ledger["decision_experiments"]}
        statuses_seen = set()
        for row in nodes:
            with self.subTest(node=row["id"]):
                statuses_seen.add(row["integration_status"])
                self.assertIn(row["integration_status"], schema["integration_statuses"])
                self.assertIn(row["state"], schema["states"])
                self.assertTrue(set(row["sources"]).issubset(schema["sources"]))
                self.assertTrue(set(row["geometries"]).issubset(schema["geometries"]))
                self.assertTrue(set(row["acquisition"]).issubset(schema["acquisitions"]))
                self.assertTrue(set(row["source_refs"]).issubset(source_ids))
                self.assertTrue(set(row["observer_sectors"]).issubset(observer_ids))
                self.assertTrue(set(row["decision_experiments"]).issubset(experiment_ids))
                self.assertIn(row["information_gain"], schema["information_gain"])
                self.assertIn(row["compute_cost"], schema["compute_cost"])
                for key in ("known_result", "boundary", "next_discriminator", "supports", "tensions", "excludes_only"):
                    self.assertTrue(row[key], key)
        self.assertEqual(statuses_seen, set(schema["integration_statuses"]))

    def test_decision_experiments_form_the_declared_attention_order(self) -> None:
        experiments = self.ledger["decision_experiments"]
        ids = [row["id"] for row in experiments]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual([row["attention_rank"] for row in experiments], list(range(1, len(experiments) + 1)))
        self.assertEqual(self.ledger["attention_order"]["decision_experiment_refs"], ids)
        self.assertTrue(self.ledger["attention_order"]["parallel"])
        nodes = {row["id"] for row in self.ledger["scientific_nodes"]}
        sources = {row["id"] for row in self.ledger["branch_sources"]}
        observers = {row["id"] for row in self.ledger["observer_sectors"]}
        dependencies = {row["id"] for row in self.ledger["dependency_groups"]}
        schema = self.ledger["frontier_schema"]
        for row in experiments:
            with self.subTest(experiment=row["id"]):
                self.assertTrue(set(row["scientific_node_refs"]).issubset(nodes))
                self.assertTrue(set(row["source_refs"]).issubset(sources))
                self.assertTrue(set(row["observer_sector_refs"]).issubset(observers))
                self.assertTrue(set(row["dependency_group_refs"]).issubset(dependencies))
                self.assertIn(row["acquisition"], schema["acquisitions"])
                self.assertIn(row["expected_model_space_reduction"], schema["information_gain"])
                self.assertIn(row["compute_cost"], schema["compute_cost"])
                self.assertTrue(row["decision_output"])

    def test_two_activation_dependency_groups_match_manifest(self) -> None:
        manifest = load_yaml(TWO_ACTIVATION_MANIFEST)
        groups = {row["id"]: row for row in self.ledger["dependency_groups"]}
        self.assertEqual(len(groups), len(self.ledger["dependency_groups"]))
        self.assertEqual({run["dependency_group"] for run in manifest["runs"]}, set(groups))
        manifest_sizes = sorted(int(run["N"]) for run in manifest["runs"])
        ledger_sizes = sorted(size for row in groups.values() for size in row["member_sizes"])
        self.assertEqual(manifest_sizes, ledger_sizes)
        metadata_in_manifest = {run["metadata"] for run in manifest["runs"]}
        metadata_in_ledger = {path for row in groups.values() for path in row["metadata_paths"]}
        self.assertEqual(metadata_in_manifest, metadata_in_ledger)
        for relative in metadata_in_ledger:
            assert_repo_relative(self, relative)
            self.assertTrue((ROOT / relative).is_file(), relative)

    def test_registered_two_activation_reanalysis_contract(self) -> None:
        rows = {row["id"]: row for row in self.registry["registered_reanalyses"]}
        row = rows["TWO_ACTIVATION_H4"]
        self.assertEqual(row["schema"], "matching-one.two-activation-h4.v1")
        for key in ("manifest", "entrypoint"):
            assert_repo_relative(self, row[key])
            self.assertTrue((ROOT / row[key]).is_file(), row[key])
        for relative in row["outputs"]:
            assert_repo_relative(self, relative)
            self.assertTrue((ROOT / relative).is_file(), relative)
        result = load_yaml(ROOT / row["outputs"][0])
        self.assertEqual(result["schema"], row["schema"])

    def test_machine_readable_frontier_has_no_permission_states(self) -> None:
        rendered = yaml.safe_dump(self.ledger, sort_keys=True).lower()
        for forbidden in ("gated", "locked"):
            with self.subTest(forbidden=forbidden):
                self.assertIsNone(re.search(rf"\b{forbidden}\b", rendered))


if __name__ == "__main__":
    unittest.main()
