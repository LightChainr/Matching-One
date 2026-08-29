from __future__ import annotations

from pathlib import Path
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "analysis" / "artifact_registry.yaml"


class ArtifactRegistryConsistencyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.registry = yaml.safe_load(REGISTRY.read_text(encoding="utf-8"))

    def _frontier_prs(self) -> set[int]:
        return {int(row["pr"]) for row in self.registry["frontier_open_as_of_2026_08_29"]}

    def _canonical_prs(self) -> set[int]:
        return {int(row["pr"]) for row in self.registry["canonical_integration_history"]}

    def _superseded_prs(self) -> set[int]:
        return {int(row["pr"]) for row in self.registry["superseded_active_paths_closed"]}

    def test_frontier_is_not_also_canonical(self) -> None:
        overlap = self._frontier_prs() & self._canonical_prs()
        self.assertEqual(overlap, set())

    def test_frontier_is_not_also_superseded(self) -> None:
        overlap = self._frontier_prs() & self._superseded_prs()
        self.assertEqual(overlap, set())

    def test_closed_reintegrated_prs_are_not_frontier(self) -> None:
        # Source PRs whose content landed through later clean integrations.
        reintegrated = {149, 150, 178, 179, 181, 183, 184}
        frontier = self._frontier_prs()
        self.assertEqual(reintegrated & frontier, set())
        canonical = self._canonical_prs()
        for pr in (178, 179, 181, 183, 184):
            self.assertIn(pr, canonical)
        superseded = self._superseded_prs()
        self.assertIn(149, superseded)
        self.assertIn(150, superseded)

    def test_every_frontier_row_has_a_classification(self) -> None:
        allowed = {
            "unique_frontier",
            "partially_integrated_with_remaining_deliverable",
            "superseded_close_recommended",
            "historical_negative_result_to_extract",
        }
        for row in self.registry["frontier_open_as_of_2026_08_29"]:
            with self.subTest(pr=row["pr"]):
                self.assertIn(row["classification"], allowed)

    def test_canonical_documents_exist(self) -> None:
        for row in self.registry["canonical_documents"]:
            path = ROOT / row["path"]
            with self.subTest(path=row["path"]):
                self.assertTrue(path.is_file(), row["path"])


if __name__ == "__main__":
    unittest.main()
