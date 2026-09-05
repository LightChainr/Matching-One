
from __future__ import annotations
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import observable_descriptor_map_audit as audit  # noqa: E402


class ObservableDescriptorMapAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = audit.build_artifact()

    def test_descriptor_enumeration_is_complete_and_unique(self) -> None:
        descriptors = audit.enumerate_descriptors()
        self.assertEqual(len(descriptors), 200)
        self.assertEqual(len({audit.descriptor_id(value) for value in descriptors}), 200)

    def test_blocked_pairs_fail_for_only_declared_reasons(self) -> None:
        self.assertEqual(self.artifact["blocked_reasons"], {
            "cannot map a scalar value to an orientation contrast": 12800,
            "no exact topology map": 26248,
        })

    def test_all_composable_paths_equal_direct_map(self) -> None:
        checks = self.artifact["composition_checks"]
        self.assertEqual(checks["checked_paths"], 5720)
        self.assertEqual(checks["failures"], [])

    def test_connected_component_partition_is_frozen(self) -> None:
        self.assertEqual(self.artifact["component_size_histogram"], {"1": 24, "2": 8, "4": 24, "8": 8})
        components = self.artifact["connected_components"]
        self.assertEqual(len(components), 64)
        self.assertEqual(sum(item["size"] for item in components), 200)

    def test_checked_in_artifacts_reproduce(self) -> None:
        checked_json = json.loads((ROOT / "results/observable-descriptor-map-audit/latest.json").read_text())
        checked_md = (ROOT / "results/observable-descriptor-map-audit/latest.md").read_text()
        self.assertEqual(checked_json, self.artifact)
        self.assertEqual(checked_md, audit.render_markdown(self.artifact))


if __name__ == "__main__":
    unittest.main()
