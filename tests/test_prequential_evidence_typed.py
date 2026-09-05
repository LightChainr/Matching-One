
from __future__ import annotations
import copy
import json
from pathlib import Path
import sys
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import score_prequential_evidence_typed as typed  # noqa: E402


class PrequentialEvidenceTypedTests(unittest.TestCase):
    @staticmethod
    def manifest() -> dict:
        return json.loads(
            (ROOT / "analysis/evidence_ledger_manifest.yaml").read_text(encoding="utf-8")
        )

    def test_typed_replay_preserves_every_frozen_ledger_field(self) -> None:
        expected = json.loads(
            (ROOT / "results/evidence-ledger/latest.json").read_text(encoding="utf-8")
        )
        scorer = mock.Mock(return_value=copy.deepcopy(expected))
        result = typed.score_manifest_typed(ROOT, self.manifest(), scorer=scorer)
        semantics = result.pop("observable_semantics")
        self.assertEqual(result, expected)
        self.assertEqual(len(semantics["blocks"]), 9)
        scorer.assert_called_once()

    def test_truthy_legacy_exact_map_is_rejected_before_score(self) -> None:
        manifest = self.manifest()
        manifest["blocks"][0]["channel"]["exact_map"] = "trust me"
        scorer = mock.Mock()
        with self.assertRaisesRegex(ValueError, "frozen canonical manifest"):
            typed.score_manifest_typed(ROOT, manifest, scorer=scorer)
        scorer.assert_not_called()


if __name__ == "__main__":
    unittest.main()
