
from __future__ import annotations
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from validate_rational_stage_b_quarantine import validate_contract  # noqa: E402


class RationalStageBQuarantineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "rational_stage_b_quarantine_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        paths = [
            cls.contract["training_source"]["path"],
            cls.contract["audit_source"]["path"],
            cls.contract["frozen_prediction"]["path"],
        ]
        cls.artifacts = {path: (ROOT / path).read_bytes() for path in paths}

    def validate(self, contract=None, artifacts=None, existing=()):
        return validate_contract(
            self.contract if contract is None else contract,
            artifacts=self.artifacts if artifacts is None else artifacts,
            path_exists=lambda path: path in set(existing),
        )

    def test_checked_in_contract_freezes_pre_reveal_inputs(self) -> None:
        result = self.validate()
        self.assertEqual(result["knowledge_cutoff"], 21)
        self.assertEqual(result["target_widths"], [22, 23, 24])
        self.assertEqual(result["training_width_count"], 21)
        self.assertEqual(result["audit_family_count"], 9)
        self.assertTrue(result["safe_to_select_models"])
        self.assertFalse(result["safe_to_score_targets"])
        self.assertFalse(result["target_data_present"])


if __name__ == "__main__":
    unittest.main()
