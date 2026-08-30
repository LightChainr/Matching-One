from __future__ import annotations

from copy import deepcopy
import hashlib
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

    def test_target_width_in_training_fails_even_with_updated_digest(self) -> None:
        contract = deepcopy(self.contract)
        artifacts = dict(self.artifacts)
        path = contract["training_source"]["path"]
        artifacts[path] += b"22,0.6,forbidden,forbidden\n"
        contract["training_source"]["sha256"] = hashlib.sha256(artifacts[path]).hexdigest()
        with self.assertRaisesRegex(ValueError, "exactly 1..21"):
            self.validate(contract, artifacts)

    def test_source_digest_and_family_drift_fail_closed(self) -> None:
        contract = deepcopy(self.contract)
        contract["audit_source"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "audit source SHA-256 mismatch"):
            self.validate(contract)

        contract = deepcopy(self.contract)
        artifacts = dict(self.artifacts)
        path = contract["audit_source"]["path"]
        artifacts[path] = artifacts[path].replace(b'Family("poly_4"', b'Family("poly_5"', 1)
        contract["audit_source"]["sha256"] = hashlib.sha256(artifacts[path]).hexdigest()
        with self.assertRaisesRegex(ValueError, "audit source family set drifted"):
            self.validate(contract, artifacts)

    def test_target_artifact_presence_fails_closed(self) -> None:
        target = self.contract["target_artifact"]["expected_path"]
        with self.assertRaisesRegex(ValueError, "target artifact exists"):
            self.validate(existing={target})

    def test_target_values_and_prediction_drift_are_rejected(self) -> None:
        contract = deepcopy(self.contract)
        contract["target_values"] = {"22": "0.6"}
        with self.assertRaisesRegex(ValueError, "forbidden target-data keys"):
            self.validate(contract)

        contract = deepcopy(self.contract)
        contract["frozen_prediction"]["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "frozen prediction SHA-256 mismatch"):
            self.validate(contract)


if __name__ == "__main__":
    unittest.main()
