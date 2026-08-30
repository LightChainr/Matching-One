from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from validate_main_ruleset_contract import (  # noqa: E402
    expanded_check_names,
    load_workflow,
    validate_contract,
)


class MainRulesetContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "main_ruleset_contract.json").read_text(encoding="utf-8")
        )
        cls.workflow = load_workflow(ROOT / ".github" / "workflows" / "ci.yml")

    def test_checked_in_contract_matches_workflow(self) -> None:
        result = validate_contract(self.contract, self.workflow)
        self.assertEqual(
            result["required_status_checks"],
            ["Python 3.9", "Python 3.11", "Python 3.13", "C++17 build and self-tests"],
        )
        self.assertTrue(result["workflow_contract_matches"])
        self.assertTrue(result["ready_for_manual_application"])
        self.assertFalse(result["hosting_side_state_checked"])
        self.assertFalse(result["hosting_side_state_modified"])

    def test_matrix_expansion_is_exact(self) -> None:
        self.assertEqual(
            expanded_check_names(self.contract["workflow_jobs"]),
            self.contract["required_status_checks"],
        )

    def test_workflow_name_and_matrix_drift_fail_closed(self) -> None:
        changed = deepcopy(self.workflow)
        changed["name"] = "other"
        with self.assertRaisesRegex(ValueError, "workflow name drifted"):
            validate_contract(self.contract, changed)
        changed = deepcopy(self.workflow)
        changed["jobs"]["python"]["strategy"]["matrix"]["python-version"] = ["3.11"]
        with self.assertRaisesRegex(ValueError, "Python version matrix drifted"):
            validate_contract(self.contract, changed)

    def test_missing_or_renamed_job_fails_closed(self) -> None:
        changed = deepcopy(self.workflow)
        del changed["jobs"]["cxx"]
        with self.assertRaisesRegex(ValueError, r"C\+\+ workflow job is missing"):
            validate_contract(self.contract, changed)
        changed = deepcopy(self.workflow)
        changed["jobs"]["python"]["name"] = "Python"
        with self.assertRaisesRegex(ValueError, "python workflow job name drifted"):
            validate_contract(self.contract, changed)

    def test_policy_drift_is_rejected(self) -> None:
        for field, value in (
            ("block_force_pushes", False),
            ("block_deletions", False),
            ("require_linear_history", True),
            ("required_approving_review_count", 1),
        ):
            changed = deepcopy(self.contract)
            changed["desired_policy"][field] = value
            with self.subTest(field=field):
                with self.assertRaisesRegex(ValueError, "desired policy drifted"):
                    validate_contract(changed, self.workflow)

    def test_hosting_state_claims_are_forbidden(self) -> None:
        changed = deepcopy(self.contract)
        changed["ruleset_id"] = 123
        with self.assertRaisesRegex(ValueError, "forbidden hosting-state claims"):
            validate_contract(changed, self.workflow)


if __name__ == "__main__":
    unittest.main()
