#!/usr/bin/env python3
"""Validate the desired local main-ruleset contract against the CI workflow."""

from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path
from typing import Any, Mapping

import yaml


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "main_ruleset_contract.json"
EXPECTED_SCHEMA = "matching-one/main-ruleset-contract/v1"
FORBIDDEN_HOSTING_CLAIMS = {
    "applied_at",
    "enabled",
    "protected",
    "ruleset_id",
    "verified_at",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _walk_forbidden(value: Any, path: str = "$") -> None:
    if isinstance(value, dict):
        bad = sorted(set(value) & FORBIDDEN_HOSTING_CLAIMS)
        _require(not bad, f"{path} contains forbidden hosting-state claims: {','.join(bad)}")
        for key, child in value.items():
            _walk_forbidden(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _walk_forbidden(child, f"{path}[{index}]")


def load_workflow(path: Path) -> Mapping[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    _require(isinstance(payload, dict), "workflow must be a mapping")
    return payload


def expanded_check_names(contract_jobs: Mapping[str, Any]) -> list[str]:
    template = contract_jobs.get("python_name_template")
    versions = contract_jobs.get("python_versions")
    cxx_name = contract_jobs.get("cxx_name")
    _require(isinstance(template, str), "python job name template must be a string")
    _require(isinstance(versions, list) and all(isinstance(item, str) for item in versions), "python versions must be strings")
    _require(isinstance(cxx_name, str), "C++ job name must be a string")
    token = "${{ matrix.python-version }}"
    _require(token in template, "python job name template lost its matrix token")
    return [template.replace(token, version) for version in versions] + [cxx_name]


def validate_contract(
    contract: Mapping[str, Any], workflow: Mapping[str, Any]
) -> dict[str, Any]:
    _walk_forbidden(contract)
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 52, "issue must be 52")
    _require(contract.get("status") == "desired_local_contract_only", "status drifted")
    _require(contract.get("target_branch") == "main", "target branch must remain main")
    _require(contract.get("workflow_name") == "ci", "workflow name contract drifted")
    _require(workflow.get("name") == contract.get("workflow_name"), "workflow name drifted")

    jobs = workflow.get("jobs")
    _require(isinstance(jobs, dict), "workflow jobs must be a mapping")
    declared_jobs = contract.get("workflow_jobs")
    _require(isinstance(declared_jobs, dict), "workflow job contract must be an object")
    python_id = declared_jobs.get("python_job_id")
    cxx_id = declared_jobs.get("cxx_job_id")
    _require(python_id in jobs, "python workflow job is missing")
    _require(cxx_id in jobs, "C++ workflow job is missing")

    python_job = jobs[python_id]
    cxx_job = jobs[cxx_id]
    _require(isinstance(python_job, dict) and isinstance(cxx_job, dict), "workflow jobs must be objects")
    _require(
        python_job.get("name") == declared_jobs.get("python_name_template"),
        "python workflow job name drifted",
    )
    _require(cxx_job.get("name") == declared_jobs.get("cxx_name"), "C++ workflow job name drifted")
    matrix = python_job.get("strategy", {}).get("matrix", {}).get("python-version")
    _require(matrix == declared_jobs.get("python_versions"), "Python version matrix drifted")

    expanded = expanded_check_names(declared_jobs)
    _require(contract.get("required_status_checks") == expanded, "required status checks drifted")

    policy = contract.get("desired_policy")
    _require(isinstance(policy, dict), "desired policy must be an object")
    expected_policy = {
        "require_pull_request": True,
        "required_approving_review_count": 0,
        "block_force_pushes": True,
        "block_deletions": True,
        "require_linear_history": False,
        "require_external_approvals": False,
    }
    _require(policy == expected_policy, "desired policy drifted")
    _require(
        contract.get("hosting_side_state") == "not_checked_or_modified",
        "hosting-side boundary drifted",
    )

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_local_ruleset_contract",
        "target_branch": "main",
        "required_status_checks": expanded,
        "desired_policy": deepcopy(expected_policy),
        "workflow_contract_matches": True,
        "ready_for_manual_application": True,
        "hosting_side_state_checked": False,
        "hosting_side_state_modified": False,
        "parent_issue": "remain open",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    workflow_path = ROOT / contract.get("workflow_path", "")
    workflow = load_workflow(workflow_path)
    print(json.dumps(validate_contract(contract, workflow), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
