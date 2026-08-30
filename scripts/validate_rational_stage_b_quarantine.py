#!/usr/bin/env python3
"""Validate the issue #27 Stage-B pre-reveal quarantine contract."""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import io
import json
from pathlib import Path
import re
from typing import Any, Callable, Mapping


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "rational_stage_b_quarantine_manifest.json"

EXPECTED_SCHEMA = "matching-one/rational-stage-b-quarantine/v1"
EXPECTED_FAMILIES = [
    "poly_0",
    "poly_1",
    "poly_2",
    "poly_3",
    "poly_4",
    "pade_1_1",
    "pade_2_1",
    "pade_1_2",
    "pade_2_2",
]
EXPECTED_METRICS = [
    "signed_error_per_width",
    "root_mean_square_error",
    "maximum_absolute_error",
    "signed_residual_trend",
]
FORBIDDEN_KEYS = {
    "actual",
    "actuals",
    "observed_targets",
    "revealed_targets",
    "target_values",
    "target_scores",
    "selected_after_reveal",
    "post_cutoff_fit",
    "post_cutoff_values",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _walk_forbidden(value: Any, path: str = "$") -> None:
    if isinstance(value, dict):
        bad = sorted(set(value) & FORBIDDEN_KEYS)
        _require(not bad, f"{path} contains forbidden target-data keys: {','.join(bad)}")
        for key, child in value.items():
            _walk_forbidden(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _walk_forbidden(child, f"{path}[{index}]")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _family_names(source: bytes) -> list[str]:
    tree = ast.parse(source.decode("utf-8"))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(isinstance(target, ast.Name) and target.id == "FAMILIES" for target in node.targets):
            continue
        _require(isinstance(node.value, (ast.List, ast.Tuple)), "FAMILIES must be a literal sequence")
        names: list[str] = []
        for element in node.value.elts:
            _require(isinstance(element, ast.Call), "FAMILIES entries must be constructor calls")
            _require(bool(element.args), "FAMILIES entry is missing its name")
            name = ast.literal_eval(element.args[0])
            _require(isinstance(name, str), "FAMILIES names must be strings")
            names.append(name)
        return names
    raise ValueError("audit source does not define FAMILIES")


def _training_widths(data: bytes) -> list[int]:
    reader = csv.DictReader(io.StringIO(data.decode("utf-8")))
    _require(
        reader.fieldnames == ["n", "value", "source", "method"],
        "training CSV header drifted",
    )
    widths: list[int] = []
    for row in reader:
        _require(all(row.get(field, "") for field in reader.fieldnames), "training CSV contains an empty field")
        try:
            widths.append(int(row["n"]))
        except ValueError as exc:
            raise ValueError("training CSV contains a non-integer width") from exc
    return widths


def _frozen_prediction_widths(data: bytes) -> list[int]:
    text = data.decode("utf-8")
    _require(re.search(r"^status:\s*preregistered\s*$", text, re.MULTILINE) is not None, "prediction is not preregistered")
    return [int(value) for value in re.findall(r"^\s*-\s+n:\s*(\d+)\s*$", text, re.MULTILINE)]


def validate_contract(
    contract: Mapping[str, Any],
    *,
    root: Path = ROOT,
    artifacts: Mapping[str, bytes] | None = None,
    path_exists: Callable[[str], bool] | None = None,
) -> dict[str, Any]:
    """Fail closed unless every pre-reveal boundary is still intact."""

    _walk_forbidden(contract)
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 27, "issue must be 27")
    _require(contract.get("status") == "pre_reveal_quarantine", "contract is not pre-reveal")

    cutoff = contract.get("knowledge_cutoff")
    targets = contract.get("target_widths")
    _require(cutoff == 21, "knowledge cutoff must remain 21")
    _require(targets == [22, 23, 24], "target widths must remain 22,23,24")
    _require(all(width > cutoff for width in targets), "target widths must be beyond the cutoff")

    def read_artifact(path: str) -> bytes:
        if artifacts is not None:
            _require(path in artifacts, f"artifact is unavailable: {path}")
            return artifacts[path]
        candidate = root / path
        _require(candidate.is_file(), f"artifact is unavailable: {path}")
        return candidate.read_bytes()

    if path_exists is None:
        path_exists = lambda path: (root / path).exists()

    training = contract.get("training_source", {})
    training_path = training.get("path")
    _require(isinstance(training_path, str), "training source path is missing")
    training_bytes = read_artifact(training_path)
    training_digest = _sha256(training_bytes)
    _require(training.get("sha256") == training_digest, "training source SHA-256 mismatch")
    _require(training.get("required_widths") == {"first": 1, "last": 21}, "training range drifted")
    widths = _training_widths(training_bytes)
    _require(widths == list(range(1, cutoff + 1)), "training widths must be exactly 1..21")
    _require(not set(widths).intersection(targets), "target width leaked into training data")

    audit = contract.get("audit_source", {})
    audit_path = audit.get("path")
    _require(isinstance(audit_path, str), "audit source path is missing")
    audit_bytes = read_artifact(audit_path)
    audit_digest = _sha256(audit_bytes)
    _require(audit.get("sha256") == audit_digest, "audit source SHA-256 mismatch")
    families = _family_names(audit_bytes)
    _require(contract.get("audit_families") == EXPECTED_FAMILIES, "contract family set drifted")
    _require(families == EXPECTED_FAMILIES, "audit source family set drifted")

    prediction = contract.get("frozen_prediction", {})
    prediction_path = prediction.get("path")
    _require(isinstance(prediction_path, str), "frozen prediction path is missing")
    prediction_bytes = read_artifact(prediction_path)
    prediction_digest = _sha256(prediction_bytes)
    _require(prediction.get("sha256") == prediction_digest, "frozen prediction SHA-256 mismatch")
    _require(_frozen_prediction_widths(prediction_bytes) == targets, "frozen prediction widths drifted")
    input_digest_line = f'input_sha256: "{training_digest}"'.encode("ascii")
    _require(input_digest_line in prediction_bytes, "prediction training digest does not match")

    _require(contract.get("frozen_metrics") == EXPECTED_METRICS, "frozen metrics drifted")
    permissions = contract.get("permissions")
    _require(
        permissions
        == {
            "model_selection_permitted": True,
            "target_scoring_permitted": False,
            "post_cutoff_refit_permitted": False,
        },
        "pre-reveal permissions drifted",
    )
    target_artifact = contract.get("target_artifact", {})
    target_path = target_artifact.get("expected_path")
    _require(isinstance(target_path, str), "target artifact path is missing")
    _require(target_artifact.get("must_be_absent") is True, "target absence guard is disabled")
    _require(not path_exists(target_path), "target artifact exists before the quarantine is lifted")

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_pre_reveal_quarantine",
        "knowledge_cutoff": cutoff,
        "target_widths": targets,
        "training_width_count": len(widths),
        "audit_family_count": len(families),
        "training_sha256": training_digest,
        "audit_sha256": audit_digest,
        "frozen_prediction_sha256": prediction_digest,
        "target_data_present": False,
        "safe_to_select_models": True,
        "safe_to_score_targets": False,
        "parent_issue": "remain open",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    print(json.dumps(validate_contract(contract), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
