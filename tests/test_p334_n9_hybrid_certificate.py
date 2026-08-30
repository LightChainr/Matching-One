import json
from pathlib import Path

import pytest

from p334_n9_hybrid_certificate import (
    DIRECT_INDICES,
    EXPECTED_FAILURES,
    HEAVY_INDICES,
    SOURCE_COMMIT,
    validate_direct_directory,
)


ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = ROOT / "results" / "p334-n9-hybrid-certificate"


def test_frozen_direct_artifacts_are_complete_and_byte_identified():
    reports = validate_direct_directory(RESULT_DIR / "direct")
    assert [report["row_index"] for report in reports] == list(DIRECT_INDICES)
    assert all(report["source_sha256"] for report in reports)


def test_checked_in_hybrid_result_closes_the_shell():
    result = json.loads((RESULT_DIR / "latest.json").read_text())
    assert result["source_commit"] == SOURCE_COMMIT
    assert result["coverage"] == {
        "complete": True,
        "expected": list(range(28)),
        "observed": list(range(28)),
    }
    assert result["summary"]["saturated_rows"] == 22
    assert result["summary"]["failed_rows"] == 6
    assert result["summary"]["failed_row_indices"] == list(EXPECTED_FAILURES)
    assert result["summary"]["failure_Hall_deficiency"] == 2160
    for index in HEAVY_INDICES:
        report = result["rows"][index]
        assert report["proof_path"] == "coarse_twin_capacitated_hall"
        assert report["flow"] == {
            "Hall_deficiency": 0,
            "maximum_flow": 45360,
            "saturates": True,
            "source_tokens": 45360,
        }


def test_hash_gate_rejects_a_mutated_direct_artifact(tmp_path):
    source = RESULT_DIR / "direct"
    for path in source.glob("*.json"):
        (tmp_path / path.name).write_bytes(path.read_bytes())
    path = tmp_path / "row0.json"
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(AssertionError, match="SHA256 mismatch"):
        validate_direct_directory(tmp_path)
