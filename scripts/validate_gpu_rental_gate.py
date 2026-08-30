#!/usr/bin/env python3
"""Validate the exact issue #30 CPU baseline and route-1 GPU rental gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "gpu_rental_gate_manifest.json"
EXPECTED_SCHEMA = "matching-one/gpu-rental-gate/v1"
FORBIDDEN_KEYS = {
    "estimated_gpu_throughput",
    "estimated_speedup",
    "gpu_elapsed_seconds",
    "gpu_paired_permutations",
    "kernel_only_speedup",
    "rental_authorized",
    "production_authorized",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _walk_forbidden(value: Any, path: str = "$") -> None:
    if isinstance(value, dict):
        bad = sorted(set(value) & FORBIDDEN_KEYS)
        _require(not bad, f"{path} contains forbidden GPU-result keys: {','.join(bad)}")
        for key, child in value.items():
            _walk_forbidden(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _walk_forbidden(child, f"{path}[{index}]")


def _fraction(value: Any, label: str) -> Fraction:
    _require(isinstance(value, str), f"{label} must be an exact fraction string")
    _require(re.fullmatch(r"-?\d+(?:/[1-9]\d*)?", value) is not None, f"{label} is not canonical")
    result = Fraction(value)
    _require(str(result) == value, f"{label} is not reduced")
    return result


def _source_bytes(contract: Mapping[str, Any], root: Path, source: bytes | None) -> bytes:
    if source is not None:
        return source
    path = contract.get("baseline_source", {}).get("path")
    _require(isinstance(path, str), "baseline source path is missing")
    candidate = root / path
    _require(candidate.is_file(), "baseline source is unavailable")
    return candidate.read_bytes()


def validate_contract(
    contract: Mapping[str, Any], *, root: Path = ROOT, source: bytes | None = None
) -> dict[str, Any]:
    """Validate the checked-in no-GPU-result state and exact CPU arithmetic."""

    _walk_forbidden(contract)
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 30, "issue must be 30")
    _require(contract.get("status") == "not_authorized", "checked-in status must be not_authorized")

    source_contract = contract.get("baseline_source", {})
    _require(
        source_contract.get("path") == "experiments/server_compute_queue_overrides_20260828_v3.yaml",
        "baseline source path drifted",
    )
    data = _source_bytes(contract, root, source)
    digest = hashlib.sha256(data).hexdigest()
    _require(source_contract.get("sha256") == digest, "baseline source SHA-256 mismatch")
    _require(
        b"10M paired permutations x five sizes in 65.2 sec on 16-core ARM" in data,
        "CPU baseline marker is absent from source",
    )
    _require(b"status: not_required_yet" in data, "source no longer records GPU as not required")

    cpu = contract.get("cpu_baseline", {})
    _require(cpu.get("pair_count") == 5, "CPU pair count drifted")
    _require(cpu.get("paired_permutations_per_pair") == 10_000_000, "CPU per-pair count drifted")
    total = cpu.get("paired_permutations_total")
    _require(total == 50_000_000, "CPU total count drifted")
    _require(total == cpu["pair_count"] * cpu["paired_permutations_per_pair"], "CPU count product is inconsistent")
    elapsed = _fraction(cpu.get("elapsed_seconds"), "CPU elapsed seconds")
    _require(elapsed == Fraction(326, 5), "CPU elapsed time drifted")
    _require(cpu.get("environment") == "16-core ARM", "CPU environment drifted")
    _require(cpu.get("exact_output_oracle_passed") is True, "CPU output oracle is not recorded as passing")

    cpu_throughput = Fraction(total, 1) / elapsed
    gate = contract.get("end_to_end_gate", {})
    minimum_speedup = _fraction(gate.get("minimum_speedup"), "minimum speedup")
    _require(minimum_speedup == 5, "minimum speedup must remain 5")
    _require(gate.get("timing_scope") == "end_to_end", "kernel-only timing is forbidden")
    _require(gate.get("output_contract_equality_required") is True, "output equality gate is disabled")
    _require(gate.get("deterministic_regression_required") is True, "determinism gate is disabled")
    minimum_gpu_throughput = cpu_throughput * minimum_speedup

    derived = contract.get("derived_exact", {})
    _require(
        _fraction(derived.get("cpu_paired_permutations_per_second"), "derived CPU throughput")
        == cpu_throughput,
        "derived CPU throughput drifted",
    )
    _require(
        _fraction(derived.get("minimum_gpu_paired_permutations_per_second"), "derived GPU threshold")
        == minimum_gpu_throughput,
        "derived GPU threshold drifted",
    )
    _require(contract.get("gpu_measurement") == {"status": "absent"}, "a GPU measurement is checked in")
    _require(
        contract.get("decision")
        == {
            "rental": "not_authorized",
            "reason": "No qualifying end-to-end GPU measurement is checked in.",
        },
        "checked-in rental decision drifted",
    )

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_cpu_baseline_only",
        "cpu_paired_permutations": total,
        "cpu_elapsed_seconds": str(elapsed),
        "cpu_paired_permutations_per_second": str(cpu_throughput),
        "minimum_speedup": str(minimum_speedup),
        "minimum_gpu_paired_permutations_per_second": str(minimum_gpu_throughput),
        "gpu_measurement_present": False,
        "rental": "not_authorized",
        "parent_issue": "remain open",
    }


def evaluate_measurement(
    contract: Mapping[str, Any],
    *,
    gpu_paired_permutations: int,
    gpu_elapsed_seconds: Fraction,
    end_to_end_timing: bool,
    output_contract_equal: bool,
    deterministic_regression_passed: bool,
) -> dict[str, Any]:
    """Evaluate, but do not persist, a future route-1 GPU measurement."""

    _require(
        isinstance(gpu_paired_permutations, int)
        and not isinstance(gpu_paired_permutations, bool)
        and gpu_paired_permutations > 0,
        "GPU permutation count must be a positive integer",
    )
    _require(isinstance(gpu_elapsed_seconds, Fraction) and gpu_elapsed_seconds > 0, "GPU elapsed time must be a positive Fraction")
    flags = (end_to_end_timing, output_contract_equal, deterministic_regression_passed)
    _require(all(isinstance(flag, bool) for flag in flags), "measurement gates must be booleans")

    cpu = contract.get("cpu_baseline", {})
    _require(cpu.get("paired_permutations_total") == 50_000_000, "CPU total count drifted")
    cpu_elapsed = _fraction(cpu.get("elapsed_seconds"), "CPU elapsed seconds")
    _require(cpu_elapsed == Fraction(326, 5), "CPU elapsed time drifted")
    cpu_throughput = Fraction(cpu["paired_permutations_total"], 1) / cpu_elapsed
    gate = contract.get("end_to_end_gate", {})
    minimum_speedup = _fraction(gate.get("minimum_speedup"), "minimum speedup")
    _require(minimum_speedup == 5, "minimum speedup must remain 5")
    _require(gate.get("timing_scope") == "end_to_end", "kernel-only timing is forbidden")
    _require(gate.get("output_contract_equality_required") is True, "output equality gate is disabled")
    _require(gate.get("deterministic_regression_required") is True, "determinism gate is disabled")
    gpu_throughput = Fraction(gpu_paired_permutations, 1) / gpu_elapsed_seconds
    speedup = gpu_throughput / cpu_throughput
    authorized = speedup >= minimum_speedup and all(flags)
    return {
        "gpu_paired_permutations_per_second": str(gpu_throughput),
        "speedup": str(speedup),
        "throughput_gate_passed": speedup >= minimum_speedup,
        "end_to_end_timing": end_to_end_timing,
        "output_contract_equal": output_contract_equal,
        "deterministic_regression_passed": deterministic_regression_passed,
        "rental": "authorized_by_route_1" if authorized else "not_authorized",
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
