#!/usr/bin/env python3
"""Rebuild the #537 witness and score, then run the independent topology oracle."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = ROOT / "results" / "p537-one-defect-gate-20260901"
KERNEL = ROOT / "experiments" / "p537-landing-matrix-preflight-20260901" / "kernel.tsv"


def run(command: list[str], cwd: Path = ROOT) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(
        command,
        cwd=cwd,
        check=True,
        text=True,
        capture_output=True,
        env=environment,
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_manifest() -> int:
    checked = 0
    manifest = HERE / "SHA256SUMS"
    for line in manifest.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        path = ROOT / relative
        observed = sha256(path)
        if observed != expected:
            raise AssertionError(f"SHA256 mismatch for {relative}: {observed} != {expected}")
        checked += 1
    return checked


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="p537-one-defect-") as directory:
        scratch = Path(directory)
        binary = scratch / "one_defect_witness_exact"
        witness = scratch / "witness.json"
        result = scratch / "result.json"
        compiler = os.environ.get("CXX", "c++")
        run(
            [
                compiler,
                "-O3",
                "-std=c++17",
                str(HERE / "one_defect_witness_exact.cpp"),
                "-o",
                str(binary),
            ]
        )
        produced = run([str(binary), "5", "0", str(KERNEL), str(witness)])
        run_record = json.loads(produced.stdout)
        if run_record["status"] != "first_observable_diagonal_edge_stop":
            raise AssertionError(run_record)
        if run_record["backgrounds_scanned"] != 12568:
            raise AssertionError(run_record)
        if witness.read_bytes() != (RESULTS / "witness.json").read_bytes():
            raise AssertionError("C++ producer did not reproduce witness.json byte for byte")

        run([sys.executable, str(HERE / "score_witness.py"), "--output", str(result)])
        if result.read_bytes() != (RESULTS / "result.json").read_bytes():
            raise AssertionError("pooled-root scorer did not reproduce result.json byte for byte")

        oracle = run([sys.executable, str(ROOT / "tests" / "test_p537_one_defect_gate.py"), "-v"])
        manifest_files = verify_manifest()
        print(
            json.dumps(
                {
                    "status": "verified",
                    "producer": "byte_identical",
                    "scorer": "byte_identical",
                    "topology_tests": 3,
                    "manifest_files": manifest_files,
                    "compiler": compiler,
                    "test_stdout": oracle.stdout.strip(),
                    "test_stderr": oracle.stderr.strip(),
                },
                sort_keys=True,
            )
        )


if __name__ == "__main__":
    main()
