#!/usr/bin/env python3
"""Rebuild the exact N25 profiles in scratch space and compare every artifact."""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def run(command: list[str], *, cwd: Path, env: dict[str, str] | None = None) -> None:
    subprocess.run(command, cwd=cwd, env=env, check=True, stdout=subprocess.PIPE, text=True)


def main() -> None:
    compiler = os.environ.get("CXX") or shutil.which("clang++") or shutil.which("g++")
    if not compiler:
        raise SystemExit("no C++ compiler found")
    inputs = [
        "regular_pair_joint_u_exact.cpp",
        "score.py",
        "validate_landing_contract.py",
        "kernel.tsv",
        "baseline-root.json",
        "baseline-axis.csv",
        "baseline-tilted.csv",
        "original-axis.csv",
        "original-tilted.csv",
    ]
    generated_profiles = [
        f"{geometry}-{mode}.csv"
        for geometry in ("axis", "tilted")
        for mode in ("clean_same", "clean_reversed", "clean_total", "all")
    ]
    with tempfile.TemporaryDirectory(prefix="matching-one-p537-") as raw:
        scratch = Path(raw)
        for name in inputs:
            shutil.copy2(HERE / name, scratch / name)
        binary = scratch / "landing_enum"
        run([compiler, "-O3", "-std=c++17", "regular_pair_joint_u_exact.cpp", "-o", str(binary)], cwd=scratch)
        for geometry, a, b in (("axis", "5", "0"), ("tilted", "4", "3")):
            for mode in ("clean_same", "clean_reversed", "clean_total", "all"):
                run(
                    [str(binary), a, b, "kernel.tsv", f"{geometry}-{mode}.csv", mode],
                    cwd=scratch,
                )
        env = os.environ.copy()
        env["MATCHING_ONE_THERMAL_GATE"] = str(
            ROOT / "experiments/p337-thermal-gate-audit-20260901/thermal_gate.py"
        )
        run(["python3", "validate_landing_contract.py"], cwd=scratch, env=env)
        run(["python3", "score.py"], cwd=scratch)
        for name in generated_profiles + ["landing_contract_validation.json", "result.json"]:
            if (scratch / name).read_bytes() != (HERE / name).read_bytes():
                raise AssertionError(f"scratch rebuild differs: {name}")

        n13_output = scratch / "n13-transition-result.json"
        run(
            ["python3", str(HERE / "n13_transition_channels.py"), "--out", str(n13_output)],
            cwd=ROOT,
        )
        if n13_output.read_bytes() != (HERE / "n13-transition-result.json").read_bytes():
            raise AssertionError("scratch rebuild differs: n13-transition-result.json")

    print(json.dumps({"status": "verified", "N25_exact_traversals": 8, "random_samples": 0}))


if __name__ == "__main__":
    main()
