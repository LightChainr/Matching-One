#!/usr/bin/env python3
"""Compile once and enumerate the frozen Gaussian pair once, without scoring."""
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]
CPP = ROOT / "scripts/p337_closed_source_finite_exact.cpp"
CONTRACT = ROOT / "analysis/p337_closed_source_finite_coupling_contract.json"
OUT = ROOT / "results/p337-closed-source-finite-coupling"


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def main():
    contract = json.loads(CONTRACT.read_text())
    jobs = [("axis", 5, 0), ("tilted", 4, 3)]
    assert contract["geometries"] == [[a, b] for _, a, b in jobs]
    assert contract["workers"] == 2 and contract["N"] == 25
    OUT.mkdir(parents=True, exist_ok=True)
    for name in ("axis.csv", "tilted.csv", "run.json"):
        if (OUT / name).exists():
            raise FileExistsError(OUT / name)
    compiler = shutil.which("clang++")
    if compiler is None:
        raise RuntimeError("clang++ is required; this runner installs nothing")
    started = utc()
    wall_start = time.perf_counter()
    receipt = {
        "schema": "matching-one.p337-closed-source-finite-coupling.run.v1",
        "started_utc": started,
        "contract_commit": "b70dc4bd2fddd7676e9536b42bf912ee00ad302f",
        "contract_sha256": sha256(CONTRACT),
        "code_commit": git("rev-parse", "HEAD"),
        "branch": git("branch", "--show-current"),
        "git_status_before": git("status", "--short"),
        "backend_source": "b70dc4bd:scripts/p337_closed_source_exact.cpp",
        "cpp_sha256": sha256(CPP),
        "runner_sha256": sha256(Path(__file__).resolve()),
        "python_executable": sys.executable,
        "python_version": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "workers": 2,
        "new_random_samples": 0,
        "scores_computed": False,
        "histogram_schema": ["k", "g", "q", "count"],
        "zero_bins": "omitted; every nonzero integer bin is retained",
    }
    with tempfile.TemporaryDirectory(prefix="p337-finite-exact-") as temp:
        binary = Path(temp) / "enumerate"
        compile_command = [compiler, "-std=c++17", "-O3", str(CPP), "-o", str(binary)]
        compile_start = time.perf_counter()
        subprocess.run(compile_command, check=True, cwd=ROOT)
        receipt["compile"] = {
            "command": compile_command,
            "elapsed_seconds": time.perf_counter() - compile_start,
            "compiler_version": subprocess.check_output(
                [compiler, "--version"], text=True
            ).strip(),
            "binary_sha256": sha256(binary),
        }

        def run_one(job):
            name, a, b = job
            output = OUT / f"{name}.csv"
            command = [str(binary), str(a), str(b), str(output)]
            job_started = utc()
            job_start = time.perf_counter()
            process = subprocess.run(command, check=True, text=True, capture_output=True)
            elapsed = time.perf_counter() - job_start
            summary = json.loads(process.stdout)
            with output.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            count = sum(int(row["count"]) for row in rows)
            if count != 2**25 or summary["configurations"] != count:
                raise RuntimeError(f"incomplete enumeration: {name}: {count}")
            return {
                "name": name, "geometry": [a, b], "N": 25,
                "command": command, "started_utc": job_started,
                "finished_utc": utc(), "wall_seconds": elapsed,
                "exit_code": process.returncode,
                "stdout": summary, "stderr": process.stderr,
                "output": str(output.relative_to(ROOT)),
                "sha256": sha256(output), "nonzero_bins": len(rows),
                "count_sum": count,
            }

        with ThreadPoolExecutor(max_workers=2) as pool:
            receipt["jobs"] = list(pool.map(run_one, jobs))
    receipt["finished_utc"] = utc()
    receipt["wall_seconds"] = time.perf_counter() - wall_start
    receipt["status"] = "complete_no_scores"
    with (OUT / "run.json").open("x") as stream:
        json.dump(receipt, stream, indent=2)
        stream.write("\n")
    print(json.dumps({
        "receipt": str((OUT / "run.json").relative_to(ROOT)),
        "code_commit": receipt["code_commit"],
        "jobs": [{key: job[key] for key in
                  ("name", "count_sum", "nonzero_bins", "wall_seconds", "sha256")}
                 for job in receipt["jobs"]],
    }, indent=2))


if __name__ == "__main__":
    main()
