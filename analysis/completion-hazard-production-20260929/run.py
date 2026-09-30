#!/usr/bin/env python3
"""Run one predeclared geometry-marked production block; Python stdlib only."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import subprocess
import tempfile
import time

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def seed(role, batch, stream):
    name = f"matching-one/completion-independent-20260929/{role}/square/L512/{batch}/{stream}"
    return int.from_bytes(hashlib.sha256(name.encode()).digest()[:8], "big")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-commit", required=True)
    ap.add_argument("--compiler", default="g++")
    ap.add_argument("--output", type=Path, default=HERE / "data")
    ap.add_argument("--budget-seconds", type=float, default=900)
    args = ap.parse_args()
    contract = json.loads((HERE / "contract.json").read_text())
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = args.output / "run.json"
    if manifest.exists():
        raise SystemExit("Existing run is not overwritten or silently resumed")
    workers = contract["workers"]
    started = time.monotonic()
    run = {
        "schema": "matching-one.completion-hazard-run.v1",
        "status": "building", "pid": os.getpid(),
        "utc_started": datetime.now(timezone.utc).isoformat(),
        "source_commit": args.source_commit,
        "source_sha256": {name: sha(HERE / name) for name in ("engine.cpp", "run.py", "contract.json")},
        "contract": contract,
        "command": shlex.join(["python3", str(Path(__file__).resolve()), *os.sys.argv[1:]]),
        "platform": platform.platform(), "machine": platform.machine(),
        "python": platform.python_version(),
        "compiler": subprocess.check_output([args.compiler, "--version"], text=True).strip(),
        "workers": workers, "budget_seconds": args.budget_seconds,
        "seed_rule": "first 64 SHA256 bits of matching-one/completion-independent-20260929/{role}/square/L512/{batch}/{stream}; stream=permutation or audit-time",
        "batches": [], "benchmarks": [],
    }
    def save():
        manifest.write_text(json.dumps(run, indent=2)+"\n")
    for name in ("cpu/cpu.cfs_quota_us", "cpu/cpu.cfs_period_us", "memory/memory.limit_in_bytes", "cpu.max", "memory.max"):
        path = Path("/sys/fs/cgroup") / name
        if path.exists():
            run.setdefault("cgroup", {})[name] = path.read_text().strip()
    save()
    try:
        with tempfile.TemporaryDirectory(prefix="completion-build-", dir=HERE) as tmp:
            tmp = Path(tmp)
            binary = tmp / "engine"
            build = [args.compiler, "-O3", "-std=c++17", str(HERE / "engine.cpp"), "-o", str(binary)]
            subprocess.run(build, check=True)
            run["build_command"] = shlex.join(build)
            run["binary_sha256"] = sha(binary)

            def job(role, batch, samples):
                stem = f"{role}-square-L512-b{batch:02d}"
                raw = tmp / (stem+".csv")
                permutation_seed, audit_seed = seed(role, batch, "permutation"), seed(role, batch, "audit-time")
                command = [str(binary), "sample", "square", str(contract["L"]), str(samples),
                           str(permutation_seed), str(audit_seed), str(contract["b"]), str(contract["c"]), str(raw)]
                begin = time.monotonic()
                log = subprocess.check_output(command, text=True).strip()
                seconds = time.monotonic()-begin
                target = args.output / (stem+".csv.gz")
                content = raw.read_bytes()
                if content.count(b"\n") != samples+1:
                    raise RuntimeError("incomplete sample output")
                with target.open("wb") as output:
                    with gzip.GzipFile(fileobj=output, mode="wb", mtime=0) as gz:
                        gz.write(content)
                print(f"{role} batch={batch}: {log}", flush=True)
                return {"file": target.name, "batch": batch, "samples": samples,
                        "seed": str(permutation_seed), "time_seed": str(audit_seed),
                        "role": role, "runtime_seconds": seconds, "engine_log": log,
                        "command": shlex.join(command), "sha256": sha(target),
                        "dependency_group": f"completion-independent-20260929/square/L512/batch{batch}" if role=="production" else "excluded-benchmark"}

            run["status"] = "benchmarking"; save()
            benchmark = job("benchmark", 0, 32)
            run["benchmarks"].append(benchmark)
            estimate = benchmark["runtime_seconds"]/32 * contract["samples_per_batch"]*contract["batches"]/workers
            run["estimated_production_wall_seconds"] = estimate
            print(f"Estimated ideal-parallel production wall: {estimate:.2f}s", flush=True)
            if estimate > args.budget_seconds:
                run["status"] = "benchmark-budget-stop"; save()
                return
            run["status"] = "running"; save()
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = [pool.submit(job, "production", b, contract["samples_per_batch"])
                           for b in range(contract["batches"])]
                for future in as_completed(futures):
                    run["batches"].append(future.result())
                    run["batches"].sort(key=lambda x: x["batch"])
                    save()
        run["status"] = "completed"
    except BaseException as error:
        run["status"] = "failed"
        run["error"] = str(error)
        raise
    finally:
        run["wall_seconds"] = time.monotonic()-started
        run["utc_finished"] = datetime.now(timezone.utc).isoformat()
        save()
    print(f"Completed {len(run['batches'])} batches in {run['wall_seconds']:.2f}s", flush=True)


if __name__ == "__main__":
    main()
