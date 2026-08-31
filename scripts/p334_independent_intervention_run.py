#!/usr/bin/env python3
"""Dispatch the frozen twenty independent batches; print operational receipts only."""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import time

SEED = 202608311920334
FREEZE = "bc0a18c207e3b09f49ea6b6af6601471114d654a"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, choices=(325, 425), required=True)
    p.add_argument("--producer-commit", required=True)
    p.add_argument("--workers", type=int, default=14)
    p.add_argument("--root", type=Path, required=True)
    a = p.parse_args()
    root = a.root.resolve()
    out = root / "results" / f"N{a.n}"
    out.mkdir(parents=True, exist_ok=False)
    src = root / "src/p334_independent_normal_intervention.cpp"
    backend = root / "src/threshold_rank_integer_period_mc.cpp"
    binary = root / f"producer_N{a.n}"
    record = {"status": "compiling", "N": a.n, "freeze_commit": FREEZE,
              "producer_commit": a.producer_commit, "seed": SEED,
              "batches": 20, "prefixes_per_batch": 25000, "paired_reps": 8,
              "workers": a.workers, "hostname": platform.node(),
              "platform": platform.platform(), "python": platform.python_version(),
              "started_unix": time.time(), "completed_batches": [],
              "source_sha256": {str(f.relative_to(root)): hashlib.sha256(f.read_bytes()).hexdigest()
                                for f in (src, backend, Path(__file__))}}
    for label, path in (("cpu_quota", "/sys/fs/cgroup/cpu/cpu.cfs_quota_us"),
                        ("cpu_period", "/sys/fs/cgroup/cpu/cpu.cfs_period_us"),
                        ("memory_limit", "/sys/fs/cgroup/memory/memory.limit_in_bytes")):
        if Path(path).exists():
            record[label] = Path(path).read_text().strip()
    def save():
        (out / "run.json").write_text(json.dumps(record, indent=2) + "\n")
    def batch(b):
        cmd = [str(binary), "--n", str(a.n), "--batch", str(b), "--prefixes", "25000",
               "--reps", "8", "--seed", str(SEED), "--code-commit", a.producer_commit,
               "--output", str(out / f"batch-{b:02}")]
        t = time.perf_counter()
        r = subprocess.run(cmd, capture_output=True, text=True)
        (out / f"batch-{b:02}.log").write_text(r.stdout + r.stderr)
        return {"batch": b, "returncode": r.returncode, "seconds": time.perf_counter()-t,
                "command": cmd, "operational_output": r.stdout, "stderr": r.stderr}
    save()
    try:
        compile_cmd = ["g++", "-O3", "-std=c++17", str(src), "-lz", "-o", str(binary)]
        r = subprocess.run(compile_cmd, capture_output=True, text=True, check=True)
        record["compile"] = {"command": compile_cmd, "stdout": r.stdout, "stderr": r.stderr,
                             "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest()}
        record["status"] = "running"
        save()
        print(f"START N={a.n} workers={a.workers} fresh_prefixes=500000", flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
            for task in concurrent.futures.as_completed([pool.submit(batch, b) for b in range(20)]):
                item = task.result()
                record["completed_batches"].append(item)
                save()
                print(f"BATCH N={a.n} batch={item['batch']} exit={item['returncode']} seconds={item['seconds']:.3f}", flush=True)
        record["status"] = "completed" if all(x["returncode"] == 0 for x in record["completed_batches"]) else "failed"
    except Exception as error:
        record["status"] = "failed"
        record["error"] = repr(error)
        if isinstance(error, subprocess.CalledProcessError):
            record["error_stdout"], record["error_stderr"] = error.stdout, error.stderr
        raise
    finally:
        record["finished_unix"] = time.time()
        record["elapsed_seconds"] = record["finished_unix"]-record["started_unix"]
        save()
    print(f"FINAL N={a.n} status={record['status']} seconds={record['elapsed_seconds']:.3f}", flush=True)
    if record["status"] != "completed":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
