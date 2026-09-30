#!/usr/bin/env python3
"""Bounded Huawei replay of one fixed mechanism readout, standard library."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OLD = HERE.parent / "safe-insertion-production-20260929/data"
SOURCES = [HERE / "engine.cpp", HERE / "run.py", HERE / "contract.json",
           HERE.parent / "exact-completion-pairs-20260929/engine.cpp",
           HERE.parent / "completion-hazard-production-20260929/engine.cpp"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    with path.open("x") as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write("\n")


def env():
    cg = Path("/sys/fs/cgroup")
    if (cg / "cpu/cpu.cfs_quota_us").exists():
        cores = int((cg / "cpu/cpu.cfs_quota_us").read_text())/int((cg / "cpu/cpu.cfs_period_us").read_text())
        memory = int((cg / "memory/memory.limit_in_bytes").read_text())
    else:
        quota, period = (cg / "cpu.max").read_text().split()
        cores, memory = int(quota)/int(period), int((cg / "memory.max").read_text())
    if cores <= 0 or memory <= 0:
        raise ValueError("unmeasured limits")
    return dict(platform=platform.platform(), machine=platform.machine(), python=platform.python_version(),
                cpu_quota=cores, memory_bytes=memory, hostname=platform.node())


def command(batch, samples, output):
    return [str(HERE / "build/engine"), "512", str(samples), str(batch["seed"]), "155385", "735", str(output)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["benchmark", "run"])
    parser.add_argument("--source-commit", required=True)
    args = parser.parse_args()
    manifest = json.loads((OLD / "run.json").read_text())
    batches = sorted(manifest["batches"], key=lambda x: x["batch"])
    assert manifest["status"] == "completed" and len(batches) == 14
    assert all(b["samples"] == 10000 for b in batches)
    environment = env()
    if args.mode == "benchmark":
        build = HERE / "build"
        build.mkdir(exist_ok=False)
        compiler = shutil.which("g++")
        if compiler:
            compile_cmd = [compiler, "-O3", "-std=c++17", str(HERE / "engine.cpp"), "-o", str(build / "engine")]
            version_cmd = [compiler, "--version"]
        else:
            toolchain = ROOT / "toolchain"
            project = toolchain / "project"
            project.mkdir(exist_ok=False)
            for src in SOURCES:
                if src.suffix != ".cpp":
                    continue
                dst = project / src.relative_to(ROOT)
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(src, dst)
            compile_cmd = ["chroot", str(toolchain), "/usr/bin/g++", "-O3", "-std=c++17",
                           "/project/analysis/geometric-source-window-20260930/engine.cpp", "-o", "/project/engine"]
            version_cmd = ["chroot", str(toolchain), "/usr/bin/g++", "--version"]
        subprocess.run(compile_cmd, check=True)
        if not compiler:
            shutil.copyfile(project / "engine", build / "engine")
            (build / "engine").chmod(0o755)
        output = build / "old-seed-32.csv"
        cmd = command(batches[0], 32, output)
        start = time.monotonic()
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        seconds = time.monotonic()-start
        with output.open() as f, gzip.open(OLD / batches[0]["file"], "rt") as g:
            new, old = csv.DictReader(f), csv.DictReader(g)
            keys = ["J1", "rank_b", "dx_b", "dy_b", "nu_b"]
            rows = list(new)
            assert len(rows) == 32
            for row in rows:
                previous = next(old)
                assert all(row[k] == previous[k] for k in keys)
        rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        workers = min(14, math.floor(environment["cpu_quota"]), int(environment["memory_bytes"]*.8//(rss*1024)))
        assert workers >= 1
        result = dict(environment=environment, command=cmd, build_command=compile_cmd,
                      compiler=subprocess.check_output(version_cmd, text=True), stdout=proc.stdout,
                      wall_seconds=seconds, children_peak_rss_kib=rss, workers=workers,
                      predicted_ideal_seconds=seconds/32*10000*math.ceil(14/workers),
                      rows=32, aligned_first5=True, scientific_weight=0, source_commit=args.source_commit,
                      binary_sha256=sha(build / "engine"))
        save(HERE / "benchmark.json", result)
        print(json.dumps(result), flush=True)
        return
    bench = json.loads((HERE / "benchmark.json").read_text())
    assert sha(HERE / "build/engine") == bench["binary_sha256"]
    workers = min(bench["workers"], math.floor(environment["cpu_quota"]), 14)
    data = HERE / "data"
    data.mkdir(exist_ok=False)
    record = dict(schema="matching-one.geometric-source-window-run.v1", status="running", pid=os.getpid(),
                  started_utc=datetime.now(timezone.utc).isoformat(), environment=environment,
                  source_commit=args.source_commit, binary_sha256=bench["binary_sha256"], workers=workers,
                  source_sha256={str(p.relative_to(ROOT)): sha(p) for p in SOURCES},
                  original_manifest_sha256=sha(OLD / "run.json"), new_random_prefixes=0, batches=[])
    save(data / "started.json", record)
    def one(batch):
        output = data / ("window-b%02d.csv" % batch["batch"])
        cmd = command(batch, 10000, output)
        begin = time.monotonic()
        process = subprocess.run(cmd, capture_output=True, text=True)
        row = dict(batch=batch["batch"], seed=batch["seed"], samples=10000, original_file=batch["file"],
                   dependency_group=batch["dependency_group"], command=cmd, returncode=process.returncode,
                   stdout=process.stdout, stderr=process.stderr, wall_seconds=time.monotonic()-begin)
        save(data / ("batch-%02d.json" % batch["batch"]), row)
        if process.returncode:
            raise RuntimeError(row)
        compressed = output.with_suffix(".csv.gz")
        with output.open("rb") as f, compressed.open("xb") as g:
            with gzip.GzipFile(filename="", mode="wb", fileobj=g, mtime=0) as z:
                shutil.copyfileobj(f, z)
        row.update(file=compressed.name, sha256=sha(compressed))
        print(json.dumps(dict(batch=row["batch"], status="completed", seconds=row["wall_seconds"])), flush=True)
        return row
    begin = time.monotonic()
    try:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for task in as_completed([pool.submit(one, b) for b in batches]):
                record["batches"].append(task.result())
        record["batches"].sort(key=lambda b: b["batch"])
        record["status"] = "completed"
    except Exception as exc:
        record.update(status="failed", error=str(exc))
        raise
    finally:
        record.update(wall_seconds=time.monotonic()-begin, finished_utc=datetime.now(timezone.utc).isoformat())
        save(data / "run.json", record)


if __name__ == "__main__":
    main()
