#!/usr/bin/env python3
"""One bounded cloud coupling experiment; standard-library runtime."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OLD = HERE.parent / "safe-insertion-production-20260929/data"
SPEC = importlib.util.spec_from_file_location("window_runtime", HERE.parent / "geometric-source-window-20260930/run.py")
UTIL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(UTIL)
sha, save, env = UTIL.sha, UTIL.save, UTIL.env
SOURCES = [HERE / name for name in ["engine.cpp", "run.py", "score.py", "contract.json"]] + [
    HERE.parent / "exact-completion-pairs-20260929/engine.cpp",
    HERE.parent / "completion-hazard-production-20260929/engine.cpp",
    HERE.parent / "geometric-source-window-20260930/run.py"]


def seed(batch):
    text = f"matching-one/safe-transport-selection-20260930/continuation/square/L512/{batch}"
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "big")


def command(batch, samples, output):
    return [str(HERE / "build/engine"), "512", str(samples), str(batch["seed"]),
            str(seed(batch["batch"])), "155385", str(output)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["benchmark", "run"])
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--toolchain", type=Path)
    args = parser.parse_args()
    original = json.loads((OLD / "run.json").read_text())
    batches = sorted(original["batches"], key=lambda x: x["batch"])
    assert original["status"] == "completed" and len(batches) == 14
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
            if args.toolchain is None:
                raise ValueError("explicit task-owned toolchain required")
            project = args.toolchain / "safe-transport-selection-project"
            project.mkdir(exist_ok=False)
            for src in SOURCES:
                if src.suffix == ".cpp":
                    dst = project / src.relative_to(ROOT)
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(src, dst)
            compile_cmd = ["chroot", str(args.toolchain), "/usr/bin/g++", "-O3", "-std=c++17",
                           "/safe-transport-selection-project/analysis/safe-transport-selection-20260930/engine.cpp",
                           "-o", "/safe-transport-selection-project/engine"]
            version_cmd = ["chroot", str(args.toolchain), "/usr/bin/g++", "--version"]
        subprocess.run(compile_cmd, check=True)
        if not compiler:
            shutil.copyfile(project / "engine", build / "engine")
            (build / "engine").chmod(0o755)
        output = build / "old-seed-32.csv"
        cmd = command(batches[0], 32, output)
        begin = time.monotonic()
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        seconds = time.monotonic()-begin
        with output.open() as f, gzip.open(OLD / batches[0]["file"], "rt") as g:
            rows = list(csv.DictReader(f))
            assert len(rows) == 32
            previous = csv.DictReader(g)
            for row in rows:
                assert row["J1"] == next(previous)["J1"]
        rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        workers = min(14, math.floor(environment["cpu_quota"]), int(environment["memory_bytes"]*.8//(rss*1024)))
        assert workers >= 1
        predicted = seconds/32*10000*math.ceil(14/workers)
        result = dict(environment=environment, command=cmd, build_command=compile_cmd,
                      compiler=subprocess.check_output(version_cmd, text=True), stdout=proc.stdout,
                      wall_seconds=seconds, children_peak_rss_kib=rss, workers=workers,
                      predicted_ideal_seconds=predicted, within_cost_budget=predicted<=1800,
                      rows=32, original_births_match=True, scientific_weight=0,
                      source_commit=args.source_commit, binary_sha256=sha(build / "engine"))
        save(HERE / "benchmark.json", result)
        print(json.dumps(result), flush=True)
        return
    bench = json.loads((HERE / "benchmark.json").read_text())
    assert bench["within_cost_budget"] and sha(HERE / "build/engine") == bench["binary_sha256"]
    workers = min(bench["workers"], math.floor(environment["cpu_quota"]), 14)
    data = HERE / "data"
    data.mkdir(exist_ok=False)
    record = dict(schema="matching-one.safe-transport-selection-run.v1", status="running", pid=os.getpid(),
                  started_utc=datetime.now(timezone.utc).isoformat(), environment=environment,
                  source_commit=args.source_commit, binary_sha256=bench["binary_sha256"], workers=workers,
                  source_sha256={str(p.relative_to(ROOT)): sha(p) for p in SOURCES},
                  original_manifest_sha256=sha(OLD / "run.json"), new_independent_prefixes=0, batches=[])
    save(data / "started.json", record)
    def one(batch):
        output = data / ("selection-b%02d.csv" % batch["batch"])
        cmd = command(batch, 10000, output)
        begin = time.monotonic()
        proc = subprocess.run(cmd, capture_output=True, text=True)
        row = dict(batch=batch["batch"], permutation_seed=batch["seed"], continuation_seed=str(seed(batch["batch"])),
                   samples=10000, original_file=batch["file"], original_sha256=batch["sha256"],
                   dependency_group=batch["dependency_group"], command=cmd, returncode=proc.returncode,
                   stdout=proc.stdout, stderr=proc.stderr, wall_seconds=time.monotonic()-begin)
        if not proc.returncode:
            compressed = output.with_suffix(".csv.gz")
            with output.open("rb") as f, compressed.open("xb") as g:
                with gzip.GzipFile(filename="", mode="wb", fileobj=g, mtime=0) as z:
                    shutil.copyfileobj(f, z)
            row.update(file=compressed.name, sha256=sha(compressed))
        save(data / ("batch-%02d.json" % batch["batch"]), row)
        if proc.returncode:
            raise RuntimeError(row)
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
