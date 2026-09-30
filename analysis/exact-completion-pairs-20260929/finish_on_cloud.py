#!/usr/bin/env python3
"""Finish only missing source batches on TV2N0X after cloud-first rerouting."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import json
import os
import platform
import shlex
import subprocess
import tempfile
import time
from pathlib import Path

from run import HERE, OLD, sha


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-commit", required=True)
    args = ap.parse_args()
    manifest = HERE / "data/run.json"
    run = json.loads(manifest.read_text())
    assert run["status"] != "completed"
    for name, digest in run["source_sha256"].items():
        assert sha(HERE / name) == digest
    original = json.loads((OLD / "data/run.json").read_text())
    assert sha(OLD / "data/run.json") == run["original_manifest_sha256"]
    finished = {entry["batch"] for entry in run["batches"]}
    for entry in run["batches"]:
        assert sha(HERE / "data" / entry["file"]) == entry["sha256"]
    remaining = [entry for entry in original["batches"] if entry["batch"] not in finished]
    assert len(remaining) == 6 and finished == set(range(8))
    receipt = dict(instance="DevEnvC_TV2N0X", instance_id="4a8d1d443419434889e49148ed0a7ba6",
                   source_commit=args.source_commit, runner_sha256=sha(Path(__file__)),
                   pid=os.getpid(), platform=platform.platform(), python=platform.python_version(),
                   compiler=subprocess.check_output(["g++", "--version"], text=True).strip(),
                   utc_started=datetime.now(timezone.utc).isoformat(), workers=6,
                   retained_local_batches=sorted(finished), completed_cloud_batches=[],
                   reason="Owner requests medium/large computation on Huawei first; local incomplete workers stopped",
                   new_random_prefixes=0)
    run["cloud_continuation"] = receipt
    def save():
        manifest.write_text(json.dumps(run, indent=2)+"\n")
    start = time.monotonic()
    try:
        with tempfile.TemporaryDirectory(prefix="exact-pair-cloud-", dir=HERE) as temp:
            temp = Path(temp)
            binary = temp / "engine"
            subprocess.run(["g++", "-O3", "-std=c++17", str(HERE / "engine.cpp"), "-o", str(binary)], check=True)
            receipt["binary_sha256"] = sha(binary)
            # Same old-prefix seed; cost check only, never extra evidence.
            benchmark = [str(binary), "512", "16", remaining[0]["seed"], str(run["contract"]["b"]), str(temp / "benchmark.csv")]
            receipt["benchmark_log"] = subprocess.check_output(benchmark, text=True).strip()
            run["status"] = "running"; save()
            def job(entry):
                batch = entry["batch"]
                raw = temp / f"exact-square-L512-b{batch:02d}.csv"
                target = HERE / "data" / (raw.name+".gz")
                assert not target.exists()
                command = [str(binary), "512", str(entry["samples"]), entry["seed"], str(run["contract"]["b"]), str(raw)]
                begin = time.monotonic()
                log = subprocess.check_output(command, text=True).strip()
                seconds = time.monotonic()-begin
                payload = raw.read_bytes()
                assert payload.count(b"\n") == entry["samples"]+1
                with target.open("xb") as out:
                    with gzip.GzipFile(fileobj=out, mode="wb", mtime=0) as gz:
                        gz.write(payload)
                print(f"batch={batch}: {log}", flush=True)
                return dict(batch=batch, samples=entry["samples"], seed=entry["seed"],
                            file=target.name, sha256=sha(target), original_file=entry["file"],
                            dependency_group=entry["dependency_group"], runtime_seconds=seconds,
                            engine_log=log, command=shlex.join(command), execution="TV2N0X")
            with ThreadPoolExecutor(max_workers=6) as pool:
                for future in as_completed([pool.submit(job, entry) for entry in remaining]):
                    entry = future.result(); run["batches"].append(entry)
                    run["batches"].sort(key=lambda e: e["batch"])
                    receipt["completed_cloud_batches"].append(entry["batch"]); save()
        run["status"] = "completed"
    except BaseException as error:
        run["status"] = "failed"; receipt["error"] = str(error)
        raise
    finally:
        receipt["wall_seconds"] = time.monotonic()-start
        run["utc_finished"] = datetime.now(timezone.utc).isoformat()
        run["timing_note"] = "Mixed execution: local successful batches0..7, cloud8..13. Owner pause and abandoned partial work excluded from successful batch timings; cloud wall includes compile/benchmark."
        save()
    print(f"Cloud continuation completed in {receipt['wall_seconds']:.3f}s", flush=True)


if __name__ == "__main__":
    main()
