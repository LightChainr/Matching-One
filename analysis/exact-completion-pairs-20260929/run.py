#!/usr/bin/env python3
"""Replay existing prefixes and integrate out safe-site probe randomness.

This is not an independent acquisition. All seed/sample identities come
from the completed safe-insertion block. No old files are overwritten.
"""
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
OLD = HERE.parent / "safe-insertion-production-20260929"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-commit", required=True)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--compiler", default="clang++")
    args = ap.parse_args()
    if not 1 <= args.workers <= 4:
        raise SystemExit("This local replay uses at most four workers")
    source = OLD / "data/run.json"
    original = json.loads(source.read_text())
    assert original["status"] == "completed"
    contract = original["contract"]
    directory = HERE / "data"
    directory.mkdir(exist_ok=True)
    manifest = directory / "run.json"
    if manifest.exists():
        raise SystemExit("Do not overwrite or silently resume an existing replay")
    started = time.monotonic()
    run = {"schema": "matching-one.exact-pair-replay.v1", "status": "building",
           "source_commit": args.source_commit, "pid": os.getpid(),
           "source_sha256": {name: sha(HERE / name) for name in ("engine.cpp", "run.py")},
           "included_geometry_sha256": sha(HERE.parent / "completion-hazard-production-20260929/engine.cpp"),
           "original_manifest": "analysis/safe-insertion-production-20260929/data/run.json",
           "original_manifest_sha256": sha(source), "original_acquisition_commit": original["source_commit"],
           "contract": contract, "new_random_prefixes": 0,
           "interpretation": "Same source prefixes, exact conditional mean of safe-insertion probes. Post-readout estimator refinement, not new independent evidence.",
           "workers": args.workers, "platform": platform.platform(), "python": platform.python_version(),
           "command": shlex.join(["python3", str(Path(__file__).resolve()), *os.sys.argv[1:]]),
           "compiler": subprocess.check_output([args.compiler, "--version"], text=True).strip(),
           "utc_started": datetime.now(timezone.utc).isoformat(), "batches": []}
    def save():
        manifest.write_text(json.dumps(run, indent=2)+"\n")
    save()
    try:
        with tempfile.TemporaryDirectory(prefix="exact-pair-build-", dir=HERE) as temp:
            temp = Path(temp)
            binary = temp / "engine"
            build = [args.compiler, "-O3", "-std=c++17", str(HERE / "engine.cpp"), "-o", str(binary)]
            subprocess.run(build, check=True)
            run["build_command"] = shlex.join(build);run["binary_sha256"] = sha(binary)
            # Physical witness control and first128 cost/source check were
            # already done once locally before this full replay, not repeated.
            def job(entry):
                batch = entry["batch"]
                raw = temp / f"exact-square-L512-b{batch:02d}.csv"
                command = [str(binary), str(contract["L"]), str(entry["samples"]),
                           entry["seed"], str(contract["b"]), str(raw)]
                begin = time.monotonic()
                log = subprocess.check_output(command, text=True).strip()
                seconds = time.monotonic()-begin
                data = raw.read_bytes()
                if data.count(b"\n") != entry["samples"]+1:
                    raise RuntimeError("incomplete replay")
                target = directory / (raw.name+".gz")
                with target.open("wb") as out:
                    with gzip.GzipFile(fileobj=out, mode="wb", mtime=0) as gz:
                        gz.write(data)
                print(f"batch={batch}: {log}", flush=True)
                return {"batch": batch, "samples": entry["samples"], "seed": entry["seed"],
                        "file": target.name, "sha256": sha(target), "original_file": entry["file"],
                        "dependency_group": entry["dependency_group"], "runtime_seconds": seconds,
                        "engine_log": log, "command": shlex.join(command)}
            run["status"] = "running";save()
            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                futures = [pool.submit(job, entry) for entry in original["batches"]]
                for future in as_completed(futures):
                    run["batches"].append(future.result())
                    run["batches"].sort(key=lambda x: x["batch"]);save()
        run["status"] = "completed"
    except BaseException as error:
        run["status"] = "failed";run["error"] = str(error)
        raise
    finally:
        run["wall_seconds"] = time.monotonic()-started
        run["utc_finished"] = datetime.now(timezone.utc).isoformat();save()
    print(f"Completed same-prefix replay in {run['wall_seconds']:.3f}s", flush=True)


if __name__ == "__main__":
    main()
