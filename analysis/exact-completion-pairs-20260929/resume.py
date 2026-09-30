#!/usr/bin/env python3
"""Finish missing whole batches after the owner's pause; retain completed rows."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import time

from run import HERE, OLD, sha


def main():
    path = HERE / "data/run.json"
    run = json.loads(path.read_text())
    if run["status"] == "completed":
        raise SystemExit("Replay already complete")
    for pid in [run["pid"]] + [v["pid"] for v in run.get("resumptions", [])]:
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            continue
        raise SystemExit(f"Recorded process {pid} is alive; inspect it rather than launching a duplicate")
    for name, digest in run["source_sha256"].items():
        assert sha(HERE / name) == digest, f"changed source: {name}"
    source = json.loads((OLD / "data/run.json").read_text())
    assert sha(OLD / "data/run.json") == run["original_manifest_sha256"]
    finished = {entry["batch"] for entry in run["batches"]}
    for entry in run["batches"]:
        assert sha(HERE / "data" / entry["file"]) == entry["sha256"]
    remaining = [entry for entry in source["batches"] if entry["batch"] not in finished]
    begin = time.monotonic()
    receipt = {"pid": os.getpid(), "utc_started": datetime.now(timezone.utc).isoformat(),
               "reason": "Owner paused the task; stopped processes no longer existed on resume",
               "retained_batches": sorted(finished), "replayed_batches": [e["batch"] for e in remaining],
               "source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
               "source_sha256": sha(Path(__file__)), "workers": 4,
               "new_random_prefixes": 0, "partial_files": "Original abandoned temporary CSVs retained; no partial batches scored"}
    run.setdefault("resumptions", []).append(receipt)
    def save():
        path.write_text(json.dumps(run, indent=2)+"\n")
    run["status"] = "running"; save()
    try:
        with tempfile.TemporaryDirectory(prefix="matching-exact-pairs-resume-") as temp:
            temp = Path(temp)
            binary = temp / "engine"
            subprocess.run(["clang++", "-O3", "-std=c++17", str(HERE / "engine.cpp"), "-o", str(binary)], check=True)
            receipt["binary_sha256"] = sha(binary)
            def job(entry):
                batch = entry["batch"]
                raw = temp / f"exact-square-L512-b{batch:02d}.csv"
                target = HERE / "data" / (raw.name+".gz")
                if target.exists():
                    raise RuntimeError(f"Unregistered output already exists: {target}")
                command = [str(binary), str(run["contract"]["L"]), str(entry["samples"]),
                           entry["seed"], str(run["contract"]["b"]), str(raw)]
                start = time.monotonic()
                log = subprocess.check_output(command, text=True).strip()
                seconds = time.monotonic()-start
                payload = raw.read_bytes()
                assert payload.count(b"\n") == entry["samples"]+1
                with target.open("xb") as out:
                    with gzip.GzipFile(fileobj=out, mode="wb", mtime=0) as gz:
                        gz.write(payload)
                print(f"batch={batch}: {log}", flush=True)
                return dict(batch=batch, samples=entry["samples"], seed=entry["seed"],
                            file=target.name, sha256=sha(target), original_file=entry["file"],
                            dependency_group=entry["dependency_group"], runtime_seconds=seconds,
                            engine_log=log, command=shlex.join(command))
            with ThreadPoolExecutor(max_workers=4) as pool:
                for future in as_completed([pool.submit(job, entry) for entry in remaining]):
                    run["batches"].append(future.result())
                    run["batches"].sort(key=lambda e: e["batch"]); save()
        run["status"] = "completed"
    except BaseException as error:
        run["status"] = "failed"; receipt["error"] = str(error)
        raise
    finally:
        receipt["wall_seconds"] = time.monotonic()-begin
        run["utc_finished"] = datetime.now(timezone.utc).isoformat()
        run["elapsed_seconds_including_owner_pause"] = (datetime.now(timezone.utc)-datetime.fromisoformat(run["utc_started"])).total_seconds()
        run["timing_note"] = "Per-batch times cover successful evaluations only; elapsed includes owner pause. Abandoned partial work is not a second sample."
        save()
    print(f"Completed remaining batches in {receipt['wall_seconds']:.3f}s", flush=True)


if __name__ == "__main__":
    main()
