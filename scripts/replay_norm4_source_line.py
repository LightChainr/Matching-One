#!/usr/bin/env python3
"""Add new winding-line/source products to the existing norm-4 marked samples."""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import replay_norm4_source_thermal as old

ROOT = old.ROOT
CPP = ROOT / "src/norm4_source_line_replay.cpp"
CONTRACT = ROOT / "analysis/norm4_source_line_contract.json"
OUTPUT = ROOT / "results/norm4-source-line"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if not 1 <= args.workers <= 6:
        raise ValueError("workers must be between1 and6")
    contract = json.loads(CONTRACT.read_text())
    receipt = OUTPUT / "run.json"
    outputs = {n: OUTPUT / "raw" / f"n{n}.csv" for n in contract["Ns"]}
    if receipt.exists() or any(path.exists() for path in outputs.values()):
        raise ValueError("line reobservation already exists; do not repeat or overwrite")
    if (contract["original_backend_commit"] != old.SOURCE_COMMIT
            or contract["original_archive_commit"] != old.METADATA_COMMIT):
        raise ValueError("line contract differs from immutable production backends")
    source_path = ROOT / contract["source_result"]
    source = json.loads(source_path.read_text())
    runs = []
    for n in contract["Ns"]:
        run = old.source_run(n)
        samples = contract["marked_permutations_by_N"][str(n)]
        interval = contract["endpoint_counter_interval" if n in (260, 340) else "cyclic_counter_interval"]
        if (run["seed"] != contract["seeds"][str(n)]
                or source["by_N"][str(n)]["marked_samples"] != samples
                or interval != [run["old_counter_interval"][0], run["old_counter_interval"][0] + samples]):
            raise ValueError(f"N{n}: source samples or counters differ from frozen line selection")
        run.update(reobserved_permutations=samples, replayed_counter_interval=interval,
                   permutations_per_batch=samples // 100,
                   batch_rule=contract["endpoint_batch_rule" if n in (260, 340) else "cyclic_batch_rule"])
        runs.append(run)
    backend_bytes = {kind: old.git_bytes(old.SOURCE_COMMIT, path) for kind, path in old.BACKENDS.items()}
    record = {
        "schema": "matching-one.norm4-source-line-run.v1", "status": "running",
        "started_utc": old.utc_now(), "contract": contract,
        "execution_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "saved_source_result": {"path": old.display_path(source_path), "sha256": old.sha(source_path.read_bytes())},
        "source_runs": runs,
        "source_backends": {kind: {"commit": old.SOURCE_COMMIT, "path": old.BACKENDS[kind], "sha256": old.sha(raw)}
                            for kind, raw in backend_bytes.items()},
        "code": [{"path": old.display_path(path), "sha256": old.sha(path.read_bytes())}
                 for path in (Path(__file__).resolve(), CPP, CONTRACT, Path(old.__file__))],
        "environment": {"python": platform.python_version(), "machine": platform.machine(),
                        "platform": platform.platform(),
                        "compiler": subprocess.check_output(["c++", "--version"], text=True).splitlines()[0]},
        "workers": args.workers, "new_samples": 0,
        "old_permutation_reobservations": sum(row["reobserved_permutations"] for row in runs),
        "new_information": "Physical rank1-line O4 and bulk source joint products, absent from the existing source profiles",
        "sample_boundary": "Same marked100k cyclic/1M endpoint samples and original union batches; not independent new evidence groups",
        "server_actions": 0, "test_suites": [], "compile": [], "runs": []}
    OUTPUT.joinpath("raw").mkdir(parents=True, exist_ok=True)
    with receipt.open("x") as handle:
        handle.write(json.dumps(record, indent=2) + "\n")

    def save():
        receipt.write_text(json.dumps(record, indent=2) + "\n")

    started = time.perf_counter()
    try:
        with tempfile.TemporaryDirectory(prefix="matching-norm4-source-line-") as directory:
            build = Path(directory)
            binaries = {}
            for kind, raw in backend_bytes.items():
                archive = build / f"archived_{kind}.cpp"
                archive.write_bytes(raw)
                binary = build / f"line-replay-{kind}"
                command = ["c++", "-O3", "-std=c++17", f'-DMATCHING_NORM4_BACKEND="{archive}"']
                if kind == "integer_period":
                    command.append("-DMATCHING_NORM4_INTEGER=1")
                command.extend([str(CPP), "-o", str(binary)])
                begin = time.perf_counter()
                process = subprocess.run(command, cwd=ROOT, check=True, text=True, capture_output=True)
                record["compile"].append({"backend_kind": kind, "command": command,
                    "elapsed_seconds": time.perf_counter() - begin, "stdout": process.stdout,
                    "stderr": process.stderr, "binary_sha256": old.sha(binary.read_bytes())})
                binaries[kind] = binary
            save()

            def replay(run):
                n = run["N"]
                command = [str(binaries[run["backend_kind"]]), str(n), str(outputs[n])]
                begin = time.perf_counter()
                process = subprocess.run(command, cwd=ROOT, check=True, text=True, capture_output=True)
                return {"N": n, "command": command, "elapsed_seconds": time.perf_counter() - begin,
                        "stdout": process.stdout, "stderr": process.stderr,
                        "output": old.display_path(outputs[n]), "sha256": old.sha(outputs[n].read_bytes()),
                        "reobserved_permutations": run["reobserved_permutations"],
                        "replayed_counter_interval": run["replayed_counter_interval"],
                        "seed": run["seed"], "new_samples": 0}

            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                for future in as_completed([pool.submit(replay, run) for run in runs]):
                    row = future.result()
                    record["runs"].append(row)
                    record["runs"].sort(key=lambda item: item["N"])
                    save()
                    print(row["stdout"].strip(), flush=True)
        record["status"] = "completed"
    except Exception as error:
        record["status"] = "failed"
        record["error"] = repr(error)
        if isinstance(error, subprocess.CalledProcessError):
            record["failure_output"] = {"stdout": error.stdout, "stderr": error.stderr}
        raise
    finally:
        record["elapsed_seconds"] = time.perf_counter() - started
        record["finished_utc"] = old.utc_now()
        save()
    print(f"Completed new line/source marks on2.4M old permutations in {record['elapsed_seconds']:.2f}s", flush=True)


if __name__ == "__main__":
    main()
