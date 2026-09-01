#!/usr/bin/env python3
"""Mark lag-one cluster-source/topology events on fixed OLD norm-4 permutations.

Default outputs are raw/n{N}.csv.gz and run.json under results/norm4-lagged-source.
Each endpoint row already combines its old1000 and added9000 analysis blocks.
No counter, seed, lag, sample-size or geometry override is exposed.
"""
from __future__ import annotations

import argparse
import gzip
import json
import platform
import shutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import replay_norm4_source_thermal as old

ROOT = old.ROOT
CPP = ROOT / "src/norm4_lagged_source_replay.cpp"
OUTPUT = ROOT / "results/norm4-lagged-source"
SOURCE_RESULT = ROOT / "results/norm4-source-endpoint-1m/latest.json"
NS = (65, 85, 130, 170, 260, 340)
COLUMNS = ["n", "a", "b", "orientation", "batch", "k", "samples",
           "event_count01", "event_count02", "event_count12",
           "sum_s_previous01", "sum_s_previous02", "sum_s_previous12"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=4, help="concurrent fixed-counter reobservations, 1..6")
    args = parser.parse_args()
    if not 1 <= args.workers <= 6:
        raise ValueError("workers must be between 1 and 6")
    receipt = OUTPUT / "run.json"
    outputs = {n: OUTPUT / "raw" / f"n{n}.csv.gz" for n in NS}
    if receipt.exists() or any(path.exists() for path in outputs.values()):
        raise ValueError("lagged output already exists; refusing to repeat or overwrite")
    if (old.SOURCE_COMMIT != "bfab0330f5f56ca4d746b45d737f1607e3d229a0"
            or old.METADATA_COMMIT != "8b26a30a785bc142a9d17bfed99a8d0e98ddc4dc"):
        raise ValueError("original production helper pins changed")
    source_bytes = SOURCE_RESULT.read_bytes()
    source = json.loads(source_bytes)
    runs = []
    for n in NS:
        run = old.source_run(n)
        endpoint = n in (260, 340)
        samples = 1000000 if endpoint else 100000
        if source["by_N"][str(n)]["marked_samples"] != samples:
            raise ValueError(f"N{n}: saved source sample count differs from fixed lagged selection")
        begin = run["old_counter_interval"][0]
        run.update(
            reobserved_permutations=samples,
            replayed_counter_interval=[begin, begin + samples],
            permutations_per_batch=samples // 100,
            batch_segments=([{"counter_start": begin, "stride": 1000, "count": 1000},
                             {"counter_start": begin + 100000, "stride": 9000, "count": 9000}]
                            if endpoint else [{"counter_start": begin, "stride": 1000, "count": 1000}]),
            batch_rule="For b=0..99, union the intervals [counter_start+stride*b, +count) from batch_segments",
            dependency_group=(f"norm4-N{n}-seed{run['seed']}-first100k" if endpoint
                              else "norm4-cyclic-seed2026104501-first100k"),
        )
        runs.append(run)
    backend_bytes = {kind: old.git_bytes(old.SOURCE_COMMIT, path) for kind, path in old.BACKENDS.items()}
    record = {
        "schema": "matching-one.norm4-lagged-source-run.v1", "status": "running",
        "started_utc": old.utc_now(),
        "execution_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_runs": runs,
        "saved_source_result": {"path": old.display_path(SOURCE_RESULT), "sha256": old.sha(source_bytes)},
        "source_backends": {kind: {"commit": old.SOURCE_COMMIT, "path": old.BACKENDS[kind], "sha256": old.sha(raw)}
                            for kind, raw in backend_bytes.items()},
        "code": [{"path": old.display_path(path), "sha256": old.sha(path.read_bytes())}
                 for path in (Path(__file__).resolve(), CPP, Path(old.__file__).resolve())],
        "environment": {"python": platform.python_version(), "machine": platform.machine(),
                        "platform": platform.platform(),
                        "compiler": subprocess.check_output(["c++", "--version"], text=True).splitlines()[0]},
        "workers": args.workers, "new_samples": 0,
        "old_permutation_reobservations": sum(run["reobserved_permutations"] for run in runs),
        "readouts": {
            "columns": COLUMNS, "K": "0..N inclusive; every row has the full batch denominator",
            "events": {"01": "K=Kminus<Kplus", "02": "K=Kminus=Kplus", "12": "K=Kplus>Kminus"},
            "source": "s_previous=CB(K-1)+CW(N-K+1), same permutation, bulk integer, no division by N",
            "K0": "all six event/source fields zero",
            "centering": "not performed in replay; estimate E[s_previous|K-1,early rank,g] from the paired saved source profiles",
            "q_reconstruction": "-samples+cumsum(event_count01+2*event_count02+event_count12)",
            "E_reconstruction": "samples+cumsum(-event_count01+event_count12)",
        },
        "compression": "CSV generated in temporary build directory, then gzip level6, empty filename and mtime=0; receipt hashes compressed and uncompressed bytes",
        "sample_boundary": "Same100k cyclic/1M endpoint permutations and source union batches; new joint marks, not new independent samples or full1.9B/1B source precision",
        "dependency": "Both directions paired; N65/85/130/170 share seed/counter batches; N260 and N340 each retain a separate original seed group. Historical group IDs do not change the recorded endpoint1M precision.",
        "server_actions": 0, "test_suites": [], "compile": [], "runs": [],
    }
    OUTPUT.joinpath("raw").mkdir(parents=True, exist_ok=True)
    with receipt.open("x") as handle:
        handle.write(json.dumps(record, indent=2) + "\n")

    def save() -> None:
        receipt.write_text(json.dumps(record, indent=2) + "\n")

    started = time.perf_counter()
    try:
        with tempfile.TemporaryDirectory(prefix="matching-norm4-lagged-source-") as directory:
            build = Path(directory)
            binaries = {}
            for kind, raw in backend_bytes.items():
                archive = build / f"archived_{kind}.cpp"
                archive.write_bytes(raw)
                binary = build / f"lagged-replay-{kind}"
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

            def replay(run: dict) -> dict:
                n = run["N"]
                csv = build / f"n{n}.csv"
                command = [str(binaries[run["backend_kind"]]), str(n), str(csv)]
                begin = time.perf_counter()
                process = subprocess.run(command, cwd=ROOT, check=True, text=True, capture_output=True)
                replay_seconds = time.perf_counter() - begin
                raw_sha256 = old.sha(csv.read_bytes())
                compression_begin = time.perf_counter()
                with csv.open("rb") as src, outputs[n].open("xb") as dst:
                    with gzip.GzipFile(filename="", mode="wb", fileobj=dst, compresslevel=6, mtime=0) as compressed:
                        shutil.copyfileobj(src, compressed)
                return {"N": n, "command": command, "elapsed_seconds": replay_seconds,
                        "compression_seconds": time.perf_counter() - compression_begin,
                        "stdout": process.stdout, "stderr": process.stderr,
                        "output": old.display_path(outputs[n]), "sha256": old.sha(outputs[n].read_bytes()),
                        "uncompressed_sha256": raw_sha256, "compression": "gzip",
                        "compressed_bytes": outputs[n].stat().st_size, "uncompressed_bytes": csv.stat().st_size,
                        "reobserved_permutations": run["reobserved_permutations"],
                        "permutations_per_batch": run["permutations_per_batch"],
                        "replayed_counter_interval": run["replayed_counter_interval"],
                        "batch_segments": run["batch_segments"], "seed": run["seed"], "new_samples": 0}

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
    print(f"Completed lag-one marks on {record['old_permutation_reobservations']} OLD permutations in "
          f"{record['elapsed_seconds']:.2f}s; new samples=0", flush=True)


if __name__ == "__main__":
    main()
