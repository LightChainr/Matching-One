#!/usr/bin/env python3
"""Reobserve only the frozen first 100k old permutations of both norm-4 chains.

The two immutable production backends are compiled separately: cyclic labels
for N65/85/130/170 and the original HNF labels for N260/340.  This is not a
fixed-p importance sample: each old permutation contributes at every K.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "analysis/p40_source_thermal_chain_candidates.json"
SOURCE_COMMIT = "bfab0330f5f56ca4d746b45d737f1607e3d229a0"
METADATA_COMMIT = "8b26a30a785bc142a9d17bfed99a8d0e98ddc4dc"
METADATA_DIRECTORY = "results/server-20260829/P154-norm4-production/raw"
CPP = ROOT / "src/norm4_source_thermal_replay.cpp"
BACKENDS = {
    "primitive": "src/threshold_rank_orientation_mc.cpp",
    "integer_period": "src/threshold_rank_integer_period_mc.cpp",
}
DESIGNS = {
    65: ([8, 1], [7, 4]),
    85: ([9, 2], [7, 6]),
    130: ([11, 3], [9, 7]),
    170: ([13, 1], [11, 7]),
    260: ([16, 2], [14, 8]),
    340: ([18, 4], [14, 12]),
}
REPLAY_SAMPLES = 100000
REPLAY_BATCHES = 100


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git_bytes(commit: str, path: str) -> bytes:
    return subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def display_path(path: Path) -> str:
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def source_run(n: int) -> dict:
    """Load the original production receipt, not a newly chosen RNG contract."""
    integer = n in (260, 340)
    seed = 2026105401 if n == 260 else 2026105402 if n == 340 else 2026104501
    first_counter = 8200000000 if integer else 5100000000
    original_samples = 1000000000 if integer else 1900000000
    suffix = "1b" if integer else "1900m"
    path = f"{METADATA_DIRECTORY}/n{n}_{suffix}.metadata.json"
    raw = git_bytes(METADATA_COMMIT, path)
    metadata = json.loads(raw)
    expected = {
        "git_commit": SOURCE_COMMIT,
        "seed": seed,
        "replica_counter_first": first_counter,
        "replica_counter_last_exclusive": first_counter + original_samples,
        "samples_per_pair": original_samples,
        "batches": 100,
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(f"N{n}: immutable production metadata differs at {key}")
    designs = [design for design in metadata["designs"] if design["N"] == n]
    if len(designs) != 1:
        raise ValueError(f"N{n}: original production geometry is not unique")
    design = designs[0]
    if (design["first"], design["second"]) != DESIGNS[n]:
        raise ValueError(f"N{n}: original lineage labels differ")
    return {
        "N": n,
        "backend_kind": "integer_period" if integer else "primitive",
        "seed": seed,
        "first": design["first"],
        "second": design["second"],
        "geometry": design,
        "old_counter_interval": [first_counter, first_counter + original_samples],
        "replayed_counter_interval": [first_counter, first_counter + REPLAY_SAMPLES],
        "old_samples_per_pair": original_samples,
        "reobserved_permutations": REPLAY_SAMPLES,
        "batches": REPLAY_BATCHES,
        "permutations_per_batch": REPLAY_SAMPLES // REPLAY_BATCHES,
        "metadata_source": {
            "commit": METADATA_COMMIT, "path": path, "sha256": sha(raw), "content": metadata,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/norm4-source-thermal")
    parser.add_argument("--workers", type=int, default=2, help="concurrent old-counter reobservations (1..6)")
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 6:
        raise ValueError("workers must be between 1 and 6")
    output_dir = args.output_dir.resolve()
    outputs = {n: output_dir / "raw" / f"n{n}.csv" for n in DESIGNS}
    receipt = output_dir / "run.json"
    if receipt.exists() or any(path.exists() for path in outputs.values()):
        raise ValueError("existing replay output: refusing to repeat or overwrite old-counter analysis")
    manifest_bytes = MANIFEST.read_bytes()
    contract = json.loads(manifest_bytes)
    manifest_runs = contract.get("runs", [])
    if sorted(run["N"] for run in manifest_runs) != sorted(DESIGNS):
        raise ValueError("manifest must contain the six original norm-4 production sizes exactly once")
    if contract["source_archive_commit"] != METADATA_COMMIT:
        raise ValueError("manifest source archive differs from the frozen production receipts")
    selection = contract["selection"]
    if (selection["selected_samples_per_N"], selection["selected_batches_per_N"],
            selection["selected_samples_per_batch"]) != (REPLAY_SAMPLES, REPLAY_BATCHES, 1000):
        raise ValueError("manifest must retain the declared old-counter 100k/100-batch subset")
    runs = [source_run(n) for n in DESIGNS]
    for run in runs:
        declared = next(item for item in manifest_runs if item["N"] == run["N"])
        for key in ("first", "second", "seed"):
            if key in declared and declared[key] != run[key]:
                raise ValueError(f"N{run['N']}: manifest differs from immutable metadata at {key}")
        expected_declared = {
            "selected_counter_interval": run["replayed_counter_interval"],
            "original_counter_interval": run["old_counter_interval"],
            "metadata_path": run["metadata_source"]["path"],
            "metadata_sha256_from_source_archive": run["metadata_source"]["sha256"],
        }
        for key, expected in expected_declared.items():
            if declared[key] != expected:
                raise ValueError(f"N{run['N']}: manifest differs from immutable input at {key}")
        for key in ("first_period_matrix", "second_period_matrix", "first_HNF", "second_HNF"):
            if key in declared and declared[key] != run["geometry"].get(key):
                raise ValueError(f"N{run['N']}: immutable quotient convention differs at {key}")
    backend_bytes = {kind: git_bytes(SOURCE_COMMIT, path) for kind, path in BACKENDS.items()}
    backends = {}
    for kind, path in BACKENDS.items():
        blob = subprocess.check_output(
            ["git", "rev-parse", f"{SOURCE_COMMIT}:{path}"], cwd=ROOT, text=True
        ).strip()
        backends[kind] = {"commit": SOURCE_COMMIT, "path": path, "git_blob": blob,
                          "sha256": sha(backend_bytes[kind])}
        profile = contract["engine_profiles"]["cyclic" if kind == "primitive" else kind]
        if (profile["engine_commit"], profile["engine_path"], profile["engine_blob"]) != (SOURCE_COMMIT, path, blob):
            raise ValueError(f"{kind}: manifest backend pin differs from the immutable compiled input")
    record = {
        "schema": "matching-one.norm4-source-thermal-run.v1",
        "status": "running",
        "started_utc": utc_now(),
        "execution_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "manifest": contract,
        "source_backends": backends,
        "source_runs": runs,
        "code": [{"path": display_path(path), "sha256": sha(path.read_bytes())}
                 for path in (Path(__file__).resolve(), CPP, MANIFEST)],
        "environment": {
            "python": platform.python_version(), "machine": platform.machine(),
            "platform": platform.platform(),
            "compiler": subprocess.check_output(["c++", "--version"], text=True).splitlines()[0],
        },
        "workers": args.workers,
        "new_samples": 0,
        "old_permutation_reobservations": 6 * REPLAY_SAMPLES,
        "configurations_per_permutation": "all K=0..N in each same-N direction; not independent samples",
        "sample_boundary": "first 100k old permutations per N, not original 1.9B/1B precision",
        "batch_boundary": "new 1000-permutation analysis blocks partition the original contiguous first 100k; not the old full-production 19M/10M batches",
        "dependency": {
            "shared_four_N": [65, 85, 130, 170],
            "shared_four_N_rule": "same seed/counters; original permutation key does not include N, so preserve joint aligned delete-one covariance",
            "same_N": "first and second share the identical original label permutation",
            "children": "N260 and N340 retain their distinct original seeds; independent domains under ordinary PRNG assumptions",
        },
        "readouts": {
            "q": "-1 + I[K>=K_minus] + I[K>=K_plus], using the original backend cross-wrapping definition",
            "E": "q*q",
            "s": "occupied black NN component count + complementary white matching (NN plus both diagonals) component count; unscaled integer",
            "downstream": "exact Binomial(N,p) conditional-K summation, not fixed-p importance weighting",
        },
        "server_actions": 0,
        "test_suites": [],
        "compile": [],
        "runs": [],
    }
    output_dir.joinpath("raw").mkdir(parents=True, exist_ok=True)
    with receipt.open("x") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")

    def save_receipt() -> None:
        receipt.write_text(json.dumps(record, indent=2) + "\n")

    started = time.perf_counter()
    try:
        with tempfile.TemporaryDirectory(prefix="matching-norm4-source-thermal-") as directory:
            build = Path(directory)
            binaries = {}
            for kind, raw in backend_bytes.items():
                archived = build / f"archived_{kind}.cpp"
                archived.write_bytes(raw)
                binary = build / f"norm4-source-thermal-{kind}"
                command = ["c++", "-O3", "-std=c++17", f'-DMATCHING_NORM4_BACKEND="{archived}"']
                if kind == "integer_period":
                    command.append("-DMATCHING_NORM4_INTEGER=1")
                command.extend([str(CPP), "-o", str(binary)])
                begin = time.perf_counter()
                process = subprocess.run(command, check=True, cwd=ROOT, text=True, capture_output=True)
                record["compile"].append({
                    "backend_kind": kind, "command": command,
                    "elapsed_seconds": time.perf_counter() - begin,
                    "binary_sha256": sha(binary.read_bytes()),
                    "stdout": process.stdout, "stderr": process.stderr,
                })
                binaries[kind] = binary
            save_receipt()

            def replay(run: dict) -> dict:
                n = run["N"]
                output = outputs[n]
                command = [str(binaries[run["backend_kind"]]), str(n), str(output)]
                begin = time.perf_counter()
                process = subprocess.run(command, check=True, cwd=ROOT, text=True, capture_output=True)
                return {
                    "N": n, "command": command, "elapsed_seconds": time.perf_counter() - begin,
                    "stdout": process.stdout, "stderr": process.stderr,
                    "output": display_path(output), "sha256": sha(output.read_bytes()),
                    "seed": run["seed"], "replayed_counter_interval": run["replayed_counter_interval"],
                    "reobserved_permutations": REPLAY_SAMPLES, "new_samples": 0,
                }

            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                futures = [executor.submit(replay, run) for run in runs]
                for future in as_completed(futures):
                    completed = future.result()
                    record["runs"].append(completed)
                    record["runs"].sort(key=lambda run: run["N"])
                    save_receipt()
                    print(completed["stdout"].strip(), flush=True)
        record["status"] = "completed"
        record["completed_utc"] = utc_now()
        record["elapsed_seconds"] = time.perf_counter() - started
        save_receipt()
    except Exception as error:
        record["status"] = "failed"
        record["failure"] = {"type": type(error).__name__, "message": str(error)}
        if isinstance(error, subprocess.CalledProcessError):
            record["failure"].update(stdout=error.stdout, stderr=error.stderr, returncode=error.returncode)
        record["elapsed_seconds"] = time.perf_counter() - started
        save_receipt()
        raise
    print(json.dumps({"status": record["status"], "runs": record["runs"], "new_samples": 0}, indent=2))


if __name__ == "__main__":
    main()
