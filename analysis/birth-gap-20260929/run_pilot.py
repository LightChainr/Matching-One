#!/usr/bin/env python3
"""Compile, exact-control, benchmark, then bounded independent joint-birth pilot.

Only Python stdlib and a C++17 compiler are required. No cloud or package install.
Each batch is a fresh mt19937_64 stream; lattice/size/batch streams are distinct.
"""
import argparse
import concurrent.futures
import datetime
import gzip
import hashlib
import json
import os
import pathlib
import platform
import re
import shlex
import subprocess
import tempfile
import time

HERE = pathlib.Path(__file__).resolve().parent


def seed_for(lattice, length, batch, role):
    tag = f"matching-one/birth-gap/20260929/{role}/{lattice}/{length}/{batch}"
    return int.from_bytes(hashlib.sha256(tag.encode()).digest()[:8], "big")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=int, nargs="+", default=[16, 32, 64, 128])
    ap.add_argument("--batches", type=int, default=8)
    ap.add_argument("--samples", type=int, default=2000)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--compiler", default="clang++", help="C++17 compiler command, e.g. g++")
    ap.add_argument("--role", default="pilot", help="production seed namespace; use a NEW value for independent blocks")
    ap.add_argument("--skip-exact-control", action="store_true")
    ap.add_argument("--exact-control-reference", help="prior exact-run reference required when skipping local exhaustive control")
    ap.add_argument("--git-base-commit", help="record source base when remote directory is not a Git checkout")
    ap.add_argument("--budget-seconds", type=float, default=240,
                    help="stop before production if benchmark-estimated wall time exceeds this bound")
    ap.add_argument("--output", type=pathlib.Path, default=HERE / "data")
    args = ap.parse_args()
    if args.workers > 14 or args.workers < 1:
        ap.error("bounded pilot permits 1..14 workers; choose from measured local/remote quota")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", args.role):
        ap.error("--role must be a filename-safe nonempty seed namespace")
    if args.skip_exact_control and not args.exact_control_reference:
        ap.error("--skip-exact-control requires --exact-control-reference")
    if min(args.sizes) < 3 or max(args.sizes) > 1024 or args.batches < 2 or args.samples < 1 or args.budget_seconds <= 0:
        ap.error("require sizes 3..1024, >=2 batches, positive samples and budget")
    args.output.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output / "run.json"
    if manifest_path.exists():
        raise SystemExit("Refusing to overwrite an existing run; choose another --output")
    compiler = shlex.split(args.compiler)
    if not compiler:
        ap.error("empty compiler command")
    base_commit = args.git_base_commit
    if base_commit is None:
        try:
            base_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE,
                                                 text=True, stderr=subprocess.DEVNULL).strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            base_commit = "unavailable-not-a-git-checkout"
    run = {"schema": "matching-one.joint-birth-run.v1",
           "utc_started": datetime.datetime.now(datetime.timezone.utc).isoformat(),
           "command": shlex.join(["python3", str(pathlib.Path(__file__).resolve()), *os.sys.argv[1:]]),
           "platform": platform.platform(), "machine": platform.machine(),
           "python": platform.python_version(), "workers": args.workers,
           "compiler": subprocess.check_output([*compiler, "--version"], text=True).strip(),
           "git_base_commit": base_commit,
           "production_role": args.role,
           "exact_control_reference": args.exact_control_reference,
           "exact_control_skipped": args.skip_exact_control,
           "estimated_wall_budget_seconds": args.budget_seconds,
           "source_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in (HERE / "engine.cpp", pathlib.Path(__file__).resolve())},
           "seed_rule": "first 64 SHA256 bits of matching-one/birth-gap/20260929/{role}/{lattice}/{L}/{batch}",
           "random_streams": "separate mt19937_64 per (role,lattice,L,batch); no paired/cross-size common random numbers",
           "sampling": "unbiased Fisher-Yates of all N sites; stop after same-filtration global rank reaches 2",
           "controls": [], "benchmarks": [], "batches": []}
    start = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="matching-birth-gap-build-") as tmp:
        binary = pathlib.Path(tmp) / "birth-gap"
        build = [*compiler, "-O3", "-std=c++17", "-DNDEBUG", str(HERE / "engine.cpp"), "-o", str(binary)]
        subprocess.run(build, check=True)
        run["build_command"] = shlex.join(build)

        def job(lattice, length, batch, role, n, mode="sample"):
            name = f"{role}-{lattice}-L{length}-b{batch:02d}.json"
            raw = pathlib.Path(tmp) / name
            seed = seed_for(lattice, length, batch, role)
            command = [str(binary), lattice, str(length), str(n), str(seed), str(raw), mode]
            log = subprocess.check_output(command, text=True).strip()
            content = raw.read_bytes()
            result = json.loads(content)
            target = args.output / (name + ".gz")
            with target.open("wb") as f:
                with gzip.GzipFile(fileobj=f, mode="wb", mtime=0) as gz:
                    gz.write(content)
            print(log, flush=True)
            return {"file": target.name, "lattice": lattice, "L": length, "batch": batch,
                    "role": role, "seed": str(seed), "samples": result["samples"],
                    "runtime_seconds": result["runtime_seconds"], "command": shlex.join(command),
                    "sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
                    "dependency_group": f"{role}/{lattice}/L{length}/batch{batch}"}

        if not args.skip_exact_control:
            for lattice in ("square", "triangular"):
                run["controls"].append(job(lattice, 3, 0, "exact-control", 0, "exact"))
        # Benchmark the largest requested cell sequentially BEFORE production.
        # These samples have separate seeds and are excluded from every estimate.
        for lattice in ("square", "triangular"):
            benchmark_role = "benchmark" if args.role == "pilot" else f"benchmark-{args.role}"
            run["benchmarks"].append(job(lattice, max(args.sizes), 0, benchmark_role, 200))
        estimated_cpu = sum(b["runtime_seconds"] / b["samples"] for b in run["benchmarks"]) * args.samples * args.batches * sum((L / max(args.sizes)) ** 2 for L in args.sizes)
        run["estimated_production_cpu_seconds"] = estimated_cpu
        print(f"largest-cell benchmark predicts roughly {estimated_cpu:.1f} aggregate CPU-seconds", flush=True)
        if estimated_cpu / args.workers > args.budget_seconds:
            run["status"] = "benchmark-budget-stop"
            manifest_path.write_text(json.dumps(run, indent=2) + "\n")
            raise SystemExit("Benchmark exceeds bounded pilot; no production started")
        jobs = [(lattice, L, b) for lattice in ("square", "triangular")
                for L in args.sizes for b in range(args.batches)]
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(job, lattice, L, b, args.role, args.samples) for lattice, L, b in jobs]
            for future in concurrent.futures.as_completed(futures):
                run["batches"].append(future.result())
                run["batches"].sort(key=lambda x: (x["lattice"], x["L"], x["batch"]))
                manifest_path.write_text(json.dumps(run, indent=2) + "\n")
    run["wall_seconds"] = time.monotonic() - start
    run["aggregate_engine_seconds"] = sum(b["runtime_seconds"] for b in run["batches"])
    run["status"] = "completed"
    run["utc_finished"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    manifest_path.write_text(json.dumps(run, indent=2) + "\n")
    print(f"Completed in {run['wall_seconds']:.2f}s; manifest {manifest_path}")


if __name__ == "__main__":
    main()
