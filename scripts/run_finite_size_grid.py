#!/usr/bin/env python3
"""Run the finite-size audit parameter grid with bounded concurrency.

The default grid contains 54 independent jobs (3 precisions x 6 starting
widths x 3 holdout lengths).  Existing, valid JSON outputs are reused, so a
server run can be resumed after interruption.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence


DEFAULT_MODELS = ("4", "4,6", "4,6,8", "4,6,8,10", "4,6,8,10,12")


@dataclass(frozen=True)
class Job:
    dps: int
    min_train: int
    holdout: int

    @property
    def name(self) -> str:
        return f"dps-{self.dps}_nmin-{self.min_train}_holdout-{self.holdout}"


@dataclass(frozen=True)
class JobResult:
    name: str
    status: str
    seconds: float
    output: str
    stdout_log: str
    stderr_log: str
    returncode: int
    error: str = ""


def build_jobs(
    dps_values: Iterable[int],
    min_train_values: Iterable[int],
    holdout_values: Iterable[int],
) -> list[Job]:
    return [
        Job(dps=dps, min_train=n_min, holdout=holdout)
        for dps, n_min, holdout in itertools.product(
            dps_values, min_train_values, holdout_values
        )
    ]


def valid_output(
    path: Path,
    job: Job,
    *,
    csv_path: Optional[Path] = None,
    models: Optional[Sequence[str]] = None,
) -> bool:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if (
            payload["dps"] != job.dps
            or payload["min_train"] != job.min_train
            or payload["holdout"] != job.holdout
            or not payload["summaries"]
        ):
            return False
        if csv_path is not None and payload.get("input") != str(csv_path):
            return False
        if models is not None and {
            str(summary["model"]) for summary in payload["summaries"]
        } != set(models):
            return False
        return all(
            int(fold["train_max"]) < int(fold["test_min"])
            for fold in payload["folds"]
        )
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return False


def command_for(
    job: Job,
    python: str,
    audit_script: Path,
    csv_path: Path,
    models: Sequence[str],
    output: Path,
) -> list[str]:
    return [
        python,
        str(audit_script),
        str(csv_path),
        "--models",
        *models,
        "--min-train",
        str(job.min_train),
        "--holdout",
        str(job.holdout),
        "--dps",
        str(job.dps),
        "--json",
        str(output),
    ]


def run_job(
    job: Job,
    *,
    python: str,
    audit_script: Path,
    csv_path: Path,
    models: Sequence[str],
    output_dir: Path,
    force: bool,
) -> JobResult:
    output = output_dir / "raw" / f"{job.name}.json"
    stdout_log = output_dir / "logs" / f"{job.name}.stdout.txt"
    stderr_log = output_dir / "logs" / f"{job.name}.stderr.txt"
    if not force and valid_output(output, job, csv_path=csv_path, models=models):
        return JobResult(
            job.name,
            "reused",
            0.0,
            str(output),
            str(stdout_log),
            str(stderr_log),
            0,
        )

    started = time.monotonic()
    command = command_for(job, python, audit_script, csv_path, models, output)
    with stdout_log.open("w", encoding="utf-8") as stdout_handle, stderr_log.open(
        "w", encoding="utf-8"
    ) as stderr_handle:
        completed = subprocess.run(
            command,
            stdout=stdout_handle,
            stderr=stderr_handle,
            check=False,
        )
    elapsed = time.monotonic() - started
    if completed.returncode != 0:
        return JobResult(
            job.name,
            "failed",
            elapsed,
            str(output),
            str(stdout_log),
            str(stderr_log),
            completed.returncode,
            "finite_size_audit exited nonzero",
        )
    if not valid_output(output, job, csv_path=csv_path, models=models):
        return JobResult(
            job.name,
            "failed",
            elapsed,
            str(output),
            str(stdout_log),
            str(stderr_log),
            0,
            "audit output is missing or invalid",
        )
    return JobResult(
        job.name,
        "completed",
        elapsed,
        str(output),
        str(stdout_log),
        str(stderr_log),
        0,
    )


def write_manifest(
    path: Path,
    *,
    csv_path: Path,
    models: Sequence[str],
    workers: int,
    jobs: Sequence[Job],
    results: Sequence[JobResult],
) -> None:
    payload = {
        "schema_version": 1,
        "input": str(csv_path),
        "models": list(models),
        "workers": workers,
        "jobs": [asdict(job) for job in jobs],
        "results": [asdict(result) for result in sorted(results, key=lambda row: row.name)],
    }
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def parse_int_list(text: str) -> tuple[int, ...]:
    try:
        values = tuple(int(part) for part in text.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from exc
    if not values or any(value <= 0 for value in values):
        raise argparse.ArgumentTypeError("values must be positive integers")
    return values


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "csv",
        type=Path,
        nargs="?",
        default=root / "data/jacobsen_2015_square_site_cylinder.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=root / "results/issue-5/grid")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--audit-script", type=Path, default=root / "scripts/finite_size_audit.py")
    parser.add_argument("--models", nargs="+", default=list(DEFAULT_MODELS))
    parser.add_argument("--dps", type=parse_int_list, default=(60, 100, 160))
    parser.add_argument("--min-train", type=parse_int_list, default=(5, 6, 7, 8, 9, 10))
    parser.add_argument("--holdout", type=parse_int_list, default=(2, 3, 4))
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--force", action="store_true", help="rerun valid existing jobs")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.workers <= 0:
        raise SystemExit("workers must be positive")
    if any(dps < 40 for dps in args.dps):
        raise SystemExit("all dps values must be at least 40")
    if not args.csv.is_file():
        raise SystemExit(f"input CSV not found: {args.csv}")
    if not args.audit_script.is_file():
        raise SystemExit(f"audit script not found: {args.audit_script}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "raw").mkdir(exist_ok=True)
    (args.output_dir / "logs").mkdir(exist_ok=True)
    jobs = build_jobs(args.dps, args.min_train, args.holdout)
    results: list[JobResult] = []
    manifest = args.output_dir / "manifest.json"

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                run_job,
                job,
                python=args.python,
                audit_script=args.audit_script,
                csv_path=args.csv,
                models=args.models,
                output_dir=args.output_dir,
                force=args.force,
            ): job
            for job in jobs
        }
        for future in as_completed(futures):
            job = futures[future]
            try:
                result = future.result()
            except Exception as exc:  # preserve other jobs and record unexpected failures
                result = JobResult(job.name, "failed", 0.0, "", "", "", -1, repr(exc))
            results.append(result)
            write_manifest(
                manifest,
                csv_path=args.csv,
                models=args.models,
                workers=args.workers,
                jobs=jobs,
                results=results,
            )
            print(f"{result.status:9} {result.name} ({result.seconds:.2f}s)", flush=True)

    failures = [result for result in results if result.status == "failed"]
    print(f"manifest: {manifest}")
    print(f"jobs: {len(results)}, failures: {len(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
