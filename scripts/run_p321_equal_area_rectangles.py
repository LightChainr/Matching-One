#!/usr/bin/env python3
"""Run one equal-area P321 rectangle campaign with aligned CRN batches."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DESIGN = ROOT / "analysis/p321_equal_area_rectangle_design.json"


def _flatten(matrix: Sequence[Sequence[int]]) -> list[str]:
    return [str(value) for row in matrix for value in row]


def rows_for_n(design: Mapping[str, Any], n: int) -> tuple[Mapping[str, Any], list[Mapping[str, Any]]]:
    rows = [row for row in design["rows"] if int(row["N"]) == n]
    square = [row for row in rows if row["aspect_ratio"] == "1/1"]
    if len(square) != 1 or len(rows) != 5:
        raise ValueError(f"design must contain one square and four rectangles at N={n}")
    rectangles = [row for row in rows if row["aspect_ratio"] != "1/1"]
    return square[0], rectangles


def run_campaign(
    *,
    binary: Path,
    design_path: Path,
    n: int,
    samples: int,
    batches: int,
    seed: int,
    replica_offset: int,
    threads: int,
    parallel_runs: int,
    git_commit: str,
    output_dir: Path,
) -> dict[str, Any]:
    if samples <= 0 or batches < 2 or samples % batches:
        raise ValueError("samples must be positive and divisible by at least two batches")
    design = json.loads(design_path.read_text(encoding="utf-8"))
    if design.get("schema") != "matching-one.p321-equal-area-rectangle-design.v1":
        raise ValueError("wrong P321 geometry design")
    square, rectangles = rows_for_n(design, n)
    output_dir.mkdir(parents=True, exist_ok=True)
    if not 1 <= parallel_runs <= 4:
        raise ValueError("parallel_runs must lie in 1..4")
    runs = []
    for rectangle in rectangles:
        label = rectangle["aspect_ratio"].replace("/", "_")
        prefix = output_dir / f"N{n}_rho_{label}"
        command = [
            str(binary),
            "--samples", str(samples),
            "--batches", str(batches),
            "--seed", str(seed),
            "--replica-offset", str(replica_offset),
            "--threads", str(threads),
            "--git-commit", git_commit,
            "--first-matrix", *_flatten(square["period_matrix_row_major"]),
            "--second-matrix", *_flatten(rectangle["period_matrix_row_major"]),
            "--first-rep", str(square["width"]), "0",
            "--second-rep", str(rectangle["width"]), "0",
            "--output-prefix", str(prefix),
        ]
        runs.append({
            "rho": rectangle["aspect_ratio"],
            "role": rectangle["role"],
            "prefix": prefix.name,
            "histogram": prefix.name + ".hist.csv",
            "moments": prefix.name + ".moments.csv",
            "metadata": prefix.name + ".metadata.json",
            "command": command,
        })
    def execute(run: Mapping[str, Any]) -> None:
        subprocess.run(run["command"], cwd=ROOT, check=True)

    if parallel_runs == 1:
        for run in runs:
            execute(run)
    else:
        with ThreadPoolExecutor(max_workers=parallel_runs) as executor:
            list(executor.map(execute, runs))
    manifest = {
        "schema": "matching-one/p321-equal-area-campaign/v1",
        "status": "local_variance_smoke" if samples <= 100_000 else "production",
        "design": str(design_path.relative_to(ROOT)) if design_path.is_relative_to(ROOT) else str(design_path),
        "N": n,
        "samples_per_shape": samples,
        "batches": batches,
        "seed": seed,
        "replica_counter_first": replica_offset,
        "replica_counter_last_exclusive": replica_offset + samples,
        "git_commit": git_commit,
        "parallel_rectangle_processes": parallel_runs,
        "square_rerun_count": len(runs),
        "runs": runs,
    }
    manifest_path = output_dir / "campaign.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--design", type=Path, default=DEFAULT_DESIGN)
    parser.add_argument("--n", type=int, required=True, choices=(144, 576, 1296))
    parser.add_argument("--samples", type=int, default=20_000)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--seed", type=int, default=32114420260830)
    parser.add_argument("--replica-offset", type=int, default=0)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--parallel-runs", type=int, default=1, choices=(1, 2, 3, 4),
        help="run distinct rectangle subprocesses concurrently",
    )
    parser.add_argument("--git-commit", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    manifest = run_campaign(
        binary=args.binary,
        design_path=args.design,
        n=args.n,
        samples=args.samples,
        batches=args.batches,
        seed=args.seed,
        replica_offset=args.replica_offset,
        threads=args.threads,
        parallel_runs=args.parallel_runs,
        git_commit=args.git_commit,
        output_dir=args.output_dir,
    )
    print(json.dumps({"campaign": str(args.output_dir / "campaign.json"), "runs": len(manifest["runs"])}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
