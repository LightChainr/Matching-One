#!/usr/bin/env python3
"""Run the existing norm-5 analyses from one threshold-rank data block.

This is deliberately a thin orchestration layer.  It does not define a new
statistic and does not change any frozen scorer.  It infers exact shared-random
counter groups from metadata, runs each existing analysis independently, and
records failures without blocking the other views unless ``--fail-on-error``
is requested.

Expected run specification:

    --run N:HISTOGRAM:MOMENTS:METADATA

for N=65,85,130,170,325,425.  The same six full-curve histograms can then feed
all compatible fixed/frozen and exploratory analyses.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import subprocess
import sys
from typing import Mapping, Sequence


REQUIRED_SIZES = (65, 85, 130, 170, 325, 425)
PRIMARY_SIZES = (65, 85, 325, 425)
TARGET_SIZES = (325, 425)


@dataclass(frozen=True)
class RunSpec:
    n: int
    histogram: Path
    moments: Path
    metadata: Path


def parse_run(text: str) -> RunSpec:
    fields = text.split(":", 3)
    if len(fields) != 4:
        raise argparse.ArgumentTypeError("run must be N:HISTOGRAM:MOMENTS:METADATA")
    try:
        n = int(fields[0])
    except ValueError as exc:
        raise argparse.ArgumentTypeError("run N must be an integer") from exc
    return RunSpec(n, Path(fields[1]), Path(fields[2]), Path(fields[3]))


def load_metadata(run: RunSpec) -> Mapping[str, object]:
    payload = json.loads(run.metadata.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"N={run.n}: metadata must be a JSON object")
    return payload


def counter_group_key(metadata: Mapping[str, object]) -> tuple[int, int, int]:
    return (
        int(metadata["seed"]),
        int(metadata["replica_counter_first"]),
        int(metadata["replica_counter_last_exclusive"]),
    )


def infer_covariance_groups(
    runs: Sequence[RunSpec], metadata: Mapping[int, Mapping[str, object]]
) -> list[list[int]]:
    groups: dict[tuple[int, int, int], list[int]] = {}
    for run in runs:
        groups.setdefault(counter_group_key(metadata[run.n]), []).append(run.n)
    return [sorted(values) for _, values in sorted(groups.items(), key=lambda item: min(item[1]))]


def validate_runs(runs: Sequence[RunSpec], *, check_files: bool = True) -> dict[int, RunSpec]:
    by_n: dict[int, RunSpec] = {}
    for run in runs:
        if run.n in by_n:
            raise ValueError(f"duplicate N={run.n}")
        by_n[run.n] = run
        if check_files:
            for path in (run.histogram, run.moments, run.metadata):
                if not path.is_file():
                    raise ValueError(f"N={run.n}: missing input {path}")
    missing = [n for n in REQUIRED_SIZES if n not in by_n]
    if missing:
        raise ValueError(f"norm-5 bundle requires sizes {REQUIRED_SIZES}; missing {missing}")
    return by_n


def build_commands(
    root: Path,
    by_n: Mapping[int, RunSpec],
    covariance_groups: Sequence[Sequence[int]],
    output_dir: Path,
    source_rank_gap_score: Path,
) -> list[tuple[str, list[str], Path]]:
    python = sys.executable
    scripts = root / "scripts"
    output_dir = output_dir.resolve()

    primary_output = output_dir / "primary_harmonic.json"
    primary = [python, str(scripts / "score_norm5_harmonic_primary_typed.py")]
    for n in PRIMARY_SIZES:
        run = by_n[n]
        primary += ["--run", f"{n}:{run.histogram}:{run.metadata}"]
    primary += ["--output", str(primary_output)]

    cocycle_output = output_dir / "intrinsic_functional_cocycle.json"
    cocycle = [python, str(scripts / "score_intrinsic_functional_cocycle_typed.py")]
    cocycle += ["--histograms", *[str(by_n[n].histogram) for n in REQUIRED_SIZES]]
    cocycle += [
        "--covariance-groups",
        *[",".join(str(n) for n in group) for group in covariance_groups],
        "--json",
        str(cocycle_output),
    ]

    modes_output = output_dir / "krawtchouk_score_modes.json"
    modes = [python, str(scripts / "threshold_score_modes.py")]
    modes += [str(by_n[n].histogram) for n in REQUIRED_SIZES]
    modes += ["--output", str(modes_output)]

    thermal_jet_output = output_dir / "thermal_jet_score.json"
    thermal_jet = [python, str(scripts / "score_norm5_thermal_jet.py")]
    for n in REQUIRED_SIZES:
        run = by_n[n]
        thermal_jet += [
            "--run",
            str(n),
            str(run.histogram),
            str(run.moments),
            str(run.metadata),
        ]
    thermal_jet += [
        "--lineage", "65", "130", "325",
        "--lineage", "85", "170", "425",
        "--output", str(thermal_jet_output),
    ]

    rank_gap_output = output_dir / "rank_gap_boundary_score.json"
    rank_gap = [
        python,
        str(scripts / "score_rank_gap_boundary_targets.py"),
        "--source-score",
        str(source_rank_gap_score),
    ]
    for n in TARGET_SIZES:
        run = by_n[n]
        rank_gap += ["--target-run", f"{n}:{run.moments}:{run.metadata}"]
    rank_gap += ["--output", str(rank_gap_output)]

    return [
        ("primary_harmonic", primary, primary_output),
        ("intrinsic_functional_cocycle", cocycle, cocycle_output),
        ("krawtchouk_score_modes", modes, modes_output),
        ("thermal_jet_score", thermal_jet, thermal_jet_output),
        ("rank_gap_boundary", rank_gap, rank_gap_output),
    ]


def run_command(name: str, command: Sequence[str], output_path: Path, *, dry_run: bool) -> dict:
    record = {
        "name": name,
        "command": list(command),
        "output": str(output_path),
    }
    if dry_run:
        record["status"] = "DRY_RUN"
        return record
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    record.update(
        {
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "status": "OK" if completed.returncode == 0 else "FAILED_CONTINUING",
            "output_exists": output_path.is_file(),
        }
    )
    return record


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", type=parse_run, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--source-rank-gap-score",
        type=Path,
        default=root / "results/server-20260828/rank-gap-thermal-window/score.json",
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--fail-on-error",
        action="store_true",
        help="return nonzero if any subanalysis fails; default is to record and continue",
    )
    args = parser.parse_args()

    by_n = validate_runs(args.run, check_files=not args.dry_run)
    metadata = (
        {n: load_metadata(by_n[n]) for n in REQUIRED_SIZES}
        if not args.dry_run
        else {
            n: {
                "seed": n,
                "replica_counter_first": 0,
                "replica_counter_last_exclusive": 1,
            }
            for n in REQUIRED_SIZES
        }
    )
    groups = infer_covariance_groups([by_n[n] for n in REQUIRED_SIZES], metadata)
    commands = build_commands(
        root,
        by_n,
        groups,
        args.output_dir,
        args.source_rank_gap_score,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = [
        run_command(name, command, output, dry_run=args.dry_run)
        for name, command, output in commands
    ]
    failed = [row["name"] for row in records if row["status"] == "FAILED_CONTINUING"]
    payload = {
        "schema": "matching-one/norm5-analysis-bundle/v1",
        "policy": (
            "thin orchestration only; frozen scorers are unchanged; one failed or "
            "inapplicable frozen analysis does not block other exploratory views"
        ),
        "size_order": list(REQUIRED_SIZES),
        "inferred_covariance_groups": groups,
        "runs": records,
        "status": "PARTIAL" if failed else ("DRY_RUN" if args.dry_run else "OK"),
        "failed_analyses": failed,
    }
    summary = args.output_dir / "bundle.json"
    summary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(summary)
    return 1 if args.fail_on_error and failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
