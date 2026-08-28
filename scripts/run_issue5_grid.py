#!/usr/bin/env python3
"""Run the Issue #5 finite-size audit grid as independent processes.

Does not change the mathematical implementation of finite_size_audit.py.
Each job compares models 4 / 4,6 / 4,6,8 / 4,6,8,10 / 4,6,8,10,12.
"""

from __future__ import annotations

import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "scripts" / "finite_size_audit.py"
CSV = ROOT / "data" / "jacobsen_2015_square_site_cylinder.csv"
RAW = ROOT / "results" / "issue-5" / "raw"
LOGS = ROOT / "results" / "issue-5" / "logs"
MODELS = ["4", "4,6", "4,6,8", "4,6,8,10", "4,6,8,10,12"]
DPS_VALUES = (60, 100, 160)
MIN_TRAIN_VALUES = (5, 6, 7, 8, 9, 10)
HOLDOUT_VALUES = (2, 3, 4)
MAX_WORKERS = 8


def job_stem(dps: int, min_train: int, holdout: int) -> str:
    return f"audit_dps{dps:03d}_nmin{min_train:02d}_h{holdout:02d}"


def run_job(spec: tuple[int, int, int]) -> dict:
    dps, min_train, holdout = spec
    stem = job_stem(dps, min_train, holdout)
    json_path = RAW / f"{stem}.json"
    stdout_path = LOGS / f"{stem}.stdout.txt"
    stderr_path = LOGS / f"{stem}.stderr.txt"
    command = [
        sys.executable,
        str(AUDIT),
        str(CSV),
        "--models",
        *MODELS,
        "--min-train",
        str(min_train),
        "--holdout",
        str(holdout),
        "--dps",
        str(dps),
        "--json",
        str(json_path),
    ]
    with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr_handle:
        completed = subprocess.run(
            command,
            cwd=ROOT,
            stdout=stdout_handle,
            stderr=stderr_handle,
            check=False,
        )
    skipped: list[str] = []
    if stdout_path.exists():
        for line in stdout_path.read_text(encoding="utf-8").splitlines():
            if line.startswith("skip model"):
                skipped.append(line)
    parse_ok = False
    n_summaries = 0
    n_folds = 0
    if json_path.exists() and completed.returncode == 0:
        try:
            payload = json.loads(json_path.read_text(encoding="utf-8"))
            n_summaries = len(payload.get("summaries", []))
            n_folds = len(payload.get("folds", []))
            parse_ok = True
        except json.JSONDecodeError:
            parse_ok = False
    ok = completed.returncode == 0 and parse_ok
    return {
        "dps": dps,
        "min_train": min_train,
        "holdout": holdout,
        "stem": stem,
        "returncode": completed.returncode,
        "json": str(json_path.relative_to(ROOT)),
        "stdout": str(stdout_path.relative_to(ROOT)),
        "stderr": str(stderr_path.relative_to(ROOT)),
        "ok": ok,
        "parse_ok": parse_ok,
        "n_summaries": n_summaries,
        "n_folds": n_folds,
        "skipped": skipped,
        "command": command,
    }


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    jobs = [
        (dps, min_train, holdout)
        for dps in DPS_VALUES
        for min_train in MIN_TRAIN_VALUES
        for holdout in HOLDOUT_VALUES
    ]
    results: list[dict] = []
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as pool:
        futures = {pool.submit(run_job, spec): spec for spec in jobs}
        for future in as_completed(futures):
            record = future.result()
            results.append(record)
            status = "ok" if record["ok"] else "FAIL"
            print(
                f"{status:4} {record['stem']} "
                f"rc={record['returncode']} summaries={record['n_summaries']} "
                f"folds={record['n_folds']} skipped={len(record['skipped'])}",
                flush=True,
            )
    results.sort(key=lambda item: (item["dps"], item["min_train"], item["holdout"]))
    n_ok = sum(1 for item in results if item["ok"])
    n_fail = len(results) - n_ok
    manifest = {
        "n_jobs": len(results),
        "n_ok": n_ok,
        "n_fail": n_fail,
        "max_workers": MAX_WORKERS,
        "models": MODELS,
        "jobs": results,
    }
    manifest_path = ROOT / "results" / "issue-5" / "grid_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"grid complete: {n_ok} ok, {n_fail} fail, {len(results)} total")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
