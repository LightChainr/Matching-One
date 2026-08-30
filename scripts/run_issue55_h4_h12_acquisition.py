#!/usr/bin/env python3
"""Strict launcher for the frozen Issue #55 H4/H12 acquisition.

The launcher reads the acquisition manifest, revalidates its design and RNG
domains, and emits or executes exactly one smoke campaign or one production
shard.  It never chooses geometries or a sample count from observed means.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from typing import Any, Dict, List, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "experiments" / "issue55_h4_h12_orthogonal_acquisition_20260830.json"
COMMIT_RE = re.compile(r"^[0-9a-fA-F]{40}$")

sys.path.insert(0, str(ROOT / "scripts"))
from rng_domain_policy import derive_size_seed  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def load_and_validate(path: Path) -> Dict[str, Any]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    _require(
        manifest.get("schema") == "matching-one/issue55-h4-h12-orthogonal-acquisition/v1",
        "unknown acquisition schema",
    )
    frozen = manifest["frozen_design"]
    design_path = ROOT / frozen["path"]
    _require(sha256(design_path) == frozen["sha256"], "frozen design manifest hash drift")
    engine = manifest["engine"]
    _require(
        sha256(ROOT / engine["path"]) == engine["sha256_at_source_commit"],
        "threshold-rank engine hash drift; freeze a reviewed engine extension instead",
    )
    source = json.loads(design_path.read_text(encoding="utf-8"))
    source_designs = {int(item["N"]): item for item in source["designs"]}
    designs = manifest["designs"]
    _require([int(item["N"]) for item in designs] == [305, 325], "design order must be N305,N325")
    for item in designs:
        n = int(item["N"])
        _require(n in source_designs, "acquisition contains an unregistered design")
        original = source_designs[n]
        for field in ("id", "first", "second", "delta_cos4", "delta_cos12"):
            _require(item[field] == original[field], "{} N{} drift".format(field, n))
        _require(item["alias_ratio"] == original["alias_ratio_cos12_over_cos4"], "alias ratio drift")
        a, b = item["first"]
        c, d = item["second"]
        _require(item["first_period_matrix"] == [[a, -b], [b, a]], "first matrix drift")
        _require(item["second_period_matrix"] == [[c, -d], [d, c]], "second matrix drift")
    rng = manifest["rng"]
    for stage in ("smoke", "production"):
        domain = rng[stage]
        for n in (305, 325):
            expected = derive_size_seed(int(rng["base_seed"]), domain["experiment_tag"], n)
            _require(int(domain["effective_seed_by_N"][str(n)]) == expected, "RNG derivation drift")
    smoke = manifest["variance_only_smoke"]
    _require(int(smoke["samples_per_design"]) == 20000, "smoke sample count drift")
    _require(int(smoke["batches"]) == 20, "smoke batch count drift")
    _require(
        int(smoke["counter_last_exclusive"]) - int(smoke["counter_first"])
        == int(smoke["samples_per_design"]),
        "smoke counter interval drift",
    )
    pilot = smoke.get("frozen_result", {})
    _require(
        sha256(ROOT / pilot["path"]) == pilot["sha256"],
        "variance-only pilot artifact hash drift",
    )
    pilot_payload = json.loads((ROOT / pilot["path"]).read_text(encoding="utf-8"))
    _require(
        pilot_payload.get("status") == "variance_only_target_means_withheld",
        "pilot artifact is not target-blind",
    )
    _require(
        pilot_payload.get("recommended_samples_per_design") == pilot["recommended_samples_per_design"],
        "pilot recommendation drift",
    )
    production = manifest["production"]
    _require(int(production["shard_count"]) == 3, "three production shards required")
    _require(int(production["batches_per_shard"]) >= 20, "too few production batches")
    total = production.get("samples_per_design")
    _require(total == pilot["recommended_samples_per_design"], "production count differs from blind freeze")
    _require(total % int(production["shard_count"]) == 0, "production count not shardable")
    _require(
        int(production["counter_last_exclusive"]) - int(production["counter_first"]) == total,
        "production counter interval drift",
    )
    return manifest


def _flatten_matrix(matrix: Sequence[Sequence[int]]) -> List[str]:
    return [str(value) for row in matrix for value in row]


def commands(
    manifest: Mapping[str, Any], stage: str, shard_index: int, binary: Path,
    output_dir: Path, threads: int, git_commit: str,
) -> List[Dict[str, Any]]:
    _require(COMMIT_RE.fullmatch(git_commit) is not None, "--git-commit must be full 40-hex")
    if stage == "smoke":
        _require(shard_index == 0, "smoke has only shard index 0")
        run = manifest["variance_only_smoke"]
        samples = int(run["samples_per_design"])
        batches = int(run["batches"])
        counter_first = int(run["counter_first"])
    else:
        run = manifest["production"]
        total = run.get("samples_per_design")
        _require(isinstance(total, int) and total > 0, "production count is not frozen")
        shard_count = int(run["shard_count"])
        _require(0 <= shard_index < shard_count, "production shard index out of range")
        _require(total % shard_count == 0, "production count is not divisible by shard count")
        samples = total // shard_count
        batches = int(run["batches_per_shard"])
        _require(samples % batches == 0, "shard samples not divisible by batches")
        counter_first = int(run["counter_first"]) + shard_index * samples
        _require(
            int(run["counter_last_exclusive"]) == int(run["counter_first"]) + total,
            "production counter interval drift",
        )
    rng = manifest["rng"][stage]
    output: List[Dict[str, Any]] = []
    for design in manifest["designs"]:
        n = int(design["N"])
        prefix = output_dir / "issue55_{}_n{}_shard{}".format(stage, n, shard_index)
        argv = [
            str(binary),
            "--samples", str(samples),
            "--batches", str(batches),
            "--seed", str(rng["effective_seed_by_N"][str(n)]),
            "--replica-offset", str(counter_first),
            "--threads", str(threads),
            "--first-matrix", *_flatten_matrix(design["first_period_matrix"]),
            "--second-matrix", *_flatten_matrix(design["second_period_matrix"]),
            "--first-rep", *map(str, design["first"]),
            "--second-rep", *map(str, design["second"]),
            "--git-commit", git_commit,
            "--output-prefix", str(prefix),
        ]
        output.append({
            "N": n,
            "stage": stage,
            "shard_index": shard_index,
            "samples": samples,
            "batches": batches,
            "seed": int(rng["effective_seed_by_N"][str(n)]),
            "counter_interval": [counter_first, counter_first + samples],
            "output_prefix": str(prefix),
            "argv": argv,
        })
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--stage", choices=("smoke", "production"), required=True)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--git-commit", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--parallel-sizes", action="store_true")
    args = parser.parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    manifest = load_and_validate(args.manifest)
    planned = commands(
        manifest, args.stage, args.shard_index, args.binary, args.output_dir,
        args.threads, args.git_commit,
    )
    print(json.dumps(planned, indent=2))
    if args.dry_run:
        return 0
    if not args.binary.is_file():
        raise ValueError("binary does not exist: {}".format(args.binary))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for item in planned:
        for suffix in (".hist.csv", ".moments.csv", ".metadata.json"):
            if Path(item["output_prefix"] + suffix).exists():
                raise ValueError("refusing to overwrite {}".format(item["output_prefix"] + suffix))
    if args.parallel_sizes:
        processes = [subprocess.Popen(item["argv"]) for item in planned]
        statuses = [process.wait() for process in processes]
        if any(status != 0 for status in statuses):
            raise subprocess.CalledProcessError(max(statuses), "parallel size acquisition")
    else:
        for item in planned:
            subprocess.run(item["argv"], check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
