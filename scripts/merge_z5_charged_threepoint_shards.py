#!/usr/bin/env python3
"""Merge disjoint frozen P250 replica shards without losing batch covariance."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from z5_charged_threepoint_mc import JOINT_REAL_ORDER, PRODUCTION_ID, summarize, write_batches


INTEGER_COLUMNS = ("batch", "replica_first", "samples")
FLOAT_COLUMNS = ("conjugacy_max_abs",) + JOINT_REAL_ORDER
DIGEST_COLUMNS = ("field_sha256", "translation_sha256")


def read_batches(path: Path) -> list[dict]:
    rows = []
    with path.open(newline="") as handle:
        for raw in csv.DictReader(handle):
            row = {name: int(raw[name]) for name in INTEGER_COLUMNS}
            row.update({name: float(raw[name]) for name in FLOAT_COLUMNS})
            row.update({name: raw[name] for name in DIGEST_COLUMNS})
            rows.append(row)
    return rows


def merge_payloads(payloads: list[dict], rows: list[dict], execution: dict) -> dict:
    if execution.get("production_id") != PRODUCTION_ID or not execution.get("production_authorized"):
        raise ValueError("execution manifest is not authorized")
    runner_commit = execution["runner_commit"]
    for payload in payloads:
        if payload.get("schema") != "matching-one/z5-charged-threepoint-response/v1":
            raise ValueError("unexpected shard response schema")
        if payload.get("production_id") != PRODUCTION_ID:
            raise ValueError("shard production id changed")
        if payload.get("manifest_runner_commit") != runner_commit:
            raise ValueError("shard runner commit changed")
        if not payload["mapping_gate"].get("passed"):
            raise ValueError("shard mapping gate failed")

    ordered = sorted(rows, key=lambda row: row["replica_first"])
    expected = execution["run"]["replica_offset"]
    for index, row in enumerate(ordered):
        if row["replica_first"] != expected:
            raise ValueError("replica shard coverage is not exact and contiguous")
        row["batch"] = index
        expected += row["samples"]
    if expected != execution["run"]["replica_offset"] + execution["run"]["samples"]:
        raise ValueError("replica shard coverage does not reach the frozen endpoint")
    if len(ordered) != execution["run"]["batches"]:
        raise ValueError("merged batch count changed")

    reference = payloads[0]
    for payload in payloads[1:]:
        if payload["mapping_gate"] != reference["mapping_gate"]:
            raise ValueError("mapping gates differ across shards")
        if payload["observable"] != reference["observable"]:
            raise ValueError("observables differ across shards")
        for key in ("p", "seed", "radius"):
            if payload["run"][key] != reference["run"][key]:
                raise ValueError(f"shard field contract differs for {key}")

    return {
        "schema": "matching-one/z5-charged-threepoint-response/v1",
        "status": "production_under_frozen_sharded_execution",
        "production_id": PRODUCTION_ID,
        "manifest_runner_commit": runner_commit,
        "issues": [250],
        "mapping_gate": reference["mapping_gate"],
        "run": {
            **execution["run"],
            "replica_last_exclusive": expected,
        },
        "observable": reference["observable"],
        "execution": {
            "schema": execution["schema"],
            "shards": [payload["run"] for payload in payloads],
            "merge_rule": "sort 10k batches by replica_first, verify exact disjoint coverage, then recompute all moments and covariance from the 100 batch sufficient statistics",
        },
        "analysis": summarize(ordered),
    }


def load_and_merge(execution_path: Path, input_dir: Path) -> tuple[dict, list[dict]]:
    execution = json.loads(execution_path.read_text())
    payloads = []
    rows = []
    for shard_path in execution["shard_manifests"]:
        shard = json.loads((execution_path.parents[1] / shard_path).read_text())
        response_path = input_dir / Path(shard["run"]["output"]).name
        batches_path = input_dir / Path(shard["run"]["batches_output"]).name
        payload = json.loads(response_path.read_text())
        expected = shard["run"]
        for key in ("samples", "batches", "workers", "p", "seed", "radius", "replica_offset"):
            if payload["run"][key] != expected[key]:
                raise ValueError(f"{response_path.name} differs from shard manifest for {key}")
        shard_rows = read_batches(batches_path)
        if len(shard_rows) != expected["batches"] or sum(row["samples"] for row in shard_rows) != expected["samples"]:
            raise ValueError(f"{batches_path.name} has wrong sufficient-statistic coverage")
        payloads.append(payload)
        rows.extend(shard_rows)
    return merge_payloads(payloads, rows, execution), rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution-manifest", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batches-output", type=Path, required=True)
    args = parser.parse_args()
    payload, _ = load_and_merge(args.execution_manifest, args.input_dir)
    # Re-read and canonicalize the rows through the same verification path.
    execution = json.loads(args.execution_manifest.read_text())
    rows = []
    for shard_path in execution["shard_manifests"]:
        shard = json.loads((args.execution_manifest.parents[1] / shard_path).read_text())
        rows.extend(read_batches(args.input_dir / Path(shard["run"]["batches_output"]).name))
    rows.sort(key=lambda row: row["replica_first"])
    for index, row in enumerate(rows):
        row["batch"] = index
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    write_batches(args.batches_output, rows)
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
