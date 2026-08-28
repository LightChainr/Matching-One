#!/usr/bin/env python3
"""Score the frozen Issue #50 N=290 third-lineage fixed-p result.

The scorer is intentionally protocol-specific and standard-library-only.  It
accepts one ``gaussian_orientation_mc`` batches/metadata pair, validates its
provenance and counter range, and scores the ``either/matching_function``
orientation contrast against the prediction frozen before the run.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shlex
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple


EXPECTED_N = 290
EXPECTED_FIRST = (13, 11)
EXPECTED_SECOND = (17, 1)
CHANNELS = ("cross", "both", "either", "direction_0", "direction_1")
TARGET_DELTA_M = -1.3765640041065354e-4
SOURCE_SE = 2.499658241589821e-5
COMMIT_RE = re.compile(r"^[0-9a-fA-F]{40}$")


def command_option(command: str, name: str) -> str:
    tokens = shlex.split(command)
    matches = [tokens[i + 1] for i, token in enumerate(tokens[:-1]) if token == name]
    if len(matches) != 1:
        raise ValueError("metadata command must contain exactly one {} option".format(name))
    return matches[0]


def load_metadata(path: Path) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        metadata = json.load(handle)
    if not isinstance(metadata, dict):
        raise ValueError("metadata must be a JSON object")
    return metadata


def validate_metadata(metadata: Mapping[str, object]) -> Dict[str, object]:
    if metadata.get("engine") != "same-N Gaussian site-orientation discovery":
        raise ValueError("unexpected engine")
    designs = metadata.get("designs")
    if not isinstance(designs, list) or len(designs) != 1:
        raise ValueError("Issue #50 requires exactly one design")
    design = designs[0]
    if not isinstance(design, dict):
        raise ValueError("invalid design metadata")
    if int(design.get("N", -1)) != EXPECTED_N:
        raise ValueError("Issue #50 requires N=290")
    if tuple(design.get("first", ())) != EXPECTED_FIRST or tuple(design.get("second", ())) != EXPECTED_SECOND:
        raise ValueError("N=290 lineage order must be (13,11) first and (17,1) second")

    source_commit = str(metadata.get("git_commit", ""))
    if not COMMIT_RE.fullmatch(source_commit):
        raise ValueError("source commit must be a complete 40-hex commit id")
    try:
        samples = int(metadata["samples_per_pair"])
        batches = int(metadata["batches"])
        seed = int(metadata["seed"])
        counter_first = int(metadata["replica_counter_first"])
        counter_last = int(metadata["replica_counter_last_exclusive"])
        p_ref = float(metadata["p_ref"])
        command = str(metadata["command"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("missing or invalid run metadata") from error
    if samples <= 0 or batches < 2 or samples % batches:
        raise ValueError("samples must be positive and divisible by batches>=2")
    if seed < 0 or counter_first < 0 or counter_last - counter_first != samples:
        raise ValueError("seed/counter range is inconsistent with samples")
    if not 0.0 < p_ref < 1.0:
        raise ValueError("invalid fixed p_ref")
    if tuple(metadata.get("channels", ())) != CHANNELS:
        raise ValueError("metadata channel list differs from the frozen engine output")

    expected_options = {
        "--n": EXPECTED_N,
        "--samples": samples,
        "--batches": batches,
        "--seed": seed,
        "--replica-offset": counter_first,
    }
    for option, expected in expected_options.items():
        try:
            actual = int(command_option(command, option))
        except ValueError as error:
            raise ValueError("command provenance mismatch for {}".format(option)) from error
        if actual != expected:
            raise ValueError("command provenance mismatch for {}".format(option))
    if command_option(command, "--git-commit").lower() != source_commit.lower():
        raise ValueError("command/source commit mismatch")
    # The formal fixed-p command records p_ref explicitly.  Reject an implicit
    # default so the scored thermal coordinate is independently auditable.
    try:
        command_p_ref = float(command_option(command, "--p-ref"))
    except ValueError as error:
        raise ValueError("command provenance mismatch for --p-ref") from error
    if not math.isclose(command_p_ref, p_ref, rel_tol=0.0, abs_tol=1e-15):
        raise ValueError("command/metadata p_ref mismatch")
    return {
        "samples": samples,
        "batches": batches,
        "seed": seed,
        "counter_first": counter_first,
        "counter_last": counter_last,
        "p_ref": p_ref,
        "source_commit": source_commit.lower(),
    }


def read_batches(path: Path, run: Mapping[str, object]) -> Dict[Tuple[str, int], Dict[str, int]]:
    required = {
        "n", "batch", "samples", "p_ref", "channel", "a1", "b1", "a2", "b2",
        "first_primal_sum", "first_matching_sum", "second_primal_sum", "second_matching_sum",
    }
    rows: Dict[Tuple[str, int], Dict[str, int]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("batch CSV missing fields: " + ", ".join(sorted(missing)))
        for raw in reader:
            if int(raw["n"]) != EXPECTED_N:
                raise ValueError("batch CSV contains an N other than 290")
            if (int(raw["a1"]), int(raw["b1"])) != EXPECTED_FIRST or (
                int(raw["a2"]), int(raw["b2"])
            ) != EXPECTED_SECOND:
                raise ValueError("batch lineage order differs from metadata/protocol")
            channel = raw["channel"]
            batch = int(raw["batch"])
            if channel not in CHANNELS:
                raise ValueError("unknown channel")
            key = (channel, batch)
            if key in rows:
                raise ValueError("duplicate channel/batch row")
            samples = int(raw["samples"])
            if samples != int(run["samples"]) // int(run["batches"]):
                raise ValueError("batch samples disagree with metadata")
            if not math.isclose(float(raw["p_ref"]), float(run["p_ref"]), rel_tol=0.0, abs_tol=1e-15):
                raise ValueError("batch p_ref disagrees with metadata")
            values = {
                "samples": samples,
                "first_primal": int(raw["first_primal_sum"]),
                "first_matching": int(raw["first_matching_sum"]),
                "second_primal": int(raw["second_primal_sum"]),
                "second_matching": int(raw["second_matching_sum"]),
            }
            if any(value < 0 or value > samples for name, value in values.items() if name != "samples"):
                raise ValueError("indicator count outside [0,samples]")
            rows[key] = values
    expected_keys = {(channel, batch) for channel in CHANNELS for batch in range(int(run["batches"]))}
    if set(rows) != expected_keys:
        raise ValueError("batch CSV lacks the complete channel/batch grid")
    return rows


def mean_se(values: Sequence[float]) -> Tuple[float, float]:
    mean = math.fsum(values) / len(values)
    variance = math.fsum((value - mean) ** 2 for value in values)
    return mean, math.sqrt(variance / (len(values) * (len(values) - 1)))


def score(rows: Mapping[Tuple[str, int], Mapping[str, int]], run: Mapping[str, object]) -> Dict[str, object]:
    batch_values: List[float] = []
    for batch in range(int(run["batches"])):
        row = rows[("either", batch)]
        samples = row["samples"]
        first_matching_function = (row["first_primal"] - row["first_matching"]) / samples
        second_matching_function = (row["second_primal"] - row["second_matching"]) / samples
        batch_values.append(first_matching_function - second_matching_function)
    child, child_se = mean_se(batch_values)
    residual = child - TARGET_DELTA_M
    combined_se = math.hypot(child_se, SOURCE_SE)
    return {
        "protocol": "Issue #50 N=290 third Gaussian lineage fixed-p prospective score",
        "N": EXPECTED_N,
        "lineage_first": list(EXPECTED_FIRST),
        "lineage_second": list(EXPECTED_SECOND),
        "channel": "either",
        "sector": "matching_function",
        "definition": "(first_primal-first_matching)-(second_primal-second_matching)",
        "source_commit": run["source_commit"],
        "samples": run["samples"],
        "batches": run["batches"],
        "seed": run["seed"],
        "replica_counter_first": run["counter_first"],
        "replica_counter_last_exclusive": run["counter_last"],
        "p_ref": run["p_ref"],
        "child_delta_M": child,
        "child_sampling_se": child_se,
        "frozen_target_delta_M": TARGET_DELTA_M,
        "frozen_source_se": SOURCE_SE,
        "target_residual_child_minus_frozen": residual,
        "target_residual_combined_se": combined_se,
        "target_residual_z": residual / combined_se,
        "zero_residual": child,
        "zero_sampling_se": child_se,
        "zero_z": child / child_se,
        "target_chi_square_1df": (residual / combined_se) ** 2,
        "zero_chi_square_1df": (child / child_se) ** 2,
    }


def write_outputs(result: Mapping[str, object], json_path: Path, csv_path: Path) -> None:
    json_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    fields = list(result)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        row = dict(result)
        row["lineage_first"] = "13+11i"
        row["lineage_second"] = "17+1i"
        writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run = validate_metadata(load_metadata(args.metadata))
    result = score(read_batches(args.batches, run), run)
    write_outputs(result, args.json, args.csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
