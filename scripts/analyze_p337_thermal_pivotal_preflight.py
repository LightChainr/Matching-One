#!/usr/bin/env python3
"""Validate and descriptively summarize the frozen P337 pivotal preflight.

This scorer does not estimate an ensemble mean, standard error, p-value,
centered covariance, moving-root response, or full J2.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


SEEDS = {32: 2026083123593201, 64: 2026083123596401}
COUNTERS = set(range(32))
PAIRS = set(range(32))
EXPECTED_PAIR_SITE_TOTAL = {32: 32 * 32 * 1024, 64: 32 * 32 * 4096}
BIT_NAMES = {
    1: "two_bridge_persistent",
    2: "shared_transition_or_merger",
    4: "kernel_preserving_topological",
    8: "kernel_changed_kernel_only_or_joint",
}

CONFIG_COLUMNS = {
    "L", "seed", "configuration_counter", "occupation_fnv1a64", "K", "q", "E",
    "sum_g16", "eligible_pairs", "nonzero_pairs",
    *(f"sum_g16_shared{i}" for i in range(5)),
    *(f"pairs_shared{i}" for i in range(5)),
}
SHELL_COLUMNS = {
    "L", "seed", "configuration_counter", "pair_ordinal", "x", "y",
    "anchor_x", "anchor_y", "orientation", "shell", "shell_lower", "shell_upper",
    "relation_mask", "carrier_mask", "sites", "eligible_both",
    "persistent_s2_count", "shared_transition_count",
    "kernel_preserving_topological_count", "kernel_only_count", "joint_count",
    "sum_g16_0", "sum_g16_1", "sum_abs_g16_0", "sum_abs_g16_1",
    "sum_delta_g16", "sum_abs_delta_g16", "sum_q0", "sum_q1",
    "sum_delta_q", "sum_abs_delta_q", "sum_E0", "sum_E1", "sum_delta_E",
    "sum_abs_delta_E", "sum_q_observable_num32", "sum_abs_q_observable_num32",
    "sum_q_kernel_num32", "sum_abs_q_kernel_num32",
    "sum_q_product_delta_num16", "sum_abs_q_product_delta_num16",
    "sum_E_observable_num32", "sum_abs_E_observable_num32",
    "sum_E_kernel_num32", "sum_abs_E_kernel_num32",
    "sum_E_product_delta_num16", "sum_abs_E_product_delta_num16",
    "q_identity_residual_num32", "E_identity_residual_num32",
    "q_identity_max_abs", "E_identity_max_abs", "sum_shared0", "sum_shared1",
}

SUMMARY_FIELDS = (
    "sites", "eligible_both", "persistent_s2_count", "shared_transition_count",
    "kernel_preserving_topological_count", "kernel_only_count", "joint_count",
    "sum_delta_g16", "sum_abs_delta_g16",
    "sum_q_observable_num32", "sum_abs_q_observable_num32",
    "sum_q_kernel_num32", "sum_abs_q_kernel_num32",
    "sum_q_product_delta_num16", "sum_abs_q_product_delta_num16",
    "sum_E_observable_num32", "sum_abs_E_observable_num32",
    "sum_E_kernel_num32", "sum_abs_E_kernel_num32",
    "sum_E_product_delta_num16", "sum_abs_E_product_delta_num16",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--shell", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-md", required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def integer(row: dict[str, str], column: str) -> int:
    text = row[column]
    if text.strip() != text or not text:
        raise ValueError(f"invalid integer in {column}: {text!r}")
    value = int(text)
    if str(value) != text:
        raise ValueError(f"noncanonical integer in {column}: {text!r}")
    return value


def read_csv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or []))
            raise ValueError(f"{path}: missing columns {missing}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path}: empty CSV")
    return rows


def expected_pair(L: int, ordinal: int) -> tuple[int, int, int, int, str]:
    r = L // 4
    anchor_ordinal, direction = divmod(ordinal, 2)
    j, i = divmod(anchor_ordinal, 4)
    ax, ay = i * r, j * r
    x = ay * L + ax
    if direction == 0:
        y, orientation = ay * L + ((ax + r) % L), "H"
    else:
        y, orientation = ((ay + r) % L) * L + ax, "V"
    return x, y, ax, ay, orientation


def expected_shell_bounds(L: int, shell: int) -> tuple[int, int]:
    if shell == 0:
        return 0, 0
    if shell < 0 or shell > 30:
        raise ValueError(f"invalid shell {shell}")
    return 1 << (shell - 1), min(L // 2, (1 << shell) - 1)


def validate_config(rows: list[dict[str, str]]) -> dict[tuple[int, int], dict[str, str]]:
    indexed: dict[tuple[int, int], dict[str, str]] = {}
    observed = defaultdict(set)
    for row in rows:
        L, seed, counter = integer(row, "L"), integer(row, "seed"), integer(row, "configuration_counter")
        if L not in SEEDS or seed != SEEDS[L] or counter not in COUNTERS:
            raise ValueError(f"unexpected config key {(L, seed, counter)}")
        key = (L, counter)
        if key in indexed:
            raise ValueError(f"duplicate config row {key}")
        indexed[key] = row
        observed[L].add(counter)
        N = L * L
        K, q, E = integer(row, "K"), integer(row, "q"), integer(row, "E")
        if not 0 <= K <= N or q not in {-1, 0, 1} or E != q * q:
            raise ValueError(f"invalid K/q/E at {key}")
        eligible = integer(row, "eligible_pairs")
        nonzero = integer(row, "nonzero_pairs")
        pair_counts = [integer(row, f"pairs_shared{i}") for i in range(5)]
        shared_sums = [integer(row, f"sum_g16_shared{i}") for i in range(5)]
        if sum(pair_counts) != eligible or not 0 <= nonzero <= eligible <= 32:
            raise ValueError(f"invalid original pair counts at {key}")
        if sum(shared_sums) != integer(row, "sum_g16"):
            raise ValueError(f"invalid original shared additivity at {key}")
    if set(indexed) != {(L, counter) for L in SEEDS for counter in COUNTERS}:
        raise ValueError("config CSV is not exactly 2x32 frozen counters")
    if len(rows) != 64:
        raise ValueError(f"expected 64 config rows, got {len(rows)}")
    return indexed


def relation_group(mask: int) -> str:
    if mask & 3:
        return "endpoint"
    if mask & 12:
        return "square_NN"
    return "external"


def blank_summary() -> dict[str, int]:
    return {field: 0 for field in SUMMARY_FIELDS}


def add_summary(target: dict[str, int], row: dict[str, str]) -> None:
    for field in SUMMARY_FIELDS:
        target[field] += integer(row, field)


def validate_shell(rows: list[dict[str, str]]) -> dict[str, Any]:
    unique_rows: set[tuple[int, ...]] = set()
    pair_sites: dict[tuple[int, int, int], int] = defaultdict(int)
    size_sites: dict[int, int] = defaultdict(int)
    row_count_by_L: dict[int, int] = defaultdict(int)
    by_L: dict[int, dict[str, int]] = defaultdict(blank_summary)
    by_shell: dict[tuple[int, int], dict[str, int]] = defaultdict(blank_summary)
    by_mask: dict[tuple[int, int], dict[str, int]] = defaultdict(blank_summary)
    by_relation: dict[tuple[int, str], dict[str, int]] = defaultdict(blank_summary)
    by_bit: dict[tuple[int, int], dict[str, int]] = defaultdict(blank_summary)

    signed_absolute_pairs = (
        ("sum_g16_0", "sum_abs_g16_0"), ("sum_g16_1", "sum_abs_g16_1"),
        ("sum_delta_g16", "sum_abs_delta_g16"),
        ("sum_delta_q", "sum_abs_delta_q"), ("sum_delta_E", "sum_abs_delta_E"),
        ("sum_q_observable_num32", "sum_abs_q_observable_num32"),
        ("sum_q_kernel_num32", "sum_abs_q_kernel_num32"),
        ("sum_q_product_delta_num16", "sum_abs_q_product_delta_num16"),
        ("sum_E_observable_num32", "sum_abs_E_observable_num32"),
        ("sum_E_kernel_num32", "sum_abs_E_kernel_num32"),
        ("sum_E_product_delta_num16", "sum_abs_E_product_delta_num16"),
    )

    for row in rows:
        L = integer(row, "L")
        seed = integer(row, "seed")
        counter = integer(row, "configuration_counter")
        pair = integer(row, "pair_ordinal")
        shell = integer(row, "shell")
        relation = integer(row, "relation_mask")
        mask = integer(row, "carrier_mask")
        if L not in SEEDS or seed != SEEDS[L] or counter not in COUNTERS or pair not in PAIRS:
            raise ValueError(f"unexpected shell key {(L, seed, counter, pair)}")
        if not 0 <= mask <= 15 or relation not in {0, 1, 2, 4, 8}:
            raise ValueError(f"invalid relation/carrier mask in {(L, counter, pair, shell)}")
        key = (L, counter, pair, shell, relation, mask)
        if key in unique_rows:
            raise ValueError(f"duplicate aggregate row {key}")
        unique_rows.add(key)

        x, y, ax, ay, orientation = expected_pair(L, pair)
        if (integer(row, "x"), integer(row, "y"), integer(row, "anchor_x"),
                integer(row, "anchor_y"), row["orientation"]) != (x, y, ax, ay, orientation):
            raise ValueError(f"pair definition mismatch at {(L, counter, pair)}")
        lower, upper = expected_shell_bounds(L, shell)
        if (integer(row, "shell_lower"), integer(row, "shell_upper")) != (lower, upper):
            raise ValueError(f"shell bounds mismatch at {key}")
        if relation in {1, 2} and shell != 0:
            raise ValueError(f"endpoint outside shell0 at {key}")
        if relation in {4, 8} and shell != 1:
            raise ValueError(f"square-NN outside shell1 at {key}")

        sites = integer(row, "sites")
        if sites <= 0 or not 0 <= integer(row, "eligible_both") <= sites:
            raise ValueError(f"invalid sites/eligible count at {key}")
        for bit, count_column in (
            (1, "persistent_s2_count"), (2, "shared_transition_count"),
            (4, "kernel_preserving_topological_count"),
        ):
            expected = sites if mask & bit else 0
            if integer(row, count_column) != expected:
                raise ValueError(f"carrier bit/count mismatch for bit {bit} at {key}")
        kernel_only, joint = integer(row, "kernel_only_count"), integer(row, "joint_count")
        if kernel_only < 0 or joint < 0 or kernel_only + joint != (sites if mask & 8 else 0):
            raise ValueError(f"kernel-only/joint partition mismatch at {key}")

        if integer(row, "sum_g16_1") - integer(row, "sum_g16_0") != integer(row, "sum_delta_g16"):
            raise ValueError(f"g delta mismatch at {key}")
        if integer(row, "sum_q1") - integer(row, "sum_q0") != integer(row, "sum_delta_q"):
            raise ValueError(f"q delta mismatch at {key}")
        if integer(row, "sum_E1") - integer(row, "sum_E0") != integer(row, "sum_delta_E"):
            raise ValueError(f"E delta mismatch at {key}")
        if (integer(row, "sum_q_observable_num32") + integer(row, "sum_q_kernel_num32") !=
                2 * integer(row, "sum_q_product_delta_num16")):
            raise ValueError(f"q midpoint identity mismatch at {key}")
        if (integer(row, "sum_E_observable_num32") + integer(row, "sum_E_kernel_num32") !=
                2 * integer(row, "sum_E_product_delta_num16")):
            raise ValueError(f"E midpoint identity mismatch at {key}")
        for column in ("q_identity_residual_num32", "E_identity_residual_num32",
                       "q_identity_max_abs", "E_identity_max_abs"):
            if integer(row, column) != 0:
                raise ValueError(f"nonzero identity control {column} at {key}")
        for signed, absolute_column in signed_absolute_pairs:
            absolute_value = integer(row, absolute_column)
            if absolute_value < 0 or absolute_value < abs(integer(row, signed)):
                raise ValueError(f"invalid pair/site absolute sum {absolute_column} at {key}")

        pair_sites[(L, counter, pair)] += sites
        size_sites[L] += sites
        row_count_by_L[L] += 1
        add_summary(by_L[L], row)
        add_summary(by_shell[(L, shell)], row)
        add_summary(by_mask[(L, mask)], row)
        add_summary(by_relation[(L, relation_group(relation))], row)
        for bit in BIT_NAMES:
            if mask & bit:
                add_summary(by_bit[(L, bit)], row)

    expected_keys = {(L, counter, pair) for L in SEEDS for counter in COUNTERS for pair in PAIRS}
    if set(pair_sites) != expected_keys:
        missing = len(expected_keys - set(pair_sites))
        raise ValueError(f"missing {missing} configuration/pair partitions")
    for key, count in pair_sites.items():
        if count != key[0] * key[0]:
            raise ValueError(f"configuration/pair {key} has {count} sites, expected {key[0]**2}")
    if dict(size_sites) != EXPECTED_PAIR_SITE_TOTAL:
        raise ValueError(f"pair/site totals mismatch: {dict(size_sites)}")

    def records(mapping: dict[Any, dict[str, int]], names: tuple[str, ...]) -> list[dict[str, Any]]:
        result = []
        for key in sorted(mapping, key=lambda value: value if isinstance(value, tuple) else (value,)):
            keys = key if isinstance(key, tuple) else (key,)
            record = dict(zip(names, keys))
            record.update(mapping[key])
            result.append(record)
        return result

    return {
        "row_count": len(rows),
        "row_count_by_L": {str(key): value for key, value in sorted(row_count_by_L.items())},
        "pair_site_total": sum(size_sites.values()),
        "pair_site_total_by_L": {str(key): value for key, value in sorted(size_sites.items())},
        "by_L": records(by_L, ("L",)),
        "by_shell": records(by_shell, ("L", "shell")),
        "by_carrier_mask": records(by_mask, ("L", "carrier_mask")),
        "by_carrier_bit": records(by_bit, ("L", "carrier_bit")),
        "by_relation": records(by_relation, ("L", "relation")),
    }


def validate_metadata(metadata: dict[str, Any], config_path: Path, shell_path: Path,
                      row_count_by_L: dict[str, int]) -> None:
    if metadata.get("schema") != "matching-one.p337-thermal-pivotal-preflight.run.v1":
        raise ValueError("unexpected metadata schema")
    if metadata.get("status") != "completed_frozen_64_configuration_replay":
        raise ValueError("producer metadata is not completed")
    if metadata.get("configurations_per_size") != 32 or metadata.get("pairs_per_configuration") != 32:
        raise ValueError("metadata count mismatch")
    if metadata.get("kernel_rows") != 1874:
        raise ValueError("metadata kernel row mismatch")
    if not isinstance(metadata.get("threads"), int) or not 1 <= metadata["threads"] <= 64:
        raise ValueError("metadata thread count mismatch")
    for timing in ("total_wall_seconds", "total_cpu_seconds"):
        if not isinstance(metadata.get(timing), (int, float)) or not math.isfinite(metadata[timing]) or metadata[timing] < 0:
            raise ValueError(f"invalid metadata {timing}")
    if not isinstance(metadata.get("peak_RSS_bytes"), int) or metadata["peak_RSS_bytes"] < 0:
        raise ValueError("invalid metadata peak_RSS_bytes")
    controls = metadata.get("tiny_controls", {})
    if controls.get("kernel_preserving_topological") is not True or controls.get("kernel_only") is not True:
        raise ValueError("producer tiny controls did not pass")
    for label, L in (("L32", 32), ("L64", 64)):
        block = metadata.get(label, {})
        if block.get("seed") != SEEDS[L] or block.get("configuration_rows") != 32:
            raise ValueError(f"metadata {label} replay mismatch")
        if block.get("shell_rows") != row_count_by_L[str(L)]:
            raise ValueError(f"metadata {label} shell-row mismatch")
        for timing in ("wall_seconds", "cpu_seconds"):
            if not isinstance(block.get(timing), (int, float)) or not math.isfinite(block[timing]) or block[timing] < 0:
                raise ValueError(f"invalid {label} {timing}")
    # The frozen runner may copy the immutable raw files from its temporary
    # execution directory into the repository.  Preserve the original path in
    # metadata, but require the artifact basename to remain unchanged.
    if Path(metadata.get("configuration_output", "")).name != config_path.name:
        raise ValueError("metadata configuration basename mismatch after relocation")
    if Path(metadata.get("shell_output", "")).name != shell_path.name:
        raise ValueError("metadata shell basename mismatch after relocation")


def fmt_scaled(value: int, denominator: int) -> str:
    return f"{value / denominator:+.9g}"


def make_report(result: dict[str, Any]) -> str:
    lines = [
        "# P337 thermal/pivotal preflight descriptive score",
        "",
        "Status: **completed descriptive preflight only**.",
        "",
        "All 64 frozen counters, 32 original pairs per configuration and 5,242,880 "
        "pair/site callbacks passed the structural checks. This report contains no "
        "significance test, ensemble centering, moving-root reconstruction, full J2, "
        "field ratio or independent evidence claim.",
        "",
        "## Validation",
        "",
        f"- configuration rows: `{result['validation']['configuration_rows']}`",
        f"- shell rows: `{result['validation']['shell_rows']}`",
        f"- pair/site callbacks: `{result['validation']['pair_site_total']}`",
        "- every configuration/pair partitions exactly N sites across complete carrier masks",
        "- all q/E midpoint residual sums and maxima are zero",
        "- pair/site absolute values were checked before descriptive aggregation",
        "",
        "## Signed and absolute primitives by L",
        "",
        "| L | sites | Δg signed/abs | q observable signed/abs | q kernel signed/abs | E observable signed/abs | E kernel signed/abs |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in result["summaries"]["by_L"]:
        lines.append(
            f"| {row['L']} | {row['sites']} | "
            f"{fmt_scaled(row['sum_delta_g16'],16)} / {row['sum_abs_delta_g16']/16:.9g} | "
            f"{fmt_scaled(row['sum_q_observable_num32'],32)} / {row['sum_abs_q_observable_num32']/32:.9g} | "
            f"{fmt_scaled(row['sum_q_kernel_num32'],32)} / {row['sum_abs_q_kernel_num32']/32:.9g} | "
            f"{fmt_scaled(row['sum_E_observable_num32'],32)} / {row['sum_abs_E_observable_num32']/32:.9g} | "
            f"{fmt_scaled(row['sum_E_kernel_num32'],32)} / {row['sum_abs_E_kernel_num32']/32:.9g} |"
        )
    lines += [
        "",
        "The JSON preserves the same fields by dyadic shell, complete carrier mask, "
        "nonexclusive carrier bit, and endpoint/square-NN/external relation. Carrier-bit "
        "views overlap by construction and must not be added as independent or exhaustive votes.",
        "",
        "## Interpretation boundary",
        "",
        "These are uncentered finite callback sums over a deterministic subset replay. "
        "They do not estimate the centered observable-pivot or kernel-pivot expectation, "
        "and they do not apply the baseline q/E jets, common root, denominator or slope terms. "
        "No absence or sign in this 64-counter report authorizes more counters or a new seed.",
        "",
        "Frozen scorer CLI:",
        "",
        "```bash",
        "python3 scripts/analyze_p337_thermal_pivotal_preflight.py \\",
        "  --config results/p337-thermal-pivotal-preflight/raw/preflight.config.csv \\",
        "  --shell results/p337-thermal-pivotal-preflight/raw/preflight.shell.csv \\",
        "  --metadata results/p337-thermal-pivotal-preflight/raw/preflight.metadata.json \\",
        "  --output-json results/p337-thermal-pivotal-preflight/latest.json \\",
        "  --output-md results/p337-thermal-pivotal-preflight/REPORT.md",
        "```",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    config_path, shell_path, metadata_path = map(Path, (args.config, args.shell, args.metadata))
    output_json, output_md = Path(args.output_json), Path(args.output_md)
    if output_json.exists() or output_md.exists():
        raise FileExistsError("output JSON/Markdown exists; refusing overwrite")
    for path in (config_path, shell_path, metadata_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    config_rows = read_csv(config_path, CONFIG_COLUMNS)
    shell_rows = read_csv(shell_path, SHELL_COLUMNS)
    with metadata_path.open(encoding="utf-8") as handle:
        metadata = json.load(handle)
    validate_config(config_rows)
    summaries = validate_shell(shell_rows)
    validate_metadata(metadata, config_path, shell_path, summaries["row_count_by_L"])

    result = {
        "schema": "matching-one.p337-thermal-pivotal-preflight.score.v1",
        "status": "completed_descriptive_preflight_only",
        "inputs": {
            "config": str(config_path), "config_sha256": sha256(config_path),
            "shell": str(shell_path), "shell_sha256": sha256(shell_path),
            "metadata": str(metadata_path), "metadata_sha256": sha256(metadata_path),
            "relocated_from": {
                "config": metadata.get("configuration_output"),
                "shell": metadata.get("shell_output"),
            },
        },
        "validation": {
            "configuration_rows": len(config_rows),
            "shell_rows": summaries["row_count"],
            "shell_rows_by_L": summaries["row_count_by_L"],
            "pair_site_total": summaries["pair_site_total"],
            "pair_site_total_by_L": summaries["pair_site_total_by_L"],
            "frozen_counters_and_seeds_exact": True,
            "each_configuration_pair_partitions_N_sites": True,
            "carrier_masks_partition_records": True,
            "midpoint_identity_residuals_zero": True,
            "pair_site_absolute_rule_checked": True,
        },
        "producer_metadata": metadata,
        "carrier_bit_names": {str(key): value for key, value in BIT_NAMES.items()},
        "summaries": {key: value for key, value in summaries.items()
                      if key not in {"row_count", "row_count_by_L", "pair_site_total", "pair_site_total_by_L"}},
        "interpretation_boundary": (
            "Uncentered descriptive signed/absolute sums from the frozen 64-counter subset; "
            "no significance, ensemble centering, full original-U J2, field identification, "
            "independent evidence, new-seed permission or production escalation."
        ),
    }
    report = make_report(result)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    output_md.write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()
