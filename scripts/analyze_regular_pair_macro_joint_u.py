#!/usr/bin/env python3
"""Score the frozen canonical macro-window joint-U pilot.

The CLI exposes four geometry roles.  The frozen producer writes one combined
CSV per size, so the same size-level path can be supplied to both roles and is
filtered without breaking batch pairing.  This program does not generate
occupations or choose a window.  Roots, thermal jets, all four original-U
response terms, delete-one covariance and the two conditional fixed-power
contrasts are recomputed from the frozen contract.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import sys
import time
from typing import Any

import numpy as np
import scipy
from scipy.stats import t as student_t


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis/regular_pair_macro_joint_u_contract.json"
SIZES = (100, 400)
GEOMETRIES = ("axis", "tilted")
STRATA = ("total", "s2", "sge3")
SOURCE_ATOMS = ("H", "qH", "EH")
TERMS = ("direct_centered", "root_motion", "slope_source", "slope_root")


class InvalidInterface(Exception):
    """A declared input or scoreability gate failed."""

    def __init__(self, gate: str, reason: str, diagnostic: Any | None = None):
        super().__init__(reason)
        self.gate = gate
        self.reason = reason
        self.diagnostic = diagnostic


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def normalized_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def parse_integer(raw: str | None, *, field: str, path: Path) -> int:
    if raw is None or raw.strip() == "":
        raise InvalidInterface("raw_schema", f"{path}: empty integer field {field}")
    try:
        return int(raw)
    except ValueError as error:
        raise InvalidInterface("raw_schema", f"{path}: noninteger {field}={raw!r}") from error


def field_index(fieldnames: list[str], path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for name in fieldnames:
        key = normalized_name(name)
        if key in result and result[key] != name:
            raise InvalidInterface("raw_schema", f"{path}: ambiguous columns {result[key]!r}/{name!r}")
        result[key] = name
    return result


def find_column(index: dict[str, str], aliases: list[str], *, required: bool = True) -> str | None:
    matches = []
    for alias in aliases:
        key = normalized_name(alias)
        if key in index and index[key] not in matches:
            matches.append(index[key])
    if len(matches) > 1:
        raise InvalidInterface("raw_schema", f"aliases select multiple columns: {matches}")
    if matches:
        return matches[0]
    if required:
        raise InvalidInterface("raw_schema", f"missing column; accepted aliases: {aliases}")
    return None


STRATUM_TOKENS = {
    "total": ("total", "all"),
    "s2": ("s2", "exact2", "exactly2", "exactlytwosharedoccupiedcomponents",
           "s2exactlytwosharedoccupiedcomponents"),
    "sge3": ("sge3", "ge3", "atleast3", "atleastthree",
             "atleastthreesharedoccupiedcomponents", "sge3atleastthreesharedoccupiedcomponents"),
}


FIELD_ALIASES = {
    "sum_B16": ("sumB16", "B16sum"),
    "sum_qB16": ("sumqB16", "qB16sum"),
    "sum_EB16": ("sumEB16", "EB16sum"),
    "eligible_pair_count": ("eligiblepaircount", "paircount", "eligiblepairs"),
    "nonzero_pair_count": ("nonzeropaircount", "nonzeropairs"),
}


def wide_aliases(stratum: str, field: str) -> list[str]:
    bases = FIELD_ALIASES[field]
    candidates: list[str] = []
    if stratum == "total":
        candidates.extend(bases)
    for token in STRATUM_TOKENS[stratum]:
        for base in bases:
            candidates.extend((f"{token}_{base}", f"{base}_{token}"))
    return candidates


def stratum_name(raw: str, path: Path) -> str:
    key = normalized_name(raw)
    for name, tokens in STRATUM_TOKENS.items():
        if key in {normalized_name(token) for token in tokens}:
            return name
    raise InvalidInterface("raw_schema", f"{path}: unknown stratum {raw!r}")


def merge_metadata(target: dict[str, Any], update: dict[str, Any], source: str) -> None:
    """Merge without hiding conflicting scalar metadata."""
    for key, value in update.items():
        if key not in target:
            target[key] = value
        elif isinstance(target[key], dict) and isinstance(value, dict):
            merge_metadata(target[key], value, source)
        elif target[key] != value:
            # Different nested receipts can legitimately have unrelated copies;
            # keep both only under an explicit provenance key.
            target.setdefault("_metadata_conflicts", []).append(
                {"key": key, "first": target[key], "second": value, "source": source}
            )


def load_sidecar_metadata(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    candidates = [
        Path(str(path) + ".metadata.json"),
        path.with_suffix(".metadata.json"),
        path.parent / "metadata.json",
        path.parent / "run.json",
    ]
    metadata: dict[str, Any] = {}
    inputs: list[dict[str, Any]] = []
    seen: set[Path] = set()
    for candidate in candidates:
        if candidate in seen or not candidate.exists() or candidate == path:
            continue
        seen.add(candidate)
        payload = candidate.read_bytes()
        try:
            obj = json.loads(payload)
        except json.JSONDecodeError:
            continue
        if isinstance(obj, dict):
            merge_metadata(metadata, obj, display_path(candidate))
            inputs.append({"path": display_path(candidate), "sha256": sha256(payload), "bytes": len(payload)})
    return metadata, inputs


def geometry_value_ok(raw: str, role: str, coordinates: tuple[int, int]) -> bool:
    key = normalized_name(raw)
    aliases = {
        "axis": {"axis", "first", "a", normalized_name(f"{coordinates[0]},{coordinates[1]}"),
                 normalized_name(f"{coordinates[0]}_{coordinates[1]}")},
        "tilted": {"tilted", "second", "oblique", "b",
                   normalized_name(f"{coordinates[0]},{coordinates[1]}"),
                   normalized_name(f"{coordinates[0]}_{coordinates[1]}")},
    }
    return key in aliases[role]


def geometry_role(raw: str, n: int, contract: dict[str, Any], path: Path) -> str:
    specification = next(
        (item for item in contract["model"]["geometries"] if int(item["N"]) == n),
        None,
    )
    if specification is None:
        raise InvalidInterface("contract", f"no geometry specification for N={n}")
    matches = [
        role for role in GEOMETRIES
        if geometry_value_ok(raw, role, tuple(specification[role]))
    ]
    if len(matches) != 1:
        raise InvalidInterface("raw_keys", f"{path}: unknown or ambiguous geometry {raw!r}")
    return matches[0]


def quotient_representatives(a: int, b: int) -> tuple[list[tuple[int, int]], dict[int, int]]:
    """Reproduce the producer's E-then-N quotient BFS exactly."""
    n = a * a + b * b

    def key(x: int, y: int) -> int:
        return n * ((a * x + b * y) % n) + ((-b * x + a * y) % n)

    representatives = [(0, 0)]
    coordinate_index = {key(0, 0): 0}
    cursor = 0
    while cursor < len(representatives):
        x, y = representatives[cursor]
        for dx, dy in ((1, 0), (0, 1)):
            item = key(x + dx, y + dy)
            if item not in coordinate_index:
                coordinate_index[item] = len(representatives)
                representatives.append((x + dx, y + dy))
        cursor += 1
    if len(representatives) != n:
        raise InvalidInterface(
            "window_geometry", f"({a},{b}) quotient BFS found {len(representatives)} rather than {n} classes"
        )
    return representatives, coordinate_index


def expected_window_rows(a: int, b: int) -> tuple[list[tuple[int, int, int, int]], bool]:
    """Return (vertex,r2,dx,dy) rows and the incident-port disjointness gate."""
    n = a * a + b * b
    representatives, coordinate_index = quotient_representatives(a, b)

    def key(x: int, y: int) -> int:
        return n * ((a * x + b * y) % n) + ((-b * x + a * y) % n)

    rows: list[tuple[int, int, int, int]] = []
    for vertex, (x, y) in enumerate(representatives[1:], start=1):
        du, dv = a * x + b * y, -b * x + a * y
        mf, nf = du // n, dv // n
        candidates = []
        for m in (mf, mf + 1):
            for k in (nf, nf + 1):
                dx, dy = x - m * a + k * b, y - m * b - k * a
                candidates.append((dx * dx + dy * dy, dx, dy))
        r2, dx, dy = min(candidates)
        if 16 * r2 >= n and 25 * r2 <= 4 * n:
            rows.append((vertex, r2, dx, dy))
    rows.sort(key=lambda item: (item[1], item[2], item[3]))

    neighbors: list[tuple[int, int, int, int]] = []
    for x, y in representatives:
        neighbors.append(tuple(
            coordinate_index[key(x + dx, y + dy)]
            for dx, dy in ((0, 1), (1, 0), (0, -1), (-1, 0))
        ))

    def ports(vertex: int) -> set[int]:
        north, east, south, west = neighbors[vertex]
        del north, east
        return {2 * vertex, 2 * vertex + 1, 2 * south, 2 * west + 1}

    origin_ports = ports(0)
    disjoint = all(not origin_ports.intersection(ports(row[0])) for row in rows)
    return rows, disjoint


def validate_window_table(path: Path, a: int, b: int) -> dict[str, Any]:
    try:
        payload = path.read_bytes()
    except OSError as error:
        raise InvalidInterface("window_table", f"cannot read frozen window table {path}: {error}") from error
    try:
        reader = csv.DictReader(payload.decode("utf-8").splitlines())
    except UnicodeDecodeError as error:
        raise InvalidInterface("window_table", f"{path}: window table is not UTF-8") from error
    index = field_index(reader.fieldnames or [], path)
    columns = {
        "ordinal": find_column(index, ["ordinal", "index"]),
        "vertex": find_column(index, ["vertex", "quotient_vertex"]),
        "r2": find_column(index, ["r2", "squared_distance"]),
        "dx": find_column(index, ["canonical_dx", "dx"]),
        "dy": find_column(index, ["canonical_dy", "dy"]),
    }
    observed: list[tuple[int, int, int, int]] = []
    for ordinal, raw in enumerate(reader):
        recorded = parse_integer(raw[columns["ordinal"]], field="ordinal", path=path)
        if recorded != ordinal:
            raise InvalidInterface("window_table", f"{path}: nonconsecutive ordinal {recorded} at row {ordinal}")
        observed.append(tuple(
            parse_integer(raw[columns[name]], field=name, path=path)
            for name in ("vertex", "r2", "dx", "dy")
        ))
    expected, disjoint = expected_window_rows(a, b)
    if observed != expected:
        mismatch = next(
            (index for index, pair in enumerate(zip(observed, expected)) if pair[0] != pair[1]),
            min(len(observed), len(expected)),
        )
        raise InvalidInterface(
            "exact_window_membership",
            f"{path}: saved table differs from the exact complete quotient window",
            {"first_mismatch": mismatch, "observed_rows": len(observed), "expected_rows": len(expected)},
        )
    if not disjoint:
        raise InvalidInterface("disjoint_edge_ports", f"{path}: an accepted pair shares an incident edge port")
    return {
        "path": display_path(path), "sha256": sha256(payload), "bytes": len(payload),
        "displacement_count": len(observed), "exact_membership": True,
        "incident_edge_ports_disjoint": True,
    }


def parse_comments(lines: list[str]) -> tuple[list[str], dict[str, Any]]:
    body, metadata = [], {}
    for line in lines:
        if not line.lstrip().startswith("#"):
            if line.strip():
                body.append(line)
            continue
        content = line.lstrip()[1:].strip()
        if not content:
            continue
        try:
            item = json.loads(content)
        except json.JSONDecodeError:
            item = None
        if isinstance(item, dict):
            merge_metadata(metadata, item, "CSV comment")
        elif "=" in content:
            key, value = content.split("=", 1)
            metadata[key.strip()] = value.strip()
    return body, metadata


def load_csv(path: Path, *, n: int, role: str, coordinates: tuple[int, int], contract: dict[str, Any]) -> dict[str, Any]:
    try:
        payload = path.read_bytes()
        text = payload.decode("utf-8")
    except OSError as error:
        raise InvalidInterface("raw_input", f"cannot read {path}: {error}") from error
    except UnicodeDecodeError as error:
        raise InvalidInterface("raw_input", f"{path}: raw CSV is not UTF-8") from error
    lines, comment_metadata = parse_comments(text.splitlines())
    reader = csv.DictReader(lines)
    names = reader.fieldnames or []
    index = field_index(names, path)
    columns = {
        "N": find_column(index, ["N", "n"]),
        "batch": find_column(index, ["batch", "batch_id"]),
        "geometry": find_column(index, ["geometry", "orientation", "direction"], required=False),
        "a": find_column(index, ["a", "geometry_a"], required=False),
        "b": find_column(index, ["b", "geometry_b"], required=False),
        "K": find_column(index, ["K", "k"]),
        "count": find_column(index, ["count", "samples", "sample_count"]),
        "sum_q": find_column(index, ["sum_q", "q_sum"]),
        "sum_E": find_column(index, ["sum_E", "sum_e", "E_sum", "e_sum"]),
        "stratum": find_column(index, ["stratum", "support_stratum"], required=False),
        "s_le1_pair_count": find_column(index, ["s_le1_pair_count", "sle1_pair_count"], required=False),
        "s_le1_nonzero": find_column(index, ["s_le1_nonzero_g16_count", "sle1_nonzero_g16_count"], required=False),
        "additivity_control": find_column(index, ["total_minus_s2_minus_sge3_B16", "B16_additivity_residual"], required=False),
    }
    long_layout = columns["stratum"] is not None
    source_columns: dict[str, dict[str, str]] = {}
    for stratum in STRATA:
        source_columns[stratum] = {}
        for field in FIELD_ALIASES:
            aliases = list(FIELD_ALIASES[field]) if long_layout else wide_aliases(stratum, field)
            source_columns[stratum][field] = find_column(index, aliases)

    batches = int(contract["pilot_sampling"]["batches_per_size"])
    counts = np.zeros((batches, n + 1), dtype=np.int64)
    baseline = np.zeros((batches, n + 1, 2), dtype=np.int64)
    source = np.zeros((batches, n + 1, len(STRATA), len(SOURCE_ATOMS)), dtype=np.int64)
    pair_counts = np.zeros((batches, n + 1, len(STRATA), 2), dtype=np.int64)
    controls = np.zeros((batches, n + 1, 3), dtype=np.int64)
    seen_base: set[tuple[int, int]] = set()
    seen_source: set[tuple[int, int, str]] = set()
    seen_controls: set[tuple[int, int, str]] = set()
    row_metadata: dict[str, Any] = {}
    metadata_aliases = {
        "window_displacement_count": ["window_displacement_count", "window_displacements", "displacement_count"],
        "anchors_per_configuration": ["anchors_per_configuration", "anchor_count", "anchors"],
        "occupation_seed": ["occupation_seed", "rng_occupation_seed"],
        "anchor_seed": ["anchor_seed", "rng_anchor_seed"],
        "source_denominator": ["source_denominator", "H_denominator"],
    }
    metadata_columns = {
        name: find_column(index, aliases, required=False)
        for name, aliases in metadata_aliases.items()
    }

    for raw in reader:
        row_n = parse_integer(raw[columns["N"]], field="N", path=path)
        if row_n != n:
            raise InvalidInterface("raw_keys", f"{path}: expected N={n}, found N={row_n}")
        if columns["geometry"]:
            row_role = geometry_role(raw[columns["geometry"]], n, contract, path)
            if row_role != role:
                # The frozen producer deliberately writes both paired geometries
                # to one size-level CSV.  Passing that same path to both CLI
                # roles must preserve, rather than split, the pairing.
                continue
        if (columns["a"] is None) != (columns["b"] is None):
            raise InvalidInterface("raw_schema", f"{path}: geometry coordinates require both a and b")
        if columns["a"] is not None:
            row_coordinates = (
                parse_integer(raw[columns["a"]], field="a", path=path),
                parse_integer(raw[columns["b"]], field="b", path=path),
            )
            if row_coordinates != coordinates:
                raise InvalidInterface(
                    "raw_keys", f"{path}: {role} coordinates {row_coordinates} differ from {coordinates}"
                )
        batch = parse_integer(raw[columns["batch"]], field="batch", path=path)
        k = parse_integer(raw[columns["K"]], field="K", path=path)
        if not 0 <= batch < batches or not 0 <= k <= n:
            raise InvalidInterface("raw_keys", f"{path}: invalid key N={row_n}, batch={batch}, K={k}")
        key = (batch, k)
        base_values = (
            parse_integer(raw[columns["count"]], field="count", path=path),
            parse_integer(raw[columns["sum_q"]], field="sum_q", path=path),
            parse_integer(raw[columns["sum_E"]], field="sum_E", path=path),
        )
        if key not in seen_base:
            counts[key] = base_values[0]
            baseline[key] = base_values[1:]
            seen_base.add(key)
        elif (counts[key], *baseline[key]) != base_values:
            raise InvalidInterface("raw_keys", f"{path}: repeated baseline differs at {key}")

        selected = (stratum_name(raw[columns["stratum"]], path),) if long_layout else STRATA
        for stratum in selected:
            source_key = (batch, k, stratum)
            if source_key in seen_source:
                raise InvalidInterface("raw_keys", f"{path}: duplicate source row {source_key}")
            si = STRATA.index(stratum)
            fields = source_columns[stratum]
            for ai, field in enumerate(("sum_B16", "sum_qB16", "sum_EB16")):
                source[batch, k, si, ai] = parse_integer(raw[fields[field]], field=f"{stratum}.{field}", path=path)
            for pi, field in enumerate(("eligible_pair_count", "nonzero_pair_count")):
                pair_counts[batch, k, si, pi] = parse_integer(raw[fields[field]], field=f"{stratum}.{field}", path=path)
            seen_source.add(source_key)

        for ci, name in enumerate(("s_le1_pair_count", "s_le1_nonzero", "additivity_control")):
            column = columns[name]
            if column and raw.get(column, "").strip() != "":
                value = parse_integer(raw[column], field=name, path=path)
                control_key = (batch, k, name)
                if control_key in seen_controls and controls[batch, k, ci] != value:
                    raise InvalidInterface("raw_controls", f"{path}: conflicting repeated control {name} at {key}")
                controls[batch, k, ci] = value
                seen_controls.add(control_key)
        for alias, column in metadata_columns.items():
            if column and raw.get(column, "").strip() != "":
                value: Any = parse_integer(raw[column], field=alias, path=path)
                if alias in row_metadata and row_metadata[alias] != value:
                    raise InvalidInterface("metadata", f"{path}: row metadata {alias} varies")
                row_metadata[alias] = value

    if not seen_base:
        raise InvalidInterface("raw_keys", f"{path}: no data rows")
    expected_base = {(batch, k) for batch in range(batches) for k in range(n + 1)}
    if seen_base != expected_base:
        missing = sorted(expected_base - seen_base)[:10]
        extra = sorted(seen_base - expected_base)[:10]
        raise InvalidInterface(
            "raw_keys", f"{path}: every batch/K cell, including zeros, is required",
            {"first_missing": missing, "first_extra": extra},
        )
    expected_source = {(batch, k, stratum) for batch, k in expected_base for stratum in STRATA}
    if seen_source != expected_source:
        missing = sorted(expected_source - seen_source)[:10]
        raise InvalidInterface("raw_keys", f"{path}: incomplete source strata; first missing {missing}")
    observed_batches = {batch for batch, _ in seen_base}
    if observed_batches != set(range(batches)):
        raise InvalidInterface("raw_keys", f"{path}: missing batches {sorted(set(range(batches))-observed_batches)}")
    expected_per_batch = int(contract["pilot_sampling"]["paired_configurations_per_batch"])
    batch_counts = counts.sum(axis=1)
    if not np.all(batch_counts == expected_per_batch):
        raise InvalidInterface("sample_counts", f"{path}: batch counts differ from {expected_per_batch}", batch_counts.tolist())
    if np.any(counts < 0) or np.any(pair_counts < 0):
        raise InvalidInterface("raw_counts", f"{path}: negative count")
    if np.any(np.abs(baseline[:, :, 0]) > baseline[:, :, 1]) or np.any(baseline[:, :, 1] > counts):
        raise InvalidInterface("raw_counts", f"{path}: q/E moments violate q in {{-1,0,1}} and E=q^2")
    if np.any(pair_counts[..., 1] > pair_counts[..., 0]):
        raise InvalidInterface("raw_counts", f"{path}: nonzero pair count exceeds eligible count")
    empty = counts == 0
    if np.any(baseline[empty] != 0) or np.any(source[empty] != 0) or np.any(pair_counts[empty] != 0):
        raise InvalidInterface("raw_counts", f"{path}: nonzero moments in an empty K cell")

    if not np.array_equal(source[:, :, 0], source[:, :, 1] + source[:, :, 2]):
        raise InvalidInterface("source_additivity", f"{path}: total B16/qB16/EB16 differs from s2+sge3")
    if columns["s_le1_nonzero"] is None or columns["additivity_control"] is None or columns["s_le1_pair_count"] is None:
        raise InvalidInterface("raw_controls", f"{path}: all three declared control columns are required")
    expected_controls = {
        (batch, k, name)
        for batch, k in expected_base
        for name in ("s_le1_pair_count", "s_le1_nonzero", "additivity_control")
    }
    if seen_controls != expected_controls:
        raise InvalidInterface("raw_controls", f"{path}: a control field is blank or missing")
    if np.any(controls[:, :, 1] != 0) or np.any(controls[:, :, 2] != 0):
        raise InvalidInterface("raw_controls", f"{path}: exact-zero or B16 additivity control failed")
    if not np.array_equal(pair_counts[:, :, 0, 0],
                          controls[:, :, 0] + pair_counts[:, :, 1, 0] + pair_counts[:, :, 2, 0]):
        raise InvalidInterface("raw_controls", f"{path}: eligible-pair support strata do not partition total")
    if not np.array_equal(pair_counts[:, :, 0, 1], pair_counts[:, :, 1, 1] + pair_counts[:, :, 2, 1]):
        raise InvalidInterface("raw_controls", f"{path}: nonzero-pair strata do not add to total")

    sidecar, sidecar_inputs = load_sidecar_metadata(path)
    metadata: dict[str, Any] = {}
    merge_metadata(metadata, comment_metadata, "CSV comments")
    merge_metadata(metadata, sidecar, "sidecars")
    merge_metadata(metadata, row_metadata, "CSV columns")
    return {
        "N": n, "role": role, "coordinates": list(coordinates), "counts": counts,
        "baseline": baseline, "source_raw": source, "pair_counts": pair_counts,
        "controls": controls, "metadata": metadata, "row_metadata": row_metadata,
        "raw_path": path,
        "input": {"path": display_path(path), "sha256": sha256(payload), "bytes": len(payload)},
        "sidecar_inputs": sidecar_inputs,
    }


def validate_metadata(block: dict[str, Any], contract: dict[str, Any]) -> dict[str, Any]:
    n, role = block["N"], block["role"]
    row_metadata = block["row_metadata"]
    a, b = block["coordinates"]
    if a * a + b * b != n:
        raise InvalidInterface("contract", f"N{n}/{role}: a^2+b^2 does not equal N")
    expected_kernel = contract["scientific_inputs"]["joint_kernel_sha256"]
    kernel_path = ROOT / contract["scientific_inputs"]["joint_kernel_path"]
    try:
        kernel_payload = kernel_path.read_bytes()
    except OSError as error:
        raise InvalidInterface("kernel_hash", f"cannot read pinned kernel {kernel_path}: {error}") from error
    kernel = sha256(kernel_payload)
    if kernel != expected_kernel.lower():
        raise InvalidInterface("kernel_hash", f"N{n}/{role}: kernel SHA256 differs")

    window_path = Path(str(block["raw_path"]) + f".{role}.window.csv")
    window = validate_window_table(window_path, a, b)
    displacement_count = int(row_metadata.get("window_displacement_count", -1))
    anchors = int(row_metadata.get("anchors_per_configuration", -1))
    batches = int(contract["pilot_sampling"]["batches_per_size"])
    per_batch = int(contract["pilot_sampling"]["paired_configurations_per_batch"])
    if displacement_count <= 0 or anchors != contract["pilot_sampling"]["anchors_per_configuration"]:
        raise InvalidInterface("metadata", f"N{n}/{role}: changed displacement or anchor count")
    if displacement_count != window["displacement_count"]:
        raise InvalidInterface(
            "window_table", f"N{n}/{role}: CSV window count differs from the exact saved table",
            {"csv": displacement_count, "table": window["displacement_count"]},
        )
    source_denominator = int(row_metadata.get("source_denominator", -1))
    expected_denominator = 16 * anchors * n
    if source_denominator != expected_denominator:
        raise InvalidInterface(
            "source_normalization", f"N{n}/{role}: source denominator differs from 16*anchors*N",
            {"observed": source_denominator, "expected": expected_denominator},
        )
    maximum_pairs = block["counts"] * anchors * displacement_count
    if np.any(block["pair_counts"][:, :, 0, 0] > maximum_pairs):
        raise InvalidInterface("raw_counts", f"N{n}/{role}: eligible-pair count exceeds anchors*window*configurations")

    expected_streams = contract["pilot_sampling"]["rng_streams"]
    occupation = row_metadata.get("occupation_seed")
    anchor = row_metadata.get("anchor_seed")
    if occupation is None or int(occupation) != expected_streams[f"N{n}_occupation_seed"]:
        raise InvalidInterface("rng", f"N{n}/{role}: occupation seed differs")
    if anchor is None or int(anchor) != expected_streams[f"N{n}_anchor_seed"]:
        raise InvalidInterface("rng", f"N{n}/{role}: anchor seed differs")
    block["sidecar_inputs"].append({
        "path": window["path"], "sha256": window["sha256"], "bytes": window["bytes"],
        "role": "exact_complete_window_table",
    })
    return {
        "kernel_path": display_path(kernel_path), "kernel_sha256": kernel,
        "window_table_sha256": window["sha256"],
        "window_table_path": window["path"], "window_displacement_count": displacement_count,
        "exact_window_membership": window["exact_membership"],
        "incident_edge_ports_disjoint": window["incident_edge_ports_disjoint"],
        "anchors_per_configuration": anchors, "batches": batches,
        "configurations_per_batch": per_batch,
        "source_denominator": source_denominator,
        "occupation_seed": int(occupation), "anchor_seed": int(anchor),
    }


def aggregate(block: dict[str, Any], omitted_batch: int | None = None) -> dict[str, np.ndarray]:
    selector = np.ones(block["counts"].shape[0], dtype=bool)
    if omitted_batch is not None:
        selector[omitted_batch] = False
    return {
        "counts": block["counts"][selector].sum(axis=0),
        "baseline": block["baseline"][selector].sum(axis=0),
        "source_raw": block["source_raw"][selector].sum(axis=0),
    }


def weighted_jets(counts: np.ndarray, sums: np.ndarray, p: float,
                  p_reference: np.longdouble, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    """Self-normalized K-likelihood jets using scaled long-double weights."""
    dtype = np.longdouble
    k = np.arange(n + 1, dtype=dtype)
    p_ld, pref_ld = dtype(p), dtype(p_reference)
    logs = k * np.log(p_ld / pref_ld) + (dtype(n) - k) * np.log((1 - p_ld) / (1 - pref_ld))
    logs -= np.max(logs)
    weights = np.exp(logs)
    l1 = k / p_ld - (dtype(n) - k) / (1 - p_ld)
    l2 = l1 * l1 - k / (p_ld * p_ld) - (dtype(n) - k) / ((1 - p_ld) ** 2)
    c = counts.astype(dtype)
    values = sums.astype(dtype)
    z = np.sum(c * weights, dtype=dtype)
    z1 = np.sum(c * weights * l1, dtype=dtype)
    z2 = np.sum(c * weights * l2, dtype=dtype)
    if not np.isfinite(z) or z <= 0:
        raise InvalidInterface("importance_weights", "nonfinite or zero importance normalizer")
    a0 = np.sum(values * weights[:, None], axis=0, dtype=dtype)
    a1 = np.sum(values * (weights * l1)[:, None], axis=0, dtype=dtype)
    a2 = np.sum(values * (weights * l2)[:, None], axis=0, dtype=dtype)
    mean = a0 / z
    first = a1 / z - mean * (z1 / z)
    second = a2 / z - mean * (z2 / z) - 2 * (z1 / z) * first
    denominator = np.sum(c * weights * weights, dtype=dtype)
    ess = z * z / denominator
    observed = counts > 0
    diagnostics = {
        "ESS": float(ess), "ESS_fraction": float(ess / np.sum(c)),
        "retained_configurations": int(counts.sum()),
        "weight_min_observed_scaled": float(np.min(weights[observed])),
        "weight_max_observed_scaled": float(np.max(weights[observed])),
        "normalizer_scaled": float(z), "mean_l1": float(z1 / z), "mean_l2": float(z2 / z),
    }
    return np.asarray(mean, dtype=float), np.asarray(first, dtype=float), np.asarray(second, dtype=float), diagnostics


def baseline_packet(data: dict[str, np.ndarray], p: float, p_reference: np.longdouble,
                    n: int) -> dict[str, Any]:
    mean, first, second, diagnostic = weighted_jets(data["counts"], data["baseline"], p, p_reference, n)
    return {
        "q": [float(mean[0]), float(first[0]), float(second[0])],
        "E": [float(mean[1]), float(first[1]), float(second[1])],
        "importance": diagnostic,
    }


def source_packet(data: dict[str, np.ndarray], p: float, p_reference: np.longdouble, n: int,
                  divisor: float, displacement_count: int) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for si, stratum in enumerate(STRATA):
        sums = data["source_raw"][:, si, :].astype(np.longdouble) / np.longdouble(divisor)
        mean, first, _, _ = weighted_jets(data["counts"], sums, p, p_reference, n)
        result[stratum] = {
            atom: [float(mean[ai]), float(first[ai])] for ai, atom in enumerate(SOURCE_ATOMS)
        }
        factor = n / displacement_count
        result[stratum]["Cbar"] = [float(mean[0] * factor), float(first[0] * factor)]
    return result


def root_bisection(pair: dict[str, dict[str, np.ndarray]], p_reference: np.longdouble, n: int,
                   bracket: tuple[float, float], tolerance: float) -> tuple[float, dict[str, Any]]:
    def pooled(p: float) -> float:
        return 0.5 * sum(baseline_packet(pair[geometry], p, p_reference, n)["q"][0]
                         for geometry in GEOMETRIES)

    lo, hi = bracket
    flo, fhi = pooled(lo), pooled(hi)
    endpoint_values = [float(flo), float(fhi)]
    if not math.isfinite(flo) or not math.isfinite(fhi) or flo * fhi >= 0:
        raise InvalidInterface("pooled_root", "pooled M does not have opposite signs at the frozen bracket", {"bracket": bracket, "values": [flo, fhi]})
    while hi - lo > tolerance:
        mid = (lo + hi) / 2
        fmid = pooled(mid)
        if not math.isfinite(fmid):
            raise InvalidInterface("pooled_root", "nonfinite pooled M during bisection")
        if fmid == 0:
            lo = hi = mid
            break
        if flo * fmid < 0:
            hi, fhi = mid, fmid
        else:
            lo, flo = mid, fmid
    return (lo + hi) / 2, {
        "bracket": list(bracket), "pooled_M_at_bracket": endpoint_values,
        "terminal_bracket": [float(lo), float(hi)],
        "terminal_pooled_M": [float(flo), float(fhi)],
        "absolute_p_tolerance": tolerance,
    }


def vector_labels() -> list[str]:
    labels = ["root_p", "D"]
    for geometry in GEOMETRIES:
        for stratum in STRATA:
            labels.extend((f"raw_Cbar.{geometry}.{stratum}", f"raw_Cbar_p.{geometry}.{stratum}"))
    for stratum in STRATA:
        labels.extend(f"{name}.{stratum}" for name in ("jM", "jY", "jM_p", "jY_p"))
    for stratum in STRATA:
        labels.extend(f"term.{name}.{stratum}" for name in TERMS)
    labels.extend(f"J2_macro.{stratum}" for stratum in STRATA)
    labels.extend(f"T_N.{stratum}" for stratum in STRATA)
    return labels


LABELS = vector_labels()


def check_linearity(rows: dict[str, dict[str, float]], context: str) -> None:
    fields = [*TERMS, "J2_macro", "T_N", "jM", "jY", "jM_p", "jY_p"]
    for field in fields:
        total = rows["total"][field]
        parts = rows["s2"][field] + rows["sge3"][field]
        tolerance = 5e-10 * max(1.0, abs(total), abs(parts))
        if not math.isfinite(total) or abs(total - parts) > tolerance:
            raise InvalidInterface("scored_additivity", f"{context}: {field} total != s2+sge3", {"total": total, "parts": parts})


def point(blocks: dict[str, dict[str, Any]], metadata: dict[str, dict[str, Any]],
          contract: dict[str, Any], n: int, omitted_batch: int | None = None) -> tuple[np.ndarray, dict[str, Any]]:
    pair = {geometry: aggregate(blocks[geometry], omitted_batch) for geometry in GEOMETRIES}
    p_ref_num, p_ref_den = contract["pilot_sampling"]["p_reference_exact"].split("/")
    p_reference = np.longdouble(int(p_ref_num)) / np.longdouble(int(p_ref_den))
    bracket = tuple(float(value) for value in contract["single_p_root_and_p_jet_estimator"]["root_bracket"])
    p0, root_diagnostic = root_bisection(pair, p_reference, n, bracket, 1e-12)
    allowed_root = (0.5905, 0.5955)
    if not allowed_root[0] <= p0 <= allowed_root[1]:
        raise InvalidInterface("root_range", f"N{n}: root outside {allowed_root}", p0)
    baseline = {geometry: baseline_packet(pair[geometry], p0, p_reference, n) for geometry in GEOMETRIES}
    ess_floor = 0.90
    for geometry in GEOMETRIES:
        if baseline[geometry]["importance"]["ESS_fraction"] < ess_floor:
            raise InvalidInterface("importance_ESS", f"N{n}/{geometry}: ESS below 0.90 retained configurations", baseline[geometry]["importance"])
    anchors = int(contract["pilot_sampling"]["anchors_per_configuration"])
    divisor = 16 * anchors * n
    sources = {
        geometry: source_packet(pair[geometry], p0, p_reference, n, divisor,
                                metadata[geometry]["window_displacement_count"])
        for geometry in GEOMETRIES
    }
    delta = float(contract["model"]["delta_cos4"].split("/")[0]) / float(contract["model"]["delta_cos4"].split("/")[1])
    D = 0.5 * sum(baseline[geometry]["q"][1] for geometry in GEOMETRIES)
    Mpp = 0.5 * sum(baseline[geometry]["q"][2] for geometry in GEOMETRIES)
    Yp = (baseline["axis"]["E"][1] - baseline["tilted"]["E"][1]) / delta
    Ypp = (baseline["axis"]["E"][2] - baseline["tilted"]["E"][2]) / delta
    if not math.isfinite(D) or D <= 0:
        raise InvalidInterface("root_slope", f"N{n}: nonpositive D", D)
    R = Yp / D
    prefactor = n ** (13 / 8) / 2
    scored: dict[str, dict[str, float]] = {}
    for stratum in STRATA:
        geometry_rows = {}
        for geometry in GEOMETRIES:
            q, qp, _ = baseline[geometry]["q"]
            e, ep, _ = baseline[geometry]["E"]
            h, hp = sources[geometry][stratum]["H"]
            qh, qhp = sources[geometry][stratum]["qH"]
            eh, ehp = sources[geometry][stratum]["EH"]
            geometry_rows[geometry] = {
                "j_q": qh - q * h,
                "j_q_p": qhp - qp * h - q * hp,
                "j_E": eh - e * h,
                "j_E_p": ehp - ep * h - e * hp,
            }
        jM = 0.5 * (geometry_rows["axis"]["j_q"] + geometry_rows["tilted"]["j_q"])
        jM_p = 0.5 * (geometry_rows["axis"]["j_q_p"] + geometry_rows["tilted"]["j_q_p"])
        jY = (geometry_rows["axis"]["j_E"] - geometry_rows["tilted"]["j_E"]) / delta
        jY_p = (geometry_rows["axis"]["j_E_p"] - geometry_rows["tilted"]["j_E_p"]) / delta
        terms = {
            "direct_centered": prefactor * jY_p / D,
            "root_motion": -prefactor * Ypp * jM / (D * D),
            "slope_source": -prefactor * R * jM_p / D,
            "slope_root": prefactor * R * Mpp * jM / (D * D),
        }
        j2 = math.fsum(terms.values())
        scored[stratum] = {
            "jM": jM, "jY": jY, "jM_p": jM_p, "jY_p": jY_p,
            **terms, "J2_macro": j2, "T_N": n * n * j2,
            "root_joint_tangent": -jM / D,
        }
    check_linearity(scored, f"N{n}/omit={omitted_batch}")

    vector_values: dict[str, float] = {"root_p": p0, "D": D}
    for geometry in GEOMETRIES:
        for stratum in STRATA:
            vector_values[f"raw_Cbar.{geometry}.{stratum}"] = sources[geometry][stratum]["Cbar"][0]
            vector_values[f"raw_Cbar_p.{geometry}.{stratum}"] = sources[geometry][stratum]["Cbar"][1]
    for stratum in STRATA:
        for name in ("jM", "jY", "jM_p", "jY_p"):
            vector_values[f"{name}.{stratum}"] = scored[stratum][name]
    for stratum in STRATA:
        for name in TERMS:
            vector_values[f"term.{name}.{stratum}"] = scored[stratum][name]
    for stratum in STRATA:
        vector_values[f"J2_macro.{stratum}"] = scored[stratum]["J2_macro"]
    for stratum in STRATA:
        vector_values[f"T_N.{stratum}"] = scored[stratum]["T_N"]
    vector = np.asarray([vector_values[label] for label in LABELS], dtype=float)
    if not np.isfinite(vector).all():
        raise InvalidInterface("finite_scores", f"N{n}: nonfinite score vector")
    diagnostic = {
        "root": root_diagnostic,
        "root_p": p0, "baseline": baseline,
        "baseline_derived": {"D": D, "M_pp": Mpp, "Y_p": Yp, "Y_pp": Ypp, "R": R, "A_N": prefactor},
        "source_packets": sources, "source_responses": scored,
        "omitted_batch": omitted_batch,
    }
    return vector, diagnostic


def analyze_size(blocks: dict[str, dict[str, Any]], metadata: dict[str, dict[str, Any]],
                 contract: dict[str, Any], n: int) -> dict[str, Any]:
    central, diagnostic = point(blocks, metadata, contract, n)
    deletes, delete_diagnostics = [], []
    batches = int(contract["pilot_sampling"]["batches_per_size"])
    for omitted in range(batches):
        vector, detail = point(blocks, metadata, contract, n, omitted)
        deletes.append(vector)
        delete_diagnostics.append({
            "omitted_batch": omitted,
            "root_p": detail["root_p"],
            "D": detail["baseline_derived"]["D"],
            "root": detail["root"],
            "importance": {
                geometry: detail["baseline"][geometry]["importance"]
                for geometry in GEOMETRIES
            },
        })
    loo = np.asarray(deletes, dtype=float)
    deviations = loo - loo.mean(axis=0)
    covariance = (batches - 1) / batches * deviations.T @ deviations
    errors = np.sqrt(np.maximum(0.0, np.diag(covariance)))
    estimates = {
        label: {"value": float(value), "se": float(se),
                "z": float(value / se) if se > 0 else None}
        for label, value, se in zip(LABELS, central, errors)
    }
    all_ess_fractions = [
        diagnostic["baseline"][geometry]["importance"]["ESS_fraction"]
        for geometry in GEOMETRIES
    ] + [
        item["importance"][geometry]["ESS_fraction"]
        for item in delete_diagnostics for geometry in GEOMETRIES
    ]
    return {
        "status": "scoreable", "labels": LABELS, "estimates": estimates,
        "central_vector": central.tolist(), "covariance": covariance.tolist(),
        "delete_one_vectors": loo.tolist(), "diagnostics": diagnostic,
        "delete_one_diagnostics": delete_diagnostics,
        "minimum_central_or_delete_one_ESS_fraction": min(all_ess_fractions),
        "delete_one_root_range": [
            min(item["root_p"] for item in delete_diagnostics),
            max(item["root_p"] for item in delete_diagnostics),
        ],
    }


def simultaneous_interval(estimate: dict[str, Any], critical: float) -> list[float]:
    return [estimate["value"] - critical * estimate["se"],
            estimate["value"] + critical * estimate["se"]]


def production_projection(value: float, se: float, n_pilot: int, critical: float) -> int | None:
    if value == 0 or not math.isfinite(value) or not math.isfinite(se):
        return None
    raw = n_pilot * (critical * se / (0.25 * abs(value))) ** 2
    return int(math.ceil(raw / 100.0) * 100)


def downstream(by_n: dict[str, dict[str, Any]], contract: dict[str, Any], critical: float) -> dict[str, Any]:
    t100 = by_n["100"]["estimates"]["T_N.total"]
    t400 = by_n["400"]["estimates"]["T_N.total"]
    i100, i400 = simultaneous_interval(t100, critical), simultaneous_interval(t400, critical)
    signs = []
    for interval in (i100, i400):
        signs.append(1 if interval[0] > 0 else (-1 if interval[1] < 0 else 0))
    status = "evaluated" if signs[0] != 0 and signs[0] == signs[1] else "not_evaluated_due_to_unresolved_or_sign_inconsistent_thermal_tail"
    result: dict[str, Any] = {
        "field_ratio_status": status,
        "simultaneous_familywise_95_intervals": {"T100": i100, "T400": i400},
        "critical_value_t_0.9875_df99": critical,
    }
    if status != "evaluated":
        result.update({"D17": None, "D21": None, "contrast_covariance": None})
        return result
    r17 = 2 ** (-5 / 4)
    r21 = 2 ** (-13 / 4)
    var100 = t100["se"] ** 2
    var400 = t400["se"] ** 2
    d17 = t400["value"] - r17 * t100["value"]
    d21 = t400["value"] - r21 * t100["value"]
    covariance = np.asarray([
        [var400 + r17 * r17 * var100, var400 + r17 * r21 * var100],
        [var400 + r17 * r21 * var100, var400 + r21 * r21 * var100],
    ])
    result.update({
        "ratios": {"r17": r17, "r21": r21},
        "D17": {"value": d17, "se": math.sqrt(covariance[0, 0])},
        "D21": {"value": d21, "se": math.sqrt(covariance[1, 1])},
        "contrast_covariance": covariance.tolist(),
        "boundary": contract["conditional_downstream_models"]["interpretation"],
    })
    return result


def report_markdown(result: dict[str, Any]) -> str:
    lines = ["# Canonical macro-window joint-U pilot", "", f"Status: **{result['status']}**.", ""]
    if result["status"] == "INVALID_INTERFACE":
        failure = result["failure"]
        lines.extend([f"Gate: `{failure['gate']}`", "", failure["reason"], "",
                      "No field contrast or production projection was evaluated.", ""])
        return "\n".join(lines)
    lines.extend([
        "The same fixed canonical Kreg, Euclidean 1/4--2/5 window and complete moving-root U functional are used at both sizes.",
        "Anchors and displacement pairs are within-configuration readouts; the 100 paired batches are the inference units.", "",
        "| N | root p | D | ESS min fraction | T_N total | simultaneous 95% interval | projected configurations |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for n in ("100", "400"):
        row = result["by_N"][n]
        est = row["estimates"]
        interval = result["downstream"]["simultaneous_familywise_95_intervals"][f"T{n}"]
        projection = result["production_projection"][n]
        lines.append(
            f"| {n} | {est['root_p']['value']:.12g} | {est['D']['value']:.8g} | "
            f"{row['minimum_central_or_delete_one_ESS_fraction']:.6f} | {est['T_N.total']['value']:+.8g} ± {est['T_N.total']['se']:.3g} | "
            f"[{interval[0]:+.8g}, {interval[1]:+.8g}] | {projection if projection is not None else 'not available'} |"
        )
    lines.extend(["", "## Complete original-U decomposition", "",
                  "| N | support | direct | root motion | source slope | root slope | J2_macro | T_N |",
                  "|---:|---|---:|---:|---:|---:|---:|---:|"])
    for n in ("100", "400"):
        estimates = result["by_N"][n]["estimates"]
        for stratum in STRATA:
            values = [estimates[f"term.{term}.{stratum}"]["value"] for term in TERMS]
            lines.append(
                f"| {n} | {stratum} | " + " | ".join(f"{value:+.8g}" for value in values)
                + f" | {estimates[f'J2_macro.{stratum}']['value']:+.8g} | {estimates[f'T_N.{stratum}']['value']:+.8g} |"
            )
    lines.extend(["", "## Conditional fixed-power contrasts", "",
                  f"Status: `{result['downstream']['field_ratio_status']}`.", ""])
    if result["downstream"]["field_ratio_status"] == "evaluated":
        for key in ("D17", "D21"):
            row = result["downstream"][key]
            lines.append(f"- {key}: {row['value']:+.8g} ± {row['se']:.3g}")
        lines.append("")
    lines.extend([
        "D17/D21 are conditional pilot contrasts only. They do not estimate an exponent, accept/reject a field model or authorize production.",
        "The s2 and sge3 rows are correlated additive coordinates of the same source, not independent evidence.",
        "Raw Cbar diagnostics, all delete-one vectors and the complete covariance are stored in the JSON.", "",
    ])
    return "\n".join(lines)


def analyze(args: argparse.Namespace) -> dict[str, Any]:
    started = time.perf_counter()
    try:
        contract_bytes = args.contract.read_bytes()
        contract = json.loads(contract_bytes)
    except (OSError, json.JSONDecodeError) as error:
        raise InvalidInterface("contract", f"cannot read contract {args.contract}: {error}") from error
    if contract.get("schema") != "matching-one.regular-pair-macro-joint-u.contract.v1":
        raise InvalidInterface("contract", "unexpected macro-joint contract schema")
    exact_reference = contract["pilot_sampling"]["p_reference_exact"].split("/")
    if (len(exact_reference) != 2
            or int(exact_reference[0]) != int(contract["pilot_sampling"]["p_reference_uint64_threshold"])
            or int(exact_reference[1]) != 2 ** 64):
        raise InvalidInterface("contract", "p_reference exact fraction and uint64 threshold disagree")
    model_by_n = {int(item["N"]): item for item in contract["model"]["geometries"]}
    paths = {
        100: {"axis": args.axis_100, "tilted": args.tilted_100},
        400: {"axis": args.axis_400, "tilted": args.tilted_400},
    }
    blocks: dict[int, dict[str, dict[str, Any]]] = {}
    metadata: dict[int, dict[str, dict[str, Any]]] = {}
    inputs: list[dict[str, Any]] = []
    for n in SIZES:
        blocks[n], metadata[n] = {}, {}
        for role in GEOMETRIES:
            coordinates = tuple(model_by_n[n][role])
            block = load_csv(paths[n][role], n=n, role=role, coordinates=coordinates, contract=contract)
            blocks[n][role] = block
            metadata[n][role] = validate_metadata(block, contract)
            inputs.extend([block["input"], *block["sidecar_inputs"]])
        if not np.array_equal(blocks[n]["axis"]["counts"], blocks[n]["tilted"]["counts"]):
            raise InvalidInterface("geometry_pairing", f"N{n}: paired geometries have different batch-by-K counts")
        # The same sampled anchor indices do not imply identical window tables;
        # each geometry retains its own pinned displacement count and hash.

    kernel_path = ROOT / contract["scientific_inputs"]["joint_kernel_path"]
    kernel_payload = kernel_path.read_bytes()
    inputs.append({
        "path": display_path(kernel_path), "sha256": sha256(kernel_payload),
        "bytes": len(kernel_payload), "role": "pinned_joint_kernel",
    })
    unique_inputs: list[dict[str, Any]] = []
    seen_inputs: set[tuple[str, str]] = set()
    for item in inputs:
        key = (item["path"], item["sha256"])
        if key not in seen_inputs:
            seen_inputs.add(key)
            unique_inputs.append(item)

    by_n = {str(n): analyze_size(blocks[n], metadata[n], contract, n) for n in SIZES}
    batches = int(contract["pilot_sampling"]["batches_per_size"])
    critical = float(student_t.ppf(0.9875, batches - 1))
    downstream_result = downstream(by_n, contract, critical)
    d_intervals = {
        str(n): simultaneous_interval(by_n[str(n)]["estimates"]["D"], critical)
        for n in SIZES
    }
    if any(interval[0] <= 0 for interval in d_intervals.values()):
        raise InvalidInterface("positive_D_interval", "D does not have a positive simultaneous pilot interval", d_intervals)

    total_labels = [f"N{n}.{label}" for n in SIZES for label in LABELS]
    dimension = len(LABELS)
    joint_covariance = np.zeros((2 * dimension, 2 * dimension), dtype=float)
    joint_covariance[:dimension, :dimension] = np.asarray(by_n["100"]["covariance"])
    joint_covariance[dimension:, dimension:] = np.asarray(by_n["400"]["covariance"])
    projection_critical = float(student_t.ppf(0.995, batches - 1))
    n_pilot = int(contract["pilot_sampling"]["paired_configurations_per_size"])
    projections = {
        str(n): production_projection(by_n[str(n)]["estimates"]["T_N.total"]["value"],
                                      by_n[str(n)]["estimates"]["T_N.total"]["se"],
                                      n_pilot, projection_critical)
        for n in SIZES
    }
    ceiling = 2_000_000
    signs = [by_n[str(n)]["estimates"]["T_N.total"]["value"] for n in SIZES]
    eligible = (all(value != 0 for value in signs) and signs[0] * signs[1] > 0
                and all(value is not None and value <= ceiling for value in projections.values()))
    upgrade_status = (
        "eligible_to_write_new_frozen_contract_no_action_authorized"
        if eligible else "stop_fixed_window_field_ratio_route_under_declared_gate"
    )
    result = {
        "schema": "matching-one.regular-pair-macro-joint-u.v1",
        "status": "completed_valid_pilot_score",
        "contract": {"path": display_path(args.contract), "sha256": sha256(contract_bytes), "content": contract},
        "inputs": unique_inputs, "input_metadata": metadata,
        "labels": LABELS, "by_N": by_n,
        "between_size_covariance": {
            "rule": "block_diagonal_disjoint_size_streams", "labels": total_labels,
            "covariance": joint_covariance.tolist(),
        },
        "downstream": downstream_result,
        "D_simultaneous_intervals": d_intervals,
        "production_projection": projections,
        "production_projection_critical_value_t_0.995_df99": projection_critical,
        "production_upgrade_contract_eligible": eligible,
        "production_upgrade_status": upgrade_status,
        "production_action_authorized": False,
        "validity_gates": {
            "kernel_and_window_hashes": "passed", "exact_window_membership": "passed",
            "disjoint_edge_port_groups": "passed", "geometry_pairing": "passed",
            "source_denominator": "passed", "s_le1_zero": "passed",
            "raw_and_scored_additivity": "passed", "all_roots_bracketed": "passed",
            "root_range": "passed", "positive_D_simultaneous_interval": "passed",
            "ESS_fraction_at_least_0.90": "passed", "all_batches_retained": "passed",
        },
        "definitions": {
            "importance_reweighting": "(p/p_reference)^K*((1-p)/(1-p_reference))^(N-K), self-normalized per geometry",
            "source_unit": "H=sum_B16/(16*anchors*N); anchors=16",
            "raw_Cbar": "Cbar=H*N/window_displacement_count",
            "uncertainty": "100 aligned paired-batch delete-one recomputations per size; complete covariance; sizes block diagonal",
            "boundary": contract["complete_original_U_estimator"]["observer_boundary"],
        },
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "machine": platform.machine()},
        "command": [sys.executable, *sys.argv],
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter() - started,
        "new_samples": 0, "sampler_runs": 0, "scientific_tests": [],
    }
    return result


def invalid_result(args: argparse.Namespace, error: InvalidInterface, started: float) -> dict[str, Any]:
    inputs = []
    for path in (args.axis_100, args.tilted_100, args.axis_400, args.tilted_400):
        if path.exists():
            payload = path.read_bytes()
            inputs.append({"path": display_path(path), "sha256": sha256(payload), "bytes": len(payload)})
    contract_payload = args.contract.read_bytes() if args.contract.exists() else b""
    return {
        "schema": "matching-one.regular-pair-macro-joint-u.v1",
        "status": "INVALID_INTERFACE",
        "failure": {"gate": error.gate, "reason": error.reason, "diagnostic": jsonable(error.diagnostic)},
        "contract": {"path": display_path(args.contract), "sha256": sha256(contract_payload) if contract_payload else None},
        "inputs": inputs, "field_ratio_status": "not_evaluated_due_to_invalid_interface",
        "production_action_authorized": False, "command": [sys.executable, *sys.argv],
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter() - started,
        "new_samples": 0, "sampler_runs": 0, "scientific_tests": [],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--axis-100", type=Path, required=True)
    parser.add_argument("--tilted-100", type=Path, required=True)
    parser.add_argument("--axis-400", type=Path, required=True)
    parser.add_argument("--tilted-400", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.output_json.exists() or args.output_md.exists():
        raise SystemExit("refusing to overwrite an existing output")
    started = time.perf_counter()
    try:
        result = analyze(args)
        exit_code = 0
    except InvalidInterface as error:
        result = invalid_result(args, error, started)
        exit_code = 2
    report = report_markdown(result)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    with args.output_json.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(jsonable(result), indent=2, allow_nan=False) + "\n")
    with args.output_md.open("x", encoding="utf-8") as handle:
        handle.write(report)
    print(json.dumps({"status": result["status"], "output_json": str(args.output_json),
                      "output_md": str(args.output_md)}, indent=2))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
