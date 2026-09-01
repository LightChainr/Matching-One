#!/usr/bin/env python3
"""Exact existing-fibre score for the frozen P537 diagonal one-defect gate."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
N = 25
RAW_FIELDS = (
    "count", "sum_q0", "sum_E0", "sum_a16_0", "sum_q0_a16_0", "sum_E0_a16_0",
    "sum_q1", "sum_E1", "sum_a16_1", "sum_q1_a16_1", "sum_E1_a16_1",
)
ORBITS = {
    "nn_other": {(-1, 0), (0, -1), (0, 1)},
    "diag1": {(-1, -1), (-1, 1), (1, -1), (1, 1)},
    "axial2": {(-2, 0), (0, -2), (2, 0), (0, 2)},
    "knight_a": {(-2, -1), (1, -2), (2, 1), (-1, 2)},
    "knight_b": {(-2, 1), (-1, -2), (2, -1), (1, 2)},
}
FAR = {
    "axis": {(-2, -2), (-2, 2), (2, -2), (2, 2)},
    "tilted": {(-3, 0), (0, -3), (0, 3), (3, 0)},
}


def load_schur(path: Path):
    spec = importlib.util.spec_from_file_location("p537_schur", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def source_component(row: dict[str, str], geometry: str) -> str:
    displacement = (int(row["y_dx"]), int(row["y_dy"]))
    for name, members in ORBITS.items():
        if displacement in members:
            return name
    if displacement in FAR[geometry]:
        return "far"
    raise ValueError(f"unmapped {geometry} displacement {displacement}")


def selected(row: dict[str, str]) -> bool:
    return (
        row["alternating_four_arm"] == "1"
        and row["rank0"] != row["rank1"]
        and row["source_absent"] == "0"
        and row["bell0"] != row["bell1"]
    )


def add(table: dict, key: tuple, row: dict[str, str]) -> None:
    packet = table.setdefault(key, [0] * len(RAW_FIELDS))
    for i, field in enumerate(RAW_FIELDS):
        packet[i] += int(row[field])


def packet_dict(values: list[int]) -> dict[str, F]:
    return {name: F(value) for name, value in zip(RAW_FIELDS, values)}


def midpoint(interval) -> F:
    return (interval.lo + interval.hi) / 2


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--axis", type=Path, required=True)
    ap.add_argument("--tilted", type=Path, required=True)
    ap.add_argument("--aggregates", type=Path, required=True)
    ap.add_argument("--baseline-axis", type=Path, required=True)
    ap.add_argument("--baseline-tilted", type=Path, required=True)
    ap.add_argument("--baseline-root", type=Path, required=True)
    ap.add_argument("--kernel", type=Path, required=True)
    ap.add_argument("--schur-module", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--require-kernel-change", action="store_true")
    args = ap.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")

    schur = load_schur(args.schur_module)
    p = schur.read_root(args.baseline_root)
    coefficients = {
        "axis": schur.read_baseline(args.baseline_axis, N),
        "tilted": schur.read_baseline(args.baseline_tilted, N),
    }
    baseline = {g: schur.baseline_packet(coefficients[g], N, p) for g in ("axis", "tilted")}
    delta = F(1152, 625)
    mt = (baseline["axis"]["q_t"] + baseline["tilted"]["q_t"]) / 2
    yt = (baseline["axis"]["e_t"] - baseline["tilted"]["e_t"]) / delta
    root_ratio = yt / mt
    c = {"axis": F(1, 1) / delta, "tilted": -F(1, 1) / delta}
    mu_h = {
        g: 2 * c[g] * baseline[g]["e"] - root_ratio * baseline[g]["q"]
        for g in ("axis", "tilted")
    }
    aggregate_rows = schur.read_aggregates(args.aggregates, N)
    source_packets, beta = schur.source_global_packets(
        aggregate_rows, baseline, N, p, a_raw_denominator=16
    )
    kernel = {}
    with args.kernel.open(newline="") as handle:
        for row in csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"):
            kernel[int(row.get("key", row.get("packed_key")))] = int(row["g16"])

    def accept(row: dict[str, str]) -> bool:
        if not selected(row):
            return False
        if not args.require_kernel_change:
            return True
        return kernel.get(int(row["bell0"]), 0) != kernel.get(int(row["bell1"]), 0)

    def score_rows(geometry: str, component: str, by_k: dict[int, list[int]]):
        signed = [schur.Interval.of(0), schur.Interval.of(0)]
        positive = [schur.Interval.of(0), schur.Interval.of(0)]
        for k, values in by_k.items():
            row = packet_dict(values)
            weight = schur.offsite_weight(k, N, p)
            mu_a = source_packets[(geometry, component)]["mu_a"]
            beta_component = beta[component]
            s_minus = k - (N - 1) * p
            for i in (0, 1):
                suffix = str(i)
                wi = (1 - p) if i == 0 else p
                ui = i - p
                si = s_minus + ui
                bi = ui * si - p * (1 - p)
                count = row["count"]
                sum_q = row[f"sum_q{suffix}"]
                sum_e = row[f"sum_E{suffix}"]
                sum_a = row[f"sum_a16_{suffix}"] / (N * 16)
                sum_qa = row[f"sum_q{suffix}_a16_{suffix}"] / (N * 16)
                sum_ea = row[f"sum_E{suffix}_a16_{suffix}"] / (N * 16)
                sum_h = 2 * c[geometry] * sum_e - root_ratio * sum_q - mu_h[geometry] * count
                sum_ha = (
                    2 * c[geometry] * (sum_ea - mu_a * sum_e)
                    - root_ratio * (sum_qa - mu_a * sum_q)
                    - mu_h[geometry] * (sum_a - mu_a * count)
                )
                # x4 restores the fixed-z C4 orbit; /2 is the geometry pool.
                signed[i] += weight * schur.eq10_state_term(
                    wi, ui, bi, sum_h, sum_ha, beta_component, 4
                ) / 2
                positive[i] += weight * wi * count * 2
        return {"P": positive, "S": signed, "total": signed[0] + signed[1]}

    sources = {"axis": args.axis, "tilted": args.tilted}
    coarse: dict[tuple, list[int]] = {}
    selected_rows = selected_mass = 0
    for geometry, path in sources.items():
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if not accept(row):
                    continue
                if not (row["outer_occupied_join"] == "1" and row["outer_vacant_join"] == "1"):
                    raise ValueError("rank-changing alternating edge lacks J_B=J_W=1")
                selected_rows += 1
                selected_mass += int(row["count"])
                component = source_component(row, geometry)
                base = (
                    geometry, f"{row['rank0']}->{row['rank1']}", component,
                    int(row["local_source_contact_mask"]), int(row["collar_corner_mask"]),
                )
                add(coarse, base + (int(row["k_minus"]),), row)

    coarse_by_base: dict[tuple, dict[int, list[int]]] = defaultdict(dict)
    for key, values in coarse.items():
        coarse_by_base[key[:-1]][key[-1]] = values
    coarse_scores = []
    first_coarse = None
    matrix = defaultdict(lambda: schur.Interval.of(0))
    matrix_by_contact = defaultdict(lambda: schur.Interval.of(0))
    corner_sums = defaultdict(lambda: schur.Interval.of(0))
    for base in sorted(coarse_by_base):
        geometry, transition, component, contact, corner = base
        scored = score_rows(geometry, component, coarse_by_base[base])
        record = {
            "geometry": geometry, "rank_transition": transition,
            "source_component": component, "contact_mask": contact,
            "corner_mask": corner, "P0": schur.interval_record(scored["P"][0]),
            "P1": schur.interval_record(scored["P"][1]),
            "S0": schur.interval_record(scored["S"][0]),
            "S1": schur.interval_record(scored["S"][1]),
            "weight": schur.interval_record(scored["total"]),
        }
        coarse_scores.append(record)
        matrix[(transition, component)] += scored["total"]
        matrix_by_contact[(contact, transition, component)] += scored["total"]
        corner_sums[corner] += scored["total"]
        if first_coarse is None and record["weight"]["excludes_zero"]:
            first_coarse = (base, record)
    if first_coarse is None:
        status = "NO_DIAGONAL_EDGE_AT_COARSE_RESOLUTION"
        first_bell_record = None
        bell_class_count = 0
    else:
        target_base, _ = first_coarse
        detailed: dict[tuple, list[int]] = {}
        geometry, transition, component, contact, corner = target_base
        for path_geometry, path in sources.items():
            if path_geometry != geometry:
                continue
            with path.open(newline="") as handle:
                for row in csv.DictReader(handle):
                    if not accept(row):
                        continue
                    base = (
                        geometry, f"{row['rank0']}->{row['rank1']}",
                        source_component(row, geometry), int(row["local_source_contact_mask"]),
                        int(row["collar_corner_mask"]),
                    )
                    if base != target_base:
                        continue
                    key = (int(row["bell0"]), int(row["bell1"]), int(row["k_minus"]))
                    add(detailed, key, row)
        by_bell: dict[tuple[int, int], dict[int, list[int]]] = defaultdict(dict)
        for key, values in detailed.items():
            by_bell[key[:2]][key[2]] = values
        bell_class_count = len(by_bell)
        first_bell_record = None
        for bell_pair in sorted(by_bell):
            scored = score_rows(geometry, component, by_bell[bell_pair])
            record = {
                "geometry": geometry, "rank_transition": transition,
                "source_component": component, "contact_mask": contact,
                "corner_mask": corner, "bell0": bell_pair[0], "bell1": bell_pair[1],
                "P0": schur.interval_record(scored["P"][0]),
                "P1": schur.interval_record(scored["P"][1]),
                "S0": schur.interval_record(scored["S"][0]),
                "S1": schur.interval_record(scored["S"][1]),
                "weight": schur.interval_record(scored["total"]),
            }
            if record["weight"]["excludes_zero"]:
                first_bell_record = record
                break
        status = (
            "TWO_INDEPENDENT_DEFECT_GAIN_REJECTED"
            if first_bell_record is not None else "COARSE_NONZERO_BUT_NO_BELL_CERTIFICATE"
        )

    if first_bell_record is not None:
        # The canonical kernel TSV is sparse; omitted Bell keys have g16=0.
        first_bell_record["g16_0"] = kernel.get(first_bell_record["bell0"], 0)
        first_bell_record["g16_1"] = kernel.get(first_bell_record["bell1"], 0)

    transitions = sorted({key[0] for key in matrix})
    components = sorted({key[1] for key in matrix})
    row_sums = {t: sum((matrix[(t, a)] for a in components), schur.Interval.of(0)) for t in transitions}
    column_sums = {a: sum((matrix[(t, a)] for t in transitions), schur.Interval.of(0)) for a in components}
    contact_masks = sorted({key[0] for key in matrix_by_contact})
    contact_sums = {
        contact: sum(
            (matrix_by_contact[(contact, t, a)] for t in transitions for a in components),
            schur.Interval.of(0),
        )
        for contact in contact_masks
    }
    contact0 = contact_sums.get(0, schur.Interval.of(0))
    contact0_decision = (
        "RADIUS_ONE_CONTACT_ONLY_CLOSURE_REJECTED"
        if contact0.lo > 0 or contact0.hi < 0 else "CONTACT_ZERO_RESIDUAL_UNRESOLVED"
    )
    payload = {
        "schema": "matching-one/p537-one-defect-diagonal-edge/v1",
        "status": status,
        "definition": (
            "alternating site flip changes rank, canonical source Bell state, and g16"
            if args.require_kernel_change else
            "alternating site flip changes both rank and canonical source Bell state"
        ) + "; full pooled-root Schur coefficients fixed globally",
        "require_kernel_change": args.require_kernel_change,
        "selected_raw_row_classes": selected_rows,
        "selected_physical_pair_fibres": selected_mass,
        "coarse_class_count": len(coarse_scores),
        "first_nonzero_coarse_class": first_coarse[1] if first_coarse else None,
        "bell_classes_inside_selected_coarse_class": bell_class_count,
        "first_nonzero_bell_transition": first_bell_record,
        "signed_mass_matrix": {
            "row_order": transitions, "column_order": components,
            "cells": [[schur.interval_record(matrix[(t, a)]) for a in components] for t in transitions],
            "S_times_1": {t: schur.interval_record(row_sums[t]) for t in transitions},
            "1T_times_S": {a: schur.interval_record(column_sums[a]) for a in components},
        },
        "contact_decomposition": {
            "bit_semantics": "bit 0/1: source cut contacts the two local occupied thermal arms; mask 0 is the radius-one no-contact residual",
            "sums": {str(contact): schur.interval_record(contact_sums[contact]) for contact in contact_masks},
            "contact0_matrix": [
                [schur.interval_record(matrix_by_contact[(0, t, a)]) for a in components]
                for t in transitions
            ],
            "decision": contact0_decision,
        },
        "corner_word_sums": {str(corner): schur.interval_record(value) for corner, value in sorted(corner_sums.items())},
        "global": {"root": schur.interval_record(p), "M_t": schur.interval_record(mt), "R": schur.interval_record(root_ratio)},
        "logic": "a nonzero aggregate Bell-transition class proves at least one physical diagonal edge has nonzero full Schur signed weight",
        "boundary": "existing exact N25 fibres; transition-class existence certificate, not a retained single-background mask",
        "inputs": {
            name: {"path": str(path), "sha256": schur.sha256(path)}
            for name, path in {
                "axis": args.axis, "tilted": args.tilted, "aggregates": args.aggregates,
                "baseline_axis": args.baseline_axis, "baseline_tilted": args.baseline_tilted,
                "baseline_root": args.baseline_root, "kernel": args.kernel,
                "schur_module": args.schur_module,
            }.items()
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": status, "output": str(args.output), "selected": selected_rows}, sort_keys=True))


if __name__ == "__main__":
    main()
