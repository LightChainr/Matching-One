#!/usr/bin/env python3
"""Compress exact site-flip fibres to the frozen P537 Schur scorer schema.

The retained landing block is deliberately named ``near_block``: it requires
an alternating local occupation mask and no contact from either source cut to
the two occupied thermal landing components.  The CSV still records the
occupied degree-branch flag and the occupied/matching landing identifications
in tau; it does not promote this finite label to an asymptotic arm event.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


VALUE_FIELDS = (
    "count",
    "sum_q0",
    "sum_E0",
    "sum_a16_0",
    "sum_q0_a16_0",
    "sum_E0_a16_0",
    "sum_q1",
    "sum_E1",
    "sum_a16_1",
    "sum_q1_a16_1",
    "sum_E1_a16_1",
)
OUT_FIELDS = (
    "count",
    "sum_q0",
    "sum_q1",
    "sum_e0",
    "sum_e1",
    "sum_a0",
    "sum_a1",
    "sum_q0a0",
    "sum_q1a1",
    "sum_e0a0",
    "sum_e1a1",
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
REPRESENTATIVE = {
    "axis": {
        "nn_other": (-1, 0), "diag1": (-1, -1), "axial2": (-2, 0),
        "knight_a": (-2, -1), "knight_b": (-2, 1), "far": (-2, -2),
    },
    "tilted": {
        "nn_other": (-1, 0), "diag1": (-1, -1), "axial2": (-2, 0),
        "knight_a": (-2, -1), "knight_b": (-2, 1), "far": (-3, 0),
    },
}


def source_component(row: dict[str, str], geometry: str) -> str:
    displacement = (int(row["y_dx"]), int(row["y_dy"]))
    for name, members in ORBITS.items():
        if displacement in members:
            return name
    if displacement in FAR[geometry]:
        return "far"
    raise ValueError(f"unmapped {geometry} displacement {displacement}")


def alpha(row: dict[str, str], geometry: str) -> str:
    component = source_component(row, geometry)
    return f"{component}:{'absent' if row['source_absent'] == '1' else 'present'}"


def tau(row: dict[str, str]) -> str:
    return "near_block:" + json.dumps(
        [int(row["rank0"]), int(row["rank1"]), int(row["arm_mask"]),
         degree_branch(row),
         int(row["occupied_landing_id0"]), int(row["occupied_landing_id1"]),
         int(row["vacant_separator_id0"]), int(row["vacant_separator_id1"])],
        separators=(",", ":"),
    )


def retained(row: dict[str, str]) -> bool:
    return row["alternating_four_arm"] == "1" and row["extra_source_port_contact"] == "0"


def degree_branch(row: dict[str, str]) -> int:
    # The first local scratch traversal used the broader legacy header.  Its
    # value already scanned the whole occupied landing component, including
    # the four port vertices; only the label changes.
    field = "occupied_component_degree_branch" if "occupied_component_degree_branch" in row else "off_port_occupied_branch"
    return int(row[field])


def add(target: dict[tuple[str, str, str, int], list[int]], key, row) -> None:
    packet = target.setdefault(key, [0] * len(VALUE_FIELDS))
    for i, field in enumerate(VALUE_FIELDS):
        packet[i] += int(row[field])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--axis", type=Path, required=True)
    parser.add_argument("--tilted", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sources = {"axis": args.axis, "tilted": args.tilted}

    landing: dict[tuple[str, str, str, str, int], list[int]] = {}
    cross_tab: dict[str, dict[tuple[int, int, int, int, int], int]] = {
        geometry: defaultdict(int) for geometry in sources
    }
    selected_alphas: set[str] = set()
    alpha_component: dict[str, str] = {}
    selected_taus: set[str] = set()
    for geometry, path in sources.items():
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if row["alternating_four_arm"] == "1":
                    cross_tab[geometry][(
                        int(row["occupied_landing_id1"] == "1"),
                        int(row["vacant_separator_id1"] == "1"),
                        degree_branch(row),
                        int(row["extra_source_port_contact"]),
                        int(row["source_absent"]),
                    )] += int(row["count"])
                if not retained(row):
                    continue
                aa = alpha(row, geometry)
                component = source_component(row, geometry)
                tt, k = tau(row), int(row["k_minus"])
                selected_alphas.add(aa)
                if aa in alpha_component and alpha_component[aa] != component:
                    raise ValueError("alpha maps to more than one fixed source component")
                alpha_component[aa] = component
                selected_taus.add(tt)
                add(landing, (geometry, tt, aa, component, k), row)

    selected_components = set(alpha_component.values())
    global_rows: dict[tuple[str, str, str, str, int], list[int]] = {}
    for geometry, path in sources.items():
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                component = source_component(row, geometry)
                displacement = (int(row["y_dx"]), int(row["y_dy"]))
                if component in selected_components and displacement == REPRESENTATIVE[geometry][component]:
                    add(global_rows, (geometry, "__GLOBAL__", "__SOURCE__", component,
                                      int(row["k_minus"])), row)

    # The scorer requires a literal common rectangle.  Missing cells are exact
    # zeros, not absent observations.  A single k=0 zero row is sufficient.
    zero = [0] * len(VALUE_FIELDS)
    for geometry in sources:
        for component in selected_components:
            if not any(key[0] == geometry and key[1] == "__GLOBAL__" and key[3] == component
                       for key in global_rows):
                global_rows[(geometry, "__GLOBAL__", "__SOURCE__", component, 0)] = zero.copy()
        for aa in selected_alphas:
            component = alpha_component[aa]
            for tt in selected_taus:
                if not any(key[:3] == (geometry, tt, aa) for key in landing):
                    landing[(geometry, tt, aa, component, 0)] = zero.copy()

    rows = global_rows | landing
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("geometry", "tau", "alpha", "source_component", "k_minus") + OUT_FIELDS)
        for (geometry, tt, aa, component, k), values in sorted(rows.items()):
            raw = dict(zip(VALUE_FIELDS, values))
            writer.writerow(
                (geometry, tt, aa, component, k,
                 raw["count"], raw["sum_q0"], raw["sum_q1"],
                 raw["sum_E0"], raw["sum_E1"], raw["sum_a16_0"], raw["sum_a16_1"],
                 raw["sum_q0_a16_0"], raw["sum_q1_a16_1"],
                 raw["sum_E0_a16_0"], raw["sum_E1_a16_1"])
            )
    cross_path = args.output.with_name("landing-cross-tab.json")
    cross_payload = {}
    for geometry, table in cross_tab.items():
        records = [
            {
                "occupied_landing_distinct": bool(key[0]),
                "vacant_cut_separators_distinct": bool(key[1]),
                "off_port_occupied_degree_branch": bool(key[2]),
                "extra_source_port_contact": bool(key[3]),
                "source_absent": bool(key[4]),
                "pair_fibres": value,
            }
            for key, value in sorted(table.items())
        ]
        strict = sum(
            item["pair_fibres"] for item in records
            if item["occupied_landing_distinct"]
            and item["vacant_cut_separators_distinct"]
            and not item["off_port_occupied_degree_branch"]
            and not item["extra_source_port_contact"]
        )
        forbidden_both_distinct = sum(
            item["pair_fibres"] for item in records
            if item["occupied_landing_distinct"] and item["vacant_cut_separators_distinct"]
        )
        landing_partition = {}
        for occupied_distinct in (False, True):
            for vacant_distinct in (False, True):
                pair_fibres = sum(
                    item["pair_fibres"] for item in records
                    if item["occupied_landing_distinct"] == occupied_distinct
                    and item["vacant_cut_separators_distinct"] == vacant_distinct
                )
                landing_partition[f"occupied_{'distinct' if occupied_distinct else 'merged'}__vacant_{'distinct' if vacant_distinct else 'merged'}"] = {
                    "pair_fibres": pair_fibres,
                    "backgrounds": pair_fibres // 23,
                    "pair_fibres_divisible_by_23": pair_fibres % 23 == 0,
                }
        cross_payload[geometry] = {
            "alternating_pair_fibres": sum(item["pair_fibres"] for item in records),
            "strict_ordinary_no_extra_pair_fibres": strict,
            "forbidden_both_distinct_pair_fibres": forbidden_both_distinct,
            "landing_partition": landing_partition,
            "records": records,
        }
    cross_path.write_text(json.dumps({
        "schema": "matching-one/p537-siteflip-landing-cross-tab/v1",
        "separator_semantics": "matching connectivity in the off-z cut graph",
        "strict_ordinary_definition": "alternating and occupied distinct and vacant distinct and no degree branch and no source-port extra contact",
        "result": "strict ordinary cell is empty; every alternating background has at least one colour merged",
        "geometries": cross_payload,
    }, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": "completed",
        "schema": "matching-one/p537-siteflip-schur-aggregate/v1",
        "landing_scope": "near_block_alternating_and_no_source_port_extra_contact",
        "alphas": len(selected_alphas),
        "source_components": len(selected_components),
        "taus": len(selected_taus),
        "rows": len(rows),
        "a_raw_denominator": 16,
        "fixed_z_orbit_multiplicity": 4,
        "output": str(args.output),
        "cross_tab": str(cross_path),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
