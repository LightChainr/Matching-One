#!/usr/bin/env python3
"""Aggregate the exact N25 radius-one collar into the first frozen 2x2 test.

Global source profiles retain all six C4 displacement components.  The landing
matrix is the direct preferred coarsening only: arm_mask=5, local source
contact zero, J_B=J_W=1, diagonal corner words summed, the two birth rows,
and axial2 absent/present columns.  Bell transitions are summed only after
their g16 sufficient statistics have been accumulated in the parent source
component.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


VALUE_FIELDS = (
    "count", "sum_q0", "sum_E0", "sum_a16_0", "sum_q0_a16_0", "sum_E0_a16_0",
    "sum_q1", "sum_E1", "sum_a16_1", "sum_q1_a16_1", "sum_E1_a16_1",
)
OUT_FIELDS = (
    "count", "sum_q0", "sum_q1", "sum_e0", "sum_e1", "sum_a0", "sum_a1",
    "sum_q0a0", "sum_q1a1", "sum_e0a0", "sum_e1a1",
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
BIRTH_ROWS = {(0, 1), (1, 2)}


def source_component(row: dict[str, str], geometry: str) -> str:
    displacement = (int(row["y_dx"]), int(row["y_dy"]))
    for name, members in ORBITS.items():
        if displacement in members:
            return name
    if displacement in FAR[geometry]:
        return "far"
    raise ValueError(f"unmapped {geometry} displacement {displacement}")


def add(target: dict, key: tuple, row: dict[str, str]) -> None:
    packet = target.setdefault(key, [0] * len(VALUE_FIELDS))
    for i, field in enumerate(VALUE_FIELDS):
        packet[i] += int(row[field])


def preferred(row: dict[str, str], geometry: str) -> bool:
    return (
        row["alternating_four_arm"] == "1"
        and row["arm_mask"] == "5"
        and row["outer_occupied_join"] == "1"
        and row["outer_vacant_join"] == "1"
        and row["local_source_contact_mask"] == "0"
        and (int(row["rank0"]), int(row["rank1"])) in BIRTH_ROWS
        and source_component(row, geometry) == "axial2"
    )


def tau(row: dict[str, str]) -> str:
    return f"collar_r1_birth:[{int(row['rank0'])},{int(row['rank1'])}]"


def alpha(row: dict[str, str]) -> str:
    return f"axial2:{'absent' if row['source_absent'] == '1' else 'present'}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--axis", type=Path, required=True)
    parser.add_argument("--tilted", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sources = {"axis": args.axis, "tilted": args.tilted}

    landing: dict[tuple[str, str, str, str, int], list[int]] = {}
    global_rows: dict[tuple[str, str, str, str, int], list[int]] = {}
    state_counts: dict[str, dict[tuple[int, int, int, int, int, int, int], int]] = {
        geometry: defaultdict(int) for geometry in sources
    }
    totals = {
        geometry: {
            "pair_fibres": 0,
            "alternating_pair_fibres": 0,
            "identity_failure_pair_fibres": 0,
            "forbidden_double_distinct_pair_fibres": 0,
            "preferred_pair_fibres": 0,
        }
        for geometry in sources
    }

    for geometry, path in sources.items():
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                component = source_component(row, geometry)
                displacement = (int(row["y_dx"]), int(row["y_dy"]))
                count = int(row["count"])
                totals[geometry]["pair_fibres"] += count
                if displacement == REPRESENTATIVE[geometry][component]:
                    add(global_rows, (geometry, "__GLOBAL__", "__SOURCE__", component,
                                      int(row["k_minus"])), row)

                if row["alternating_four_arm"] == "1":
                    totals[geometry]["alternating_pair_fibres"] += count
                    jb, jw = int(row["outer_occupied_join"]), int(row["outer_vacant_join"])
                    r0, r1 = int(row["rank0"]), int(row["rank1"])
                    if r1-r0 != jb+jw-1:
                        totals[geometry]["identity_failure_pair_fibres"] += count
                    if jb == 0 and jw == 0:
                        totals[geometry]["forbidden_double_distinct_pair_fibres"] += count
                    state_counts[geometry][(
                        int(row["collar_corner_mask"]), jb, jw,
                        int(row["local_source_contact_mask"]),
                        r0, r1, int(row["source_absent"]),
                    )] += count

                if preferred(row, geometry):
                    totals[geometry]["preferred_pair_fibres"] += count
                    add(landing, (geometry, tau(row), alpha(row), "axial2",
                                  int(row["k_minus"])), row)

    if any(packet["identity_failure_pair_fibres"] for packet in totals.values()):
        raise ValueError("producer emitted a finite-collar identity failure")
    if any(packet["forbidden_double_distinct_pair_fibres"] for packet in totals.values()):
        raise ValueError("producer emitted forbidden J_B=J_W=0")

    geometries = tuple(sources)
    taus = tuple(f"collar_r1_birth:[{r0},{r1}]" for r0, r1 in sorted(BIRTH_ROWS))
    alphas = ("axial2:absent", "axial2:present")
    zero = [0] * len(VALUE_FIELDS)
    for geometry in geometries:
        for tt in taus:
            for aa in alphas:
                if not any(key[:3] == (geometry, tt, aa) for key in landing):
                    landing[(geometry, tt, aa, "axial2", 0)] = zero.copy()

    rows = global_rows | landing
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("geometry", "tau", "alpha", "source_component", "k_minus") + OUT_FIELDS)
        for (geometry, tt, aa, component, k), values in sorted(rows.items()):
            raw = dict(zip(VALUE_FIELDS, values))
            writer.writerow((
                geometry, tt, aa, component, k,
                raw["count"], raw["sum_q0"], raw["sum_q1"], raw["sum_E0"], raw["sum_E1"],
                raw["sum_a16_0"], raw["sum_a16_1"], raw["sum_q0_a16_0"], raw["sum_q1_a16_1"],
                raw["sum_E0_a16_0"], raw["sum_E1_a16_1"],
            ))

    state_path = args.output.with_name("collar-state-counts.json")
    state_path.write_text(json.dumps({
        "schema": "matching-one/p537-finite-collar-state-counts/v1",
        "collar": "B_inf(z,1)",
        "corner_bit_order": ["NE", "SE", "SW", "NW"],
        "arm_order": ["B_N", "W_E", "B_S", "W_W"],
        "outer_identity": "rank1-rank0=J_B+J_W-1",
        "preferred_coarsening": {
            "arm_mask": 5,
            "outer_occupied_join": 1,
            "outer_vacant_join": 1,
            "local_source_contact_mask": 0,
            "corner_masks": "summed",
            "rows": [[0, 1], [1, 2]],
            "columns": ["axial2:absent", "axial2:present"],
        },
        "geometries": {
            geometry: {
                **totals[geometry],
                "records": [
                    {
                        "corner_mask": key[0], "J_B": key[1], "J_W": key[2],
                        "local_source_contact_mask": key[3], "rank0": key[4], "rank1": key[5],
                        "source_absent": bool(key[6]), "pair_fibres": value,
                    }
                    for key, value in sorted(state_counts[geometry].items())
                ],
            }
            for geometry in geometries
        },
    }, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": "completed",
        "schema": "matching-one/p537-finite-collar-schur-aggregate/v1",
        "global_source_components": sorted(REPRESENTATIVE["axis"]),
        "landing_source_component": "axial2",
        "taus": list(taus),
        "alphas": list(alphas),
        "rows": len(rows),
        "a_raw_denominator": 16,
        "fixed_z_orbit_multiplicity": 4,
        "output": str(args.output),
        "state_counts": str(state_path),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
