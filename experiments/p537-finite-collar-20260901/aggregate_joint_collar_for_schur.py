#!/usr/bin/env python3
"""Aggregate the frozen P537 collar by simultaneous x+y+z incidence sector.

The joint producer already applies the frozen radius-one filter.  This program
does not add a descriptor or choose a sector after scoring: it partitions the
same coarse four cells by the first-occurrence canonical partition of the
ordered physical terminals ``x[N,E,S,W], y[N,E,S,W], z[N,E,S,W]``.  Diagonal
corner words remain summed.  The integer sufficient statistics are required
to sum exactly back to the predeclared coarse aggregate before any Schur score
is evaluated.

The sector is frozen before either source column is selected: x, y and z are
vacant in the common base, then z is toggled to record the transported
partition.  The same pair of identities labels both actual y=0 and y=1
columns.  This makes a same-sector 2x2 test meaningful instead of comparing
the source-absent sentinel to an unrelated ordinary-source partition.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
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
JOINT_FIELDS = {
    "k_minus", "rank0", "rank1", "collar_corner_mask", "source_absent",
    "base_arm_mask", "source_port_occupied_mask", "y_dx", "y_dy",
    "joint0", "joint1", "terminal_incidence0",
    "terminal_incidence1", "z_source_roles0", "z_source_roles1", *VALUE_FIELDS,
}
ROWS = ((0, 1), (1, 2))
ALPHAS = ("axial2:absent", "axial2:present")
GEOMETRIES = ("axis", "tilted")
GLOBAL = "__GLOBAL__"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def packet(row: dict[str, str], fields: tuple[str, ...] = VALUE_FIELDS) -> list[int]:
    return [int(row[field]) for field in fields]


def add(target: dict, key: tuple, values: list[int]) -> None:
    current = target.setdefault(key, [0] * len(values))
    for i, value in enumerate(values):
        current[i] += value


def tau(row: dict[str, str]) -> str:
    return f"collar_r1_birth:[{int(row['rank0'])},{int(row['rank1'])}]"


def alpha(row: dict[str, str]) -> str:
    return f"axial2:{'absent' if row['source_absent'] == '1' else 'present'}"


def signature(row: dict[str, str]) -> str:
    return (
        f"d={int(row['y_dx'])},{int(row['y_dy'])};"
        f"arm={int(row['base_arm_mask'])};"
        f"ports={int(row['source_port_occupied_mask'])};"
        f"j0={int(row['joint0'])};j1={int(row['joint1'])}"
    )


def sector_id(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()[:20]


def read_coarse(path: Path) -> tuple[list[dict[str, str]], dict[tuple, list[int]]]:
    globals_: list[dict[str, str]] = []
    landing: dict[tuple, list[int]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (
                row["geometry"], row["tau"], row["alpha"],
                row["source_component"], int(row["k_minus"]),
            )
            values = packet(row, OUT_FIELDS)
            if row["tau"] == GLOBAL:
                globals_.append(row)
            else:
                if key in landing:
                    raise ValueError(f"duplicate coarse key {key}")
                landing[key] = values
    return globals_, landing


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--axis-joint", type=Path, required=True)
    parser.add_argument("--tilted-joint", type=Path, required=True)
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    args = parser.parse_args()

    global_rows, expected_coarse = read_coarse(args.coarse)
    joint: dict[tuple[str, str, str, str, str, int], list[int]] = {}
    coarsened: dict[tuple[str, str, str, str, int], list[int]] = {}
    support: dict[str, dict[tuple[str, str, str], int]] = defaultdict(lambda: defaultdict(int))
    metadata: dict[str, dict[str, object]] = {}
    sources = {"axis": args.axis_joint, "tilted": args.tilted_joint}

    for geometry, path in sources.items():
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            missing = JOINT_FIELDS - set(reader.fieldnames or ())
            if missing:
                raise ValueError(f"{path}: missing fields {sorted(missing)}")
            for row in reader:
                text = signature(row)
                sector = sector_id(text)
                meta = {
                    "signature": text,
                    "y_displacement": [int(row["y_dx"]), int(row["y_dy"])],
                    "base_arm_mask": int(row["base_arm_mask"]),
                    "source_port_occupied_mask": int(row["source_port_occupied_mask"]),
                    "joint0": int(row["joint0"]),
                    "joint1": int(row["joint1"]),
                    "terminal_incidence0": int(row["terminal_incidence0"]),
                    "terminal_incidence1": int(row["terminal_incidence1"]),
                    "z_source_roles0": int(row["z_source_roles0"]),
                    "z_source_roles1": int(row["z_source_roles1"]),
                }
                previous = metadata.setdefault(sector, meta)
                if previous != meta:
                    raise ValueError(f"sector hash collision/inconsistent metadata {sector}")
                tt, aa, k = tau(row), alpha(row), int(row["k_minus"])
                if tt not in {"collar_r1_birth:[0,1]", "collar_r1_birth:[1,2]"}:
                    raise ValueError(f"joint producer emitted non-birth row {tt}")
                values = packet(row)
                add(joint, (geometry, sector, tt, aa, "axial2", k), values)
                add(coarsened, (geometry, tt, aa, "axial2", k), values)
                support[sector][(geometry, tt, aa)] += values[0]

    # The new partition must be an exact integer refinement of every old landing row.
    zero_out = [0] * len(OUT_FIELDS)
    normalized_expected = {key: value for key, value in expected_coarse.items() if value != zero_out}
    normalized_actual: dict[tuple, list[int]] = {}
    for key, values in coarsened.items():
        raw = dict(zip(VALUE_FIELDS, values))
        normalized_actual[key] = [
            raw["count"], raw["sum_q0"], raw["sum_q1"], raw["sum_E0"], raw["sum_E1"],
            raw["sum_a16_0"], raw["sum_a16_1"], raw["sum_q0_a16_0"], raw["sum_q1_a16_1"],
            raw["sum_E0_a16_0"], raw["sum_E1_a16_1"],
        ]
    normalized_actual = {key: value for key, value in normalized_actual.items() if value != zero_out}
    if normalized_actual != normalized_expected:
        missing = sorted(set(normalized_expected) - set(normalized_actual))[:3]
        extra = sorted(set(normalized_actual) - set(normalized_expected))[:3]
        mismatch = [key for key in normalized_expected.keys() & normalized_actual.keys()
                    if normalized_expected[key] != normalized_actual[key]][:3]
        raise ValueError(f"joint/coarse mismatch missing={missing} extra={extra} values={mismatch}")

    taus = ("collar_r1_birth:[0,1]", "collar_r1_birth:[1,2]")
    for sector in metadata:
        for geometry in GEOMETRIES:
            for tt in taus:
                for aa in ALPHAS:
                    # Each input exhausts the common counterfactual base.  A
                    # missing actual cell is therefore a proved zero in this
                    # geometry/sector, not an imputed coordinate.
                    if support[sector].get((geometry, tt, aa), 0) == 0:
                        joint[(geometry, sector, tt, aa, "axial2", 0)] = [0] * len(VALUE_FIELDS)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("joint_sector", "geometry", "tau", "alpha", "source_component", "k_minus") + OUT_FIELDS)
        for row in sorted(global_rows, key=lambda x: (
            x["geometry"], x["tau"], x["alpha"], x["source_component"], int(x["k_minus"])
        )):
            writer.writerow((GLOBAL,) + tuple(row[field] for field in (
                "geometry", "tau", "alpha", "source_component", "k_minus", *OUT_FIELDS
            )))
        for (geometry, sector, tt, aa, component, k), values in sorted(joint.items()):
            raw = dict(zip(VALUE_FIELDS, values))
            writer.writerow((
                sector, geometry, tt, aa, component, k,
                raw["count"], raw["sum_q0"], raw["sum_q1"], raw["sum_E0"], raw["sum_E1"],
                raw["sum_a16_0"], raw["sum_a16_1"], raw["sum_q0_a16_0"], raw["sum_q1_a16_1"],
                raw["sum_E0_a16_0"], raw["sum_E1_a16_1"],
            ))

    cell_names = [(geometry, tt, aa) for geometry in GEOMETRIES for tt in taus for aa in ALPHAS]
    records = []
    for sector in sorted(metadata):
        counts = {"|".join(key): support[sector].get(key, 0) for key in cell_names}
        records.append({
            "sector": sector,
            **metadata[sector],
            "cell_pair_fibres": counts,
            "full_axis_tilted_rectangle": all(counts.values()),
            "pooled_four_cell_rectangle": all(
                sum(support[sector].get((geometry, tt, aa), 0) for geometry in GEOMETRIES)
                for tt in taus for aa in ALPHAS
            ),
        })
    payload = {
        "schema": "matching-one/p537-finite-collar-joint-index/v1",
        "status": "exact_integer_refinement_reproduces_coarse",
        "joint_key": (
            "first-occurrence canon_global(x4+y4+z4), ordered N/E/S/W; "
            "common x=y=z=0 base and transported z=1 partition"
        ),
        "column_transport": (
            "the frozen base identity is shared by actual y=0 ordinary-source and "
            "actual y=1 source-absent cells"
        ),
        "coarse_reproduction": True,
        "coarse_sha256": sha256(args.coarse),
        "joint_input_sha256": {geometry: sha256(path) for geometry, path in sources.items()},
        "joint_aggregate_sha256": sha256(args.output),
        "sectors": records,
        "sector_count": len(records),
        "full_axis_tilted_rectangle_count": sum(row["full_axis_tilted_rectangle"] for row in records),
        "pooled_four_cell_rectangle_count": sum(row["pooled_four_cell_rectangle"] for row in records),
        "boundary": (
            "Exact finite N25 two-geometry refinement of the frozen r1 collar only. "
            "Zero-filled landing cells are exhaustive-enumeration zeros. This does not "
            "supply thermal jets, quantitative landing margins, cross-size persistence, "
            "or a macroscopic four-arm theorem."
        ),
    }
    args.index.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": payload["status"],
        "sectors": payload["sector_count"],
        "full_rectangles": payload["full_axis_tilted_rectangle_count"],
        "pooled_rectangles": payload["pooled_four_cell_rectangle_count"],
        "output": str(args.output),
        "index": str(args.index),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
