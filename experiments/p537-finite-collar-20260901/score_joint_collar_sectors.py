#!/usr/bin/env python3
"""Score the frozen P537 12-port joint-incidence sectors.

This is a bounded refinement of the already frozen finite-collar Schur test.
The global root, thermal projection, source means and beta coefficients remain
those of the parent axial2 component.  Only the four landing cells are split
by the counterfactual x+y+z joint-incidence sector.  Consequently the integer
sector packets must coarsen exactly to the committed four-cell aggregate, and
the sector matrices add before taking the determinant.

The determinant of that sum is decomposed into same-sector determinants plus
cross-sector terms.  A gate-quality same-sector witness is reported only when
all four pooled row/column cells have positive exact support; the stronger
``full_axis_tilted_rectangle`` flag requires all eight geometry-specific cells.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import sys
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
BASE_SCORER = HERE.parent / "p537-siteflip-landing-20260901" / "score_siteflip_schur.py"
BASE_SCORER_DISPLAY = BASE_SCORER.relative_to(HERE.parent.parent)
SPEC = importlib.util.spec_from_file_location("p537_base_schur", BASE_SCORER)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load base scorer {BASE_SCORER}")
base = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = base
SPEC.loader.exec_module(base)


GLOBAL = "__GLOBAL__"
ROWS = ("collar_r1_birth:[0,1]", "collar_r1_birth:[1,2]")
COLS = ("axial2:absent", "axial2:present")
GEOMETRIES = ("axis", "tilted")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def interval_from_record(record: dict[str, object]) -> base.Interval:
    return base.Interval(F(str(record["lower"])), F(str(record["upper"])))


def interval_midpoint(value: base.Interval) -> F:
    return (value.lo + value.hi) / 2


def determinant(matrix: list[list[base.Interval]]) -> base.Interval:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def read_rows(path: Path, n: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    required = set(base.REQUIRED_FIELDS) | {"joint_sector"}
    seen: set[tuple[str, str, str, str, str, int]] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing fields {sorted(missing)}")
        for line_number, raw in enumerate(reader, start=2):
            sector = raw["joint_sector"].strip()
            geometry = raw["geometry"].strip()
            tau = raw["tau"].strip()
            alpha = raw["alpha"].strip()
            component = raw["source_component"].strip()
            k = int(raw["k_minus"])
            if geometry not in GEOMETRIES:
                raise ValueError(f"{path}:{line_number}: unknown geometry {geometry!r}")
            if not sector or not tau or not alpha or not component:
                raise ValueError(f"{path}:{line_number}: empty sector/tau/alpha/component")
            if not 0 <= k < n:
                raise ValueError(f"{path}:{line_number}: k_minus outside [0,{n - 1}]")
            key = (sector, geometry, tau, alpha, component, k)
            if key in seen:
                raise ValueError(f"{path}:{line_number}: duplicate key {key}")
            seen.add(key)
            row: dict[str, object] = {
                "sector": sector,
                "geometry": geometry,
                "tau": tau,
                "alpha": alpha,
                "component": component,
                "k": k,
            }
            for field in base.REQUIRED_FIELDS - {
                "geometry", "tau", "alpha", "source_component", "k_minus"
            }:
                row[field] = F(raw[field])
            if row["count"] < 0:
                raise ValueError(f"{path}:{line_number}: negative count")
            rows.append(row)
    return rows


def score(args: argparse.Namespace) -> dict[str, object]:
    n = args.n
    delta = base.parse_fraction(args.delta)
    p = base.read_root(args.baseline_root)
    rows = read_rows(args.aggregates, n)
    global_rows = [row for row in rows if row["tau"] == GLOBAL]
    landing_rows = [row for row in rows if row["tau"] != GLOBAL]
    if any(row["sector"] != GLOBAL for row in global_rows):
        raise ValueError("global source rows must use the __GLOBAL__ sector")
    if any(row["sector"] == GLOBAL for row in landing_rows):
        raise ValueError("landing rows must use a concrete joint sector")

    coefficients = {
        "axis": base.read_baseline(args.baseline_axis, n),
        "tilted": base.read_baseline(args.baseline_tilted, n),
    }
    baseline = {
        geometry: base.baseline_packet(coefficients[geometry], n, p)
        for geometry in GEOMETRIES
    }
    mt = (baseline["axis"]["q_t"] + baseline["tilted"]["q_t"]) / 2
    yt = (baseline["axis"]["e_t"] - baseline["tilted"]["e_t"]) / delta
    r = yt / mt
    c = {"axis": F(1, 1) / delta, "tilted": -F(1, 1) / delta}
    mu_h = {
        geometry: 2 * c[geometry] * baseline[geometry]["e"]
        - r * baseline[geometry]["q"]
        for geometry in GEOMETRIES
    }
    source_packets, beta = base.source_global_packets(
        global_rows, baseline, n, p, args.a_raw_denominator
    )
    # The exact root interval has very large rational endpoints.  Cache the 25
    # repeated Bernoulli weights rather than exponentiating them once for every
    # sector packet; this changes no arithmetic or grouping contract.
    offsite_weights = {k: base.offsite_weight(k, n, p) for k in range(n)}

    geometry_cells: dict[tuple[str, str, str, str], base.Interval] = defaultdict(
        lambda: base.Interval.of(0)
    )
    sectors: set[str] = set()
    for row in landing_rows:
        sector = str(row["sector"])
        geometry = str(row["geometry"])
        tau = str(row["tau"])
        alpha = str(row["alpha"])
        component = str(row["component"])
        if tau not in ROWS or alpha not in COLS or component != "axial2":
            raise ValueError(f"unexpected landing label {(tau, alpha, component)}")
        sectors.add(sector)
        k = int(row["k"])
        weight = offsite_weights[k]
        mu_a = source_packets[(geometry, component)]["mu_a"]
        beta_component = beta[component]
        s_minus = k - (n - 1) * p
        for i in (0, 1):
            wi = (1 - p) if i == 0 else p
            ui = i - p
            si = s_minus + ui
            bi = ui * si - p * (1 - p)
            suffix = str(i)
            count = row["count"]
            sum_q = base.allocated(row, f"sum_q{suffix}", n, args.a_raw_denominator)
            sum_e = base.allocated(row, f"sum_e{suffix}", n, args.a_raw_denominator)
            sum_a = base.allocated(row, f"sum_a{suffix}", n, args.a_raw_denominator)
            sum_qa = base.allocated(row, f"sum_q{suffix}a{suffix}", n, args.a_raw_denominator)
            sum_ea = base.allocated(row, f"sum_e{suffix}a{suffix}", n, args.a_raw_denominator)
            sum_h = 2 * c[geometry] * sum_e - r * sum_q - mu_h[geometry] * count
            sum_ha = (
                2 * c[geometry] * (sum_ea - mu_a * sum_e)
                - r * (sum_qa - mu_a * sum_q)
                - mu_h[geometry] * (sum_a - mu_a * count)
            )
            geometry_cells[(sector, geometry, tau, alpha)] += weight * base.eq10_state_term(
                wi, ui, bi, sum_h, sum_ha, beta_component,
                args.z_orbit_multiplicity,
            )

    index = json.loads(args.index.read_text())
    if not index.get("coarse_reproduction"):
        raise ValueError("joint index does not certify exact coarse reproduction")
    metadata = {row["sector"]: row for row in index["sectors"]}
    if sectors != set(metadata):
        raise ValueError(
            f"aggregate/index sector mismatch aggregate_only={sorted(sectors-set(metadata))[:3]} "
            f"index_only={sorted(set(metadata)-sectors)[:3]}"
        )

    matrices: dict[str, list[list[base.Interval]]] = {}
    determinants: dict[str, base.Interval] = {}
    matrix_sum = [[base.Interval.of(0) for _ in COLS] for _ in ROWS]
    same_sector_det_sum = base.Interval.of(0)
    records: list[dict[str, object]] = []
    for sector in sorted(sectors):
        matrix = [
            [
                (
                    geometry_cells[(sector, "axis", tau, alpha)]
                    + geometry_cells[(sector, "tilted", tau, alpha)]
                ) / 2
                for alpha in COLS
            ]
            for tau in ROWS
        ]
        det = determinant(matrix)
        matrices[sector] = matrix
        determinants[sector] = det
        same_sector_det_sum += det
        for i in range(2):
            for j in range(2):
                matrix_sum[i][j] += matrix[i][j]
        det_record = base.interval_record(det)
        records.append({
            "sector": sector,
            "full_axis_tilted_rectangle": bool(metadata[sector]["full_axis_tilted_rectangle"]),
            "pooled_four_cell_rectangle": bool(metadata[sector]["pooled_four_cell_rectangle"]),
            "determinant": det_record,
            "P4_Schur": [[base.interval_record(cell) for cell in row] for row in matrix],
        })

    joint_sum_det = determinant(matrix_sum)
    cross_sector = joint_sum_det - same_sector_det_sum

    coarse_args = argparse.Namespace(
        aggregates=args.coarse,
        baseline_axis=args.baseline_axis,
        baseline_tilted=args.baseline_tilted,
        baseline_root=args.baseline_root,
        n=n,
        delta=args.delta,
        a_raw_denominator=args.a_raw_denominator,
        z_orbit_multiplicity=args.z_orbit_multiplicity,
        first_nonzero_only=False,
    )
    coarse_score = base.score(coarse_args)
    coarse_matrix = [
        [interval_from_record(cell) for cell in row]
        for row in coarse_score["matrix"]["P4_Schur"]
    ]
    coarse_det = determinant(coarse_matrix)
    cell_comparison = []
    for i, tau in enumerate(ROWS):
        for j, alpha in enumerate(COLS):
            joint_cell, coarse_cell = matrix_sum[i][j], coarse_matrix[i][j]
            cell_comparison.append({
                "tau": tau,
                "alpha": alpha,
                "intervals_overlap": not (
                    joint_cell.hi < coarse_cell.lo or coarse_cell.hi < joint_cell.lo
                ),
                "joint_midpoint": float(interval_midpoint(joint_cell)),
                "coarse_midpoint": float(interval_midpoint(coarse_cell)),
                "midpoint_residual": float(
                    interval_midpoint(joint_cell) - interval_midpoint(coarse_cell)
                ),
            })

    complete_pooled = [
        row for row in records if row["pooled_four_cell_rectangle"]
    ]
    complete_full = [
        row for row in records if row["full_axis_tilted_rectangle"]
    ]
    first_pooled_nonzero = next(
        (row for row in complete_pooled if row["determinant"]["excludes_zero"]), None
    )
    first_full_nonzero = next(
        (row for row in complete_full if row["determinant"]["excludes_zero"]), None
    )
    top = sorted(
        records,
        key=lambda row: abs(float(row["determinant"]["midpoint"])),
        reverse=True,
    )[:12]

    args.sector_output.parent.mkdir(parents=True, exist_ok=True)
    with args.sector_output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow((
            "sector", "full_axis_tilted_rectangle", "pooled_four_cell_rectangle",
            "det_lower", "det_upper", "det_midpoint", "det_excludes_zero",
            "p4_01_absent", "p4_01_present", "p4_12_absent", "p4_12_present",
        ))
        for record in records:
            det_record = record["determinant"]
            p4 = record["P4_Schur"]
            writer.writerow((
                record["sector"], int(record["full_axis_tilted_rectangle"]),
                int(record["pooled_four_cell_rectangle"]), det_record["lower"],
                det_record["upper"], det_record["midpoint"],
                int(det_record["excludes_zero"]), p4[0][0]["midpoint"],
                p4[0][1]["midpoint"], p4[1][0]["midpoint"],
                p4[1][1]["midpoint"],
            ))

    coarse_mid = interval_midpoint(coarse_det)
    cross_mid = interval_midpoint(cross_sector)
    same_mid = interval_midpoint(same_sector_det_sum)
    status = (
        "exact_nonzero_full_joint_sector"
        if first_full_nonzero is not None
        else "exact_nonzero_pooled_joint_sector"
        if first_pooled_nonzero is not None
        else "no_exact_nonzero_complete_joint_sector"
    )
    payload = {
        "schema": "matching-one/p537-finite-collar-joint-score/v1",
        "status": status,
        "N": n,
        "root_p": base.interval_record(p),
        "row_order": list(ROWS),
        "column_order": list(COLS),
        "sector_count": len(records),
        "full_axis_tilted_rectangle_count": len(complete_full),
        "pooled_four_cell_rectangle_count": len(complete_pooled),
        "exact_nonzero_sector_count": sum(
            bool(row["determinant"]["excludes_zero"]) for row in records
        ),
        "exact_nonzero_pooled_rectangle_count": sum(
            bool(row["determinant"]["excludes_zero"]) for row in complete_pooled
        ),
        "exact_nonzero_full_rectangle_count": sum(
            bool(row["determinant"]["excludes_zero"]) for row in complete_full
        ),
        "first_exact_nonzero_pooled_rectangle": first_pooled_nonzero,
        "first_exact_nonzero_full_rectangle": first_full_nonzero,
        "top_absolute_sector_determinants": top,
        "coarse_reproduction": {
            "integer_packets_exact": True,
            "matrix_cell_comparison": cell_comparison,
            "coarse_determinant": base.interval_record(coarse_det),
            "determinant_of_joint_matrix_sum": base.interval_record(joint_sum_det),
        },
        "sector_mixing_decomposition": {
            "identity": "det(sum_s L_s)=sum_s det(L_s)+cross_sector_terms",
            "same_sector_determinant_sum": base.interval_record(same_sector_det_sum),
            "cross_sector_terms": base.interval_record(cross_sector),
            "same_sector_midpoint_fraction_of_coarse": float(same_mid / coarse_mid),
            "cross_sector_midpoint_fraction_of_coarse": float(cross_mid / coarse_mid),
        },
        "inputs": {
            "joint_aggregates": {"path": str(args.aggregates), "sha256": sha256(args.aggregates)},
            "joint_index": {"path": str(args.index), "sha256": sha256(args.index)},
            "coarse_aggregates": {"path": str(args.coarse), "sha256": sha256(args.coarse)},
            "baseline_axis": {"path": str(args.baseline_axis), "sha256": sha256(args.baseline_axis)},
            "baseline_tilted": {"path": str(args.baseline_tilted), "sha256": sha256(args.baseline_tilted)},
            "baseline_root": {"path": str(args.baseline_root), "sha256": sha256(args.baseline_root)},
            "base_scorer": {"path": str(BASE_SCORER_DISPLAY), "sha256": sha256(BASE_SCORER)},
            "sector_scores": {"path": str(args.sector_output), "sha256": sha256(args.sector_output)},
        },
        "boundary": (
            "Exact finite N25 value-level Schur refinement for two geometries. "
            "It does not provide thermal-jet independence, quantitative landing "
            "margins, cross-size persistence, or an asymptotic arm-event theorem."
        ),
    }
    return payload


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--aggregates", type=Path, required=True)
    result.add_argument("--index", type=Path, required=True)
    result.add_argument("--coarse", type=Path, required=True)
    result.add_argument("--baseline-axis", type=Path, required=True)
    result.add_argument("--baseline-tilted", type=Path, required=True)
    result.add_argument("--baseline-root", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--sector-output", type=Path, required=True)
    result.add_argument("--n", type=int, default=25)
    result.add_argument("--delta", default="1152/625")
    result.add_argument("--a-raw-denominator", type=int, default=16)
    result.add_argument("--z-orbit-multiplicity", type=int, default=4)
    return result


def main() -> None:
    args = parser().parse_args()
    if args.output.exists() or args.sector_output.exists():
        raise FileExistsError("refusing to overwrite output")
    payload = score(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "status": payload["status"],
        "sectors": payload["sector_count"],
        "pooled_rectangles": payload["pooled_four_cell_rectangle_count"],
        "full_rectangles": payload["full_axis_tilted_rectangle_count"],
        "output": str(args.output),
    }, sort_keys=True))


if __name__ == "__main__":
    main()
