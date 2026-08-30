#!/usr/bin/env python3
"""Medium-Gaussian intervention-sensitive pilot for Issues 250 and 333.

This transports the exact adaptive D/J rule to one N325 and one N425
Gaussian torus.  Undefined partial rectangles contribute zero to the primary
unconditional response and are reported separately.  Every defined row is
paired with its exact complement plus NN/matching-hand involution.
"""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from math import cos, floor, pi, sin, sqrt
from pathlib import Path
import platform
import statistics
import sys
import time
from typing import Sequence

try:
    from integer_period_torus import (
        IntegerHomologyUnionFind,
        gaussian_integer_torus,
    )
except ModuleNotFoundError:
    from scripts.integer_period_torus import (
        IntegerHomologyUnionFind,
        gaussian_integer_torus,
    )


MASK64 = (1 << 64) - 1
P_REF = 0.592746050790
GEOMETRIES = {
    "N325": {
        "a": 17,
        "b": 6,
        "N": 325,
        "period_matrix": [[17, -6], [6, 17]],
        "seed": 25033332520260830,
        "replica_offset": 10_250_000_000,
    },
    "N425": {
        "a": 16,
        "b": 13,
        "N": 425,
        "period_matrix": [[16, -13], [13, 16]],
        "seed": 25033342520260830,
        "replica_offset": 10_260_000_000,
    },
}


def splitmix64(value: int) -> int:
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return (value ^ (value >> 31)) & MASK64


def counter_uniform(seed: int, replica: int, site: int, n: int) -> float:
    key = (
        seed
        ^ splitmix64(n)
        ^ splitmix64(replica + 0xD1B54A32D192ED03)
        ^ splitmix64(site + 0x94D049BB133111EB)
    )
    return (splitmix64(key) >> 11) * (2.0**-53)


def transform_offset(orientation: int, point: tuple[int, int]) -> tuple[int, int]:
    """Four rotations followed by their reflected partners."""

    if not 0 <= orientation < 8:
        raise ValueError("orientation must be in 0..7")
    x, y = point
    if orientation >= 4:
        y = -y
    for _ in range(orientation % 4):
        x, y = -y, x
    return x, y


class TypedIndex:
    def __init__(self, oracle: "MediumOracle", active: Sequence[bool], hand: bool):
        self.oracle = oracle
        self.active = tuple(bool(value) for value in active)
        self.hand = hand
        self.black_graph = "matching" if hand else "NN"
        self.white_graph = "NN" if hand else "matching"
        self.black = self._build(self.active, self.black_graph)
        white = tuple(not value for value in self.active)
        self.white = self._build(white, self.white_graph)
        self._component_cache: dict[tuple[bool, int], dict] = {}

    def _build(self, enabled: Sequence[bool], graph: str):
        union = IntegerHomologyUnionFind(
            self.oracle.geometry.n, self.oracle.geometry.periods
        )
        edges = (
            self.oracle.geometry.matching_edges
            if graph == "matching"
            else self.oracle.geometry.primal_edges
        )
        for edge in edges:
            if enabled[edge.i] and enabled[edge.j]:
                union.add_edge(edge.i, edge.j, edge.dx, edge.dy)
        return union

    def component(self, occupied: bool, anchor: int) -> dict | None:
        key = occupied, anchor
        if key in self._component_cache:
            return self._component_cache[key]
        if self.active[anchor] != occupied:
            return None
        union = self.black if occupied else self.white
        root, _, _ = union.find(anchor)
        vertices = tuple(
            vertex
            for vertex, value in enumerate(self.active)
            if value == occupied and union.find(vertex)[0] == root
        )
        exact = union.component(anchor)
        result = {
            "id": root,
            "vertices": vertices,
            "vertex_set": frozenset(vertices),
            "size": len(vertices),
            "rank": exact.rank,
            "basis": tuple(tuple(vector) for vector in exact.basis),
            "colour": "occupied" if occupied else "vacant",
            "graph": self.black_graph if occupied else self.white_graph,
        }
        self._component_cache[key] = result
        return result

    def leg(self, root: int) -> int:
        occupied = self.active[root]
        component = self.component(occupied, root)
        if component is None:
            raise AssertionError("root is absent from its typed colour component")
        value = int(component["rank"] == 1)
        return value if occupied else -value


class MediumOracle:
    def __init__(self, geometry_id: str):
        spec = GEOMETRIES[geometry_id]
        self.geometry_id = geometry_id
        self.a = int(spec["a"])
        self.b = int(spec["b"])
        self.geometry = gaussian_integer_torus(self.a, self.b)
        if self.geometry.n != spec["N"]:
            raise AssertionError("frozen Gaussian order changed")
        self.neighbours = {
            "NN": self._neighbours(self.geometry.primal_edges),
            "matching": self._neighbours(self.geometry.matching_edges),
        }
        self.distance_from_origin = self._distance_table()

    def _neighbours(self, edges) -> tuple[tuple[int, ...], ...]:
        rows = [set() for _ in range(self.geometry.n)]
        for edge in edges:
            rows[edge.i].add(edge.j)
            rows[edge.j].add(edge.i)
        return tuple(tuple(sorted(row)) for row in rows)

    def _distance_table(self) -> tuple[int, ...]:
        """Exact nearest-image squared distance for the orthogonal Gaussian basis."""

        n = self.geometry.n
        output = []
        for x, y in self.geometry.coordinates:
            coefficient_0 = (self.a * x + self.b * y) / n
            coefficient_1 = (-self.b * x + self.a * y) / n
            base_0 = floor(coefficient_0)
            base_1 = floor(coefficient_1)
            candidates = []
            for u in range(base_0 - 2, base_0 + 3):
                for v in range(base_1 - 2, base_1 + 3):
                    dx = x - (self.a * u - self.b * v)
                    dy = y - (self.b * u + self.a * v)
                    candidates.append(dx * dx + dy * dy)
            output.append(min(candidates))
        return tuple(output)

    def distance_squared(self, first: int, second: int) -> int:
        x0, y0 = self.geometry.coordinates[first]
        x1, y1 = self.geometry.coordinates[second]
        difference = self.geometry.vertex((x0 - x1, y0 - y1))
        return self.distance_from_origin[difference]

    def marks(self, seed: int, replica: int) -> tuple[int, int, int, int]:
        orientation = replica % 8
        translation = splitmix64(seed ^ splitmix64(replica)) % self.geometry.n
        tx, ty = self.geometry.coordinates[translation]
        offsets = ((0, 0), (1, 1), (0, 1))
        marked = []
        for offset in offsets:
            dx, dy = transform_offset(orientation, offset)
            marked.append(self.geometry.vertex((tx + dx, ty + dy)))
        if len(set(marked)) != 3:
            raise AssertionError("medium marked triple collided")
        return marked[0], marked[1], marked[2], orientation

    def support(
        self,
        index: TypedIndex,
        anchor_D: int,
        anchor_J: int,
        landing: int,
        operation: str,
    ) -> dict | None:
        occupied = operation == "D"
        source_anchor = anchor_D if occupied else anchor_J
        target_anchor = anchor_J if occupied else anchor_D
        source = index.component(occupied, source_anchor)
        target = index.component(not occupied, target_anchor)
        if source is None or target is None or source["rank"] < 1:
            return None
        target_graph = target["graph"]
        target_vertices = target["vertex_set"]
        transferable = []
        for vertex in source["vertices"]:
            if vertex in (anchor_D, anchor_J):
                continue
            if any(
                neighbour in target_vertices
                for neighbour in self.neighbours[target_graph][vertex]
            ):
                transferable.append(vertex)
        if not transferable:
            return None
        distances = [self.distance_squared(vertex, landing) for vertex in transferable]
        minimum = min(distances)
        minimizers = tuple(
            vertex
            for vertex, distance in zip(transferable, distances)
            if distance == minimum
        )
        if len(minimizers) != 1:
            return None
        return {
            "site": minimizers[0],
            "source_rank": source["rank"],
            "source_size": source["size"],
            "source_basis": source["basis"],
            "source_graph": source["graph"],
            "target_rank": target["rank"],
            "target_size": target["size"],
            "target_basis": target["basis"],
            "target_graph": target["graph"],
            "distance_squared": minimum,
            "minimizer_count": 1,
        }

    @staticmethod
    def flipped(active: Sequence[bool], site: int, occupied: bool) -> list[bool]:
        changed = list(active)
        changed[site] = occupied
        return changed

    def rectangle(
        self,
        active: Sequence[bool],
        anchor_D: int,
        anchor_J: int,
        landing: int,
        *,
        hand: bool,
    ) -> dict | None:
        if not active[anchor_D] or active[anchor_J]:
            return None
        base = TypedIndex(self, active, hand)
        support_D = self.support(base, anchor_D, anchor_J, landing, "D")
        support_J = self.support(base, anchor_D, anchor_J, landing, "J")
        if support_D is None or support_J is None:
            return None
        if support_D["site"] == support_J["site"]:
            raise AssertionError("opposite-colour base supports collided")
        field_D = self.flipped(active, support_D["site"], False)
        field_J = self.flipped(active, support_J["site"], True)
        index_D = TypedIndex(self, field_D, hand)
        index_J = TypedIndex(self, field_J, hand)
        support_J_after_D = self.support(index_D, anchor_D, anchor_J, landing, "J")
        support_D_after_J = self.support(index_J, anchor_D, anchor_J, landing, "D")
        if support_J_after_D is None or support_D_after_J is None:
            return None
        field_DJ = self.flipped(field_D, support_J_after_D["site"], True)
        field_JD = self.flipped(field_J, support_D_after_J["site"], False)
        index_DJ = TypedIndex(self, field_DJ, hand)
        index_JD = TypedIndex(self, field_JD, hand)
        responses = {
            "L_D": index_D.leg(landing),
            "L_J": index_J.leg(landing),
            "L_DJ": index_DJ.leg(landing),
            "L_JD": index_JD.leg(landing),
        }
        return {
            "supports": {
                "D0": support_D,
                "J0": support_J,
                "J_after_D": support_J_after_D,
                "D_after_J": support_D_after_J,
            },
            "fields": {
                "D": tuple(field_D),
                "J": tuple(field_J),
                "DJ": tuple(field_DJ),
                "JD": tuple(field_JD),
            },
            "responses": responses,
            "R_plus": (
                responses["L_D"]
                + responses["L_J"]
                - responses["L_DJ"]
                - responses["L_JD"]
            ),
            "R_minus": responses["L_DJ"] - responses["L_JD"],
        }


def empty_batch(batch: int, first_replica: int, samples: int) -> dict:
    row = {
        "batch": batch,
        "first_replica": first_replica,
        "samples": samples,
        "defined": 0,
        "sum_Rminus": 0,
        "sum_Rminus2": 0,
        "sum_Rplus": 0,
        "sum_Rplus2": 0,
        "nonzero_Rminus": 0,
        "typed_defined_mismatch": 0,
        "typed_Rminus_residual_max": 0,
        "typed_Rplus_sum_max": 0,
        "typed_support_mismatch": 0,
        "fixed_support_order_null_failures": 0,
    }
    for orientation in range(8):
        row[f"o{orientation}_samples"] = 0
        row[f"o{orientation}_defined"] = 0
        row[f"o{orientation}_sum_Rminus"] = 0
        row[f"o{orientation}_sum_Rplus"] = 0
    for value in range(-2, 3):
        row[f"hist_Rminus_{value}"] = 0
    return row


def run_batch(task: tuple[str, int, int, int, float, int]) -> dict:
    geometry_id, batch, first_replica, samples, probability, seed = task
    oracle = MediumOracle(geometry_id)
    row = empty_batch(batch, first_replica, samples)
    for replica in range(first_replica, first_replica + samples):
        field = [
            counter_uniform(seed, replica, site, oracle.geometry.n) < probability
            for site in range(oracle.geometry.n)
        ]
        anchor_D, anchor_J, landing, orientation = oracle.marks(seed, replica)
        row[f"o{orientation}_samples"] += 1
        primary = oracle.rectangle(
            field, anchor_D, anchor_J, landing, hand=False
        )
        dual = oracle.rectangle(
            [not value for value in field],
            anchor_J,
            anchor_D,
            landing,
            hand=True,
        )
        if (primary is None) != (dual is None):
            row["typed_defined_mismatch"] += 1
            continue
        if primary is None:
            continue

        row["defined"] += 1
        row[f"o{orientation}_defined"] += 1
        rminus = int(primary["R_minus"])
        rplus = int(primary["R_plus"])
        row["sum_Rminus"] += rminus
        row["sum_Rminus2"] += rminus * rminus
        row["sum_Rplus"] += rplus
        row["sum_Rplus2"] += rplus * rplus
        row["nonzero_Rminus"] += int(rminus != 0)
        row[f"o{orientation}_sum_Rminus"] += rminus
        row[f"o{orientation}_sum_Rplus"] += rplus
        row[f"hist_Rminus_{rminus}"] += 1

        row["typed_Rminus_residual_max"] = max(
            row["typed_Rminus_residual_max"], abs(rminus - int(dual["R_minus"]))
        )
        row["typed_Rplus_sum_max"] = max(
            row["typed_Rplus_sum_max"], abs(rplus + int(dual["R_plus"]))
        )
        primary_supports = primary["supports"]
        dual_supports = dual["supports"]
        support_pairs = (
            ("D0", "J0"),
            ("J0", "D0"),
            ("J_after_D", "D_after_J"),
            ("D_after_J", "J_after_D"),
        )
        if any(
            primary_supports[left]["site"] != dual_supports[right]["site"]
            for left, right in support_pairs
        ):
            row["typed_support_mismatch"] += 1

        # Freeze the two base supports.  Their ordered final fields must agree;
        # this is the exact fixed-support order-null control.
        fixed_DJ = oracle.flipped(
            oracle.flipped(field, primary_supports["D0"]["site"], False),
            primary_supports["J0"]["site"],
            True,
        )
        fixed_JD = oracle.flipped(
            oracle.flipped(field, primary_supports["J0"]["site"], True),
            primary_supports["D0"]["site"],
            False,
        )
        if fixed_DJ != fixed_JD:
            row["fixed_support_order_null_failures"] += 1
    return row


def mean_se(values: Sequence[float]) -> dict:
    mean = statistics.fmean(values)
    standard_error = (
        statistics.stdev(values) / sqrt(len(values)) if len(values) > 1 else 0.0
    )
    return {
        "value": mean,
        "standard_error": standard_error,
        "z": mean / standard_error if standard_error else None,
    }


def ratio_by_batch(rows: Sequence[dict], numerator: str, denominator: str) -> dict:
    values = [
        row[numerator] / row[denominator]
        for row in rows
        if row[denominator]
    ]
    return mean_se(values) | {"batches_used": len(values)}


def projection_summary(rows: Sequence[dict], source: str) -> dict:
    batch_modes: dict[str, list[float]] = {
        "scalar": [],
        "rotation_q1_re": [],
        "rotation_q1_im": [],
        "rotation_q2": [],
        "reflection_odd": [],
    }
    for row in rows:
        values = []
        for orientation in range(8):
            denominator = row[f"o{orientation}_samples"]
            numerator = (
                row[f"o{orientation}_sum_Rminus"]
                if source == "Rminus"
                else row[f"o{orientation}_defined"]
            )
            values.append(numerator / denominator)
        rotation_averages = [0.5 * (values[k] + values[k + 4]) for k in range(4)]
        reflection_differences = [
            0.5 * (values[k] - values[k + 4]) for k in range(4)
        ]
        q1 = sum(
            rotation_averages[k] * complex(cos(-2 * pi * k / 4), sin(-2 * pi * k / 4))
            for k in range(4)
        ) / 4
        q2 = sum(rotation_averages[k] * ((-1) ** k) for k in range(4)) / 4
        batch_modes["scalar"].append(sum(rotation_averages) / 4)
        batch_modes["rotation_q1_re"].append(q1.real)
        batch_modes["rotation_q1_im"].append(q1.imag)
        batch_modes["rotation_q2"].append(q2)
        batch_modes["reflection_odd"].append(sum(reflection_differences) / 4)
    return {name: mean_se(values) for name, values in batch_modes.items()}


def summarize(
    geometry_id: str,
    rows: Sequence[dict],
    probability: float,
    seed: int,
    replica_offset: int,
    elapsed: float,
) -> dict:
    total_samples = sum(row["samples"] for row in rows)
    total_defined = sum(row["defined"] for row in rows)
    histogram = {
        str(value): sum(row[f"hist_Rminus_{value}"] for row in rows)
        for value in range(-2, 3)
    }
    return {
        "schema": "matching-one/p250-adaptive-order-medium-pilot/v1",
        "status": "completed_intervention_sensitive_pilot",
        "geometry": {
            "id": geometry_id,
            **GEOMETRIES[geometry_id],
        },
        "run": {
            "samples": total_samples,
            "batches": len(rows),
            "p": probability,
            "seed": seed,
            "replica_offset": replica_offset,
            "elapsed_seconds": elapsed,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        },
        "estimands": {
            "defined_probability": ratio_by_batch(rows, "defined", "samples"),
            "Rminus_unconditional": ratio_by_batch(rows, "sum_Rminus", "samples"),
            "Rminus_given_defined": ratio_by_batch(rows, "sum_Rminus", "defined"),
            "Rplus_unconditional": ratio_by_batch(rows, "sum_Rplus", "samples"),
            "Rplus_given_defined": ratio_by_batch(rows, "sum_Rplus", "defined"),
            "nonzero_Rminus_given_defined": ratio_by_batch(
                rows, "nonzero_Rminus", "defined"
            ),
        },
        "counts": {
            "defined": total_defined,
            "Rminus_histogram_on_defined": histogram,
        },
        "typed_involution_controls": {
            "defined_mismatch": sum(row["typed_defined_mismatch"] for row in rows),
            "support_mismatch": sum(row["typed_support_mismatch"] for row in rows),
            "max_abs_Rminus_residual": max(row["typed_Rminus_residual_max"] for row in rows),
            "max_abs_Rplus_sum": max(row["typed_Rplus_sum_max"] for row in rows),
        },
        "fixed_support_order_null": {
            "failures": sum(row["fixed_support_order_null_failures"] for row in rows),
            "expected": "DJ and JD final fields coincide when D0/J0 supports are frozen",
        },
        "projection_leakage": {
            "Rminus_unconditional_D4_modes": projection_summary(rows, "Rminus"),
            "defined_probability_D4_modes": projection_summary(rows, "defined"),
            "interpretation": (
                "q1/q2 and reflection-odd modes are leakage controls; the scalar is the "
                "intervention response, not a CFT projection"
            ),
        },
        "claim_boundary": (
            "A nonzero result shows only that the frozen adaptive intervention remains "
            "order-sensitive on this finite Gaussian torus. It is not a CFT-field, "
            "universality, Jordan, or scaling-exponent claim."
        ),
    }


def write_batches(path: Path, rows: Sequence[dict]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", choices=tuple(GEOMETRIES), required=True)
    parser.add_argument("--samples", type=int, default=20_000)
    parser.add_argument("--batches", type=int, default=40)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--p", type=float, default=P_REF)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--replica-offset", type=int)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batches-output", type=Path)
    parser.add_argument("--git-commit", default="unknown")
    parser.add_argument("--environment", default="unknown")
    args = parser.parse_args()
    if args.samples <= 0 or args.batches < 8 or args.samples % args.batches:
        raise ValueError("samples must be positive and divisible by batches>=8")
    spec = GEOMETRIES[args.geometry]
    seed = int(spec["seed"] if args.seed is None else args.seed)
    replica_offset = int(
        spec["replica_offset"] if args.replica_offset is None else args.replica_offset
    )
    per_batch = args.samples // args.batches
    tasks = [
        (
            args.geometry,
            batch,
            replica_offset + batch * per_batch,
            per_batch,
            args.p,
            seed,
        )
        for batch in range(args.batches)
    ]
    start = time.perf_counter()
    if args.workers == 1:
        rows = [run_batch(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            rows = list(pool.map(run_batch, tasks))
    elapsed = time.perf_counter() - start
    rows.sort(key=lambda row: row["batch"])
    batches_path = args.batches_output or args.output.with_suffix(".batches.csv")
    batch_sha256 = write_batches(batches_path, rows)
    payload = summarize(
        args.geometry, rows, args.p, seed, replica_offset, elapsed
    )
    payload["run"].update(
        {
            "workers": args.workers,
            "git_commit": args.git_commit,
            "environment": args.environment,
            "batches_output": str(batches_path),
            "batches_sha256": batch_sha256,
        }
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
