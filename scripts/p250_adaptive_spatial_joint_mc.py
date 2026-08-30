#!/usr/bin/env python3
"""N505 adaptive-order response with a sampled Z/101 spatial coordinate."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

from p250_adaptive_order_pilot import MediumOracle, TypedIndex, transform_offset
from z5_projective_leg_cross_scale_mc import (
    CHILD_ORDER, PARENT_GEOMETRY, contexts, counter_uniform as cover_uniform,
)
from norm5_chiral_fixedp_mc import P_FIXED, splitmix64


GROUP_ORDER = 101
HANDS = ("plus", "minus")
RESIDUES = tuple(range(1, GROUP_ORDER))
SUPPORT_NAMES = ("D0", "J0", "J_after_D", "D_after_J")
RESPONSE_NAMES = ("L_D", "L_J", "L_DJ", "L_JD")
MASK64 = (1 << 64) - 1


def parent_residue(point: tuple[int, int]) -> int:
    return (point[0] - 10 * point[1]) % GROUP_ORDER


RESIDUE_COORDINATE = {
    parent_residue(tuple(point)): tuple(point) for point in PARENT_GEOMETRY.coordinates
}
if set(RESIDUE_COORDINATE) != set(range(GROUP_ORDER)):
    raise AssertionError("parent CRT residue map is not bijective")


class CoverOracle(MediumOracle):
    def __init__(self, hand: str):
        context = dict(zip(HANDS, contexts()))[hand]
        self.geometry_id = f"N505-{hand}"
        self.a, self.b = ((19, 12) if hand == "plus" else (21, -8))
        self.geometry = context.geometry
        self.neighbours = {
            "NN": self._neighbours(self.geometry.primal_edges),
            "matching": self._neighbours(self.geometry.matching_edges),
        }
        self.distance_from_origin = self._distance_table()
        self.field_to_vertex = tuple(
            self.geometry.vertex(point) for point in context.field_coordinates
        )
        inverse = [-1] * CHILD_ORDER
        for field_id, vertex in enumerate(self.field_to_vertex):
            if inverse[vertex] != -1:
                raise AssertionError("cover field labels are not injective")
            inverse[vertex] = field_id
        if any(value < 0 for value in inverse):
            raise AssertionError("cover field labels are not exhaustive")
        self.vertex_to_field = tuple(inverse)

    def active_from_field(self, field):
        active = [False] * CHILD_ORDER
        for field_id, value in enumerate(field):
            active[self.field_to_vertex[field_id]] = bool(value)
        return active


def support_with_reason(
    oracle: CoverOracle,
    index: TypedIndex,
    anchor_D: int,
    anchor_J: int,
    landing: int,
    operation: str,
):
    occupied = operation == "D"
    source_anchor = anchor_D if occupied else anchor_J
    target_anchor = anchor_J if occupied else anchor_D
    source = index.component(occupied, source_anchor)
    target = index.component(not occupied, target_anchor)
    if source is None or target is None:
        return None, "anchor_state"
    if source["rank"] < 1:
        return None, "source_rank0"
    target_vertices = target["vertex_set"]
    transferable = []
    for vertex in source["vertices"]:
        if vertex in (anchor_D, anchor_J):
            continue
        if any(neighbour in target_vertices
               for neighbour in oracle.neighbours[target["graph"]][vertex]):
            transferable.append(vertex)
    if not transferable:
        return None, "no_transfer"
    distances = [oracle.distance_squared(vertex, landing) for vertex in transferable]
    minimum = min(distances)
    minimizers = [vertex for vertex, distance in zip(transferable, distances)
                  if distance == minimum]
    if len(minimizers) != 1:
        return None, "tie"
    return {
        "site": minimizers[0],
        "site_field": oracle.vertex_to_field[minimizers[0]],
        "source_id": source["id"], "target_id": target["id"],
        "source_rank": source["rank"], "target_rank": target["rank"],
        "source_size": source["size"], "target_size": target["size"],
        "source_basis": source["basis"], "target_basis": target["basis"],
        "source_graph": source["graph"], "target_graph": target["graph"],
        "distance_squared": minimum, "minimizer_count": 1,
    }, "defined"


def rectangle_with_reason(
    oracle: CoverOracle, active, anchor_D: int, anchor_J: int, landing: int,
    *, hand: bool,
):
    if anchor_D == anchor_J or anchor_D == landing or anchor_J == landing:
        return None, "marked_collision"
    if not active[anchor_D] or active[anchor_J]:
        return None, "anchor_state"
    base = TypedIndex(oracle, active, hand)
    support_D, reason_D = support_with_reason(
        oracle, base, anchor_D, anchor_J, landing, "D")
    support_J, reason_J = support_with_reason(
        oracle, base, anchor_D, anchor_J, landing, "J")
    if support_D is None or support_J is None:
        reasons = (reason_D, reason_J)
        return None, "tie" if "tie" in reasons else "base_support"
    field_D = oracle.flipped(active, support_D["site"], False)
    field_J = oracle.flipped(active, support_J["site"], True)
    index_D, index_J = TypedIndex(oracle, field_D, hand), TypedIndex(oracle, field_J, hand)
    support_JD, reason_JD = support_with_reason(
        oracle, index_D, anchor_D, anchor_J, landing, "J")
    support_DJ, reason_DJ = support_with_reason(
        oracle, index_J, anchor_D, anchor_J, landing, "D")
    if support_JD is None or support_DJ is None:
        reasons = (reason_JD, reason_DJ)
        return None, "tie" if "tie" in reasons else "ordered_support"
    field_DJ = oracle.flipped(field_D, support_JD["site"], True)
    field_JD = oracle.flipped(field_J, support_DJ["site"], False)
    responses = {
        "L_D": index_D.leg(landing), "L_J": index_J.leg(landing),
        "L_DJ": TypedIndex(oracle, field_DJ, hand).leg(landing),
        "L_JD": TypedIndex(oracle, field_JD, hand).leg(landing),
    }
    payload = {
        "supports": {"D0": support_D, "J0": support_J,
                     "J_after_D": support_JD, "D_after_J": support_DJ},
        "responses": responses,
        "R_plus": responses["L_D"] + responses["L_J"] -
                  responses["L_DJ"] - responses["L_JD"],
        "R_minus": responses["L_DJ"] - responses["L_JD"],
    }
    return payload, "defined"


def marks(oracles, seed: int, replica: int):
    translation = splitmix64(seed ^ splitmix64(replica ^ 0x250505ADA001)) % GROUP_ORDER
    fiber = splitmix64(seed ^ splitmix64(replica ^ 0x250505ADA002)) % 5
    raw_residue = 1 + splitmix64(seed ^ splitmix64(replica ^ 0x250505ADA003)) % 100
    orientation = replica % 4
    tcoord = tuple(PARENT_GEOMETRY.coordinates[translation])
    displacement = transform_offset(orientation, RESIDUE_COORDINATE[raw_residue])
    actual_residue = parent_residue(displacement)
    target_parent = PARENT_GEOMETRY.vertex(
        (tcoord[0] + displacement[0], tcoord[1] + displacement[1]))
    output = {}
    for hand, oracle in oracles.items():
        anchor_D = oracle.field_to_vertex[5 * translation + fiber]
        ax, ay = oracle.geometry.coordinates[anchor_D]
        j_offset = transform_offset(orientation, (1, 1))
        anchor_J = oracle.geometry.vertex((ax + j_offset[0], ay + j_offset[1]))
        landing = oracle.field_to_vertex[5 * target_parent + fiber]
        output[hand] = (anchor_D, anchor_J, landing)
    return output, int(actual_residue), int(translation), int(fiber), orientation


def empty_batch(batch: int, first_replica: int, samples: int):
    row = {"batch": batch, "first_replica": first_replica, "samples": samples,
           "typed_defined_mismatch": 0, "typed_support_mismatch": 0,
           "typed_Rminus_residual_max": 0, "typed_Rplus_sum_max": 0}
    for hand in HANDS:
        row[f"{hand}_defined"] = 0; row[f"{hand}_ties"] = 0
        row[f"{hand}_sum_Rminus"] = 0; row[f"{hand}_sum_Rplus"] = 0
        for residue in RESIDUES:
            prefix = f"{hand}_j{residue}_"
            for field in ("samples", "defined", "ties", "sum_Rminus", "sum_Rplus",
                          "sum_L_D", "sum_L_J", "sum_L_DJ", "sum_L_JD"):
                row[prefix + field] = 0
    return row


def support_raw(prefix: str, support: dict, row: dict):
    row[f"{prefix}_site_field"] = support["site_field"]
    row[f"{prefix}_source_id"] = support["source_id"]
    row[f"{prefix}_target_id"] = support["target_id"]
    row[f"{prefix}_source_rank"] = support["source_rank"]
    row[f"{prefix}_target_rank"] = support["target_rank"]
    row[f"{prefix}_source_basis"] = repr(support["source_basis"])
    row[f"{prefix}_target_basis"] = repr(support["target_basis"])


def run_batch(task):
    batch, start, samples, p, seed = task
    oracles = {hand: CoverOracle(hand) for hand in HANDS}
    row = empty_batch(batch, start, samples)
    raw = []
    field_digest, mark_digest = hashlib.sha256(), hashlib.sha256()
    for replica in range(start, start + samples):
        field = [cover_uniform(seed, replica, site) < p for site in range(CHILD_ORDER)]
        field_digest.update(bytes(field))
        marked, residue, translation, fiber, orientation = marks(oracles, seed, replica)
        mark_digest.update(bytes((residue, translation, fiber, orientation)))
        for hand, oracle in oracles.items():
            prefix = f"{hand}_j{residue}_"
            row[prefix + "samples"] += 1
            active = oracle.active_from_field(field)
            anchor_D, anchor_J, landing = marked[hand]
            primary, reason = rectangle_with_reason(
                oracle, active, anchor_D, anchor_J, landing, hand=False)
            dual, dual_reason = rectangle_with_reason(
                oracle, [not value for value in active], anchor_J, anchor_D,
                landing, hand=True)
            if (primary is None) != (dual is None):
                row["typed_defined_mismatch"] += 1
            if primary is None:
                if reason == "tie":
                    row[f"{hand}_ties"] += 1; row[prefix + "ties"] += 1
                continue
            row[f"{hand}_defined"] += 1; row[prefix + "defined"] += 1
            for name in RESPONSE_NAMES:
                row[prefix + "sum_" + name] += primary["responses"][name]
            row[f"{hand}_sum_Rminus"] += primary["R_minus"]
            row[f"{hand}_sum_Rplus"] += primary["R_plus"]
            row[prefix + "sum_Rminus"] += primary["R_minus"]
            row[prefix + "sum_Rplus"] += primary["R_plus"]
            if dual is not None:
                row["typed_Rminus_residual_max"] = max(
                    row["typed_Rminus_residual_max"],
                    abs(primary["R_minus"] - dual["R_minus"]))
                row["typed_Rplus_sum_max"] = max(
                    row["typed_Rplus_sum_max"],
                    abs(primary["R_plus"] + dual["R_plus"]))
                pairs = (("D0", "J0"), ("J0", "D0"),
                         ("J_after_D", "D_after_J"),
                         ("D_after_J", "J_after_D"))
                row["typed_support_mismatch"] += int(any(
                    primary["supports"][left]["site_field"] !=
                    dual["supports"][right]["site_field"]
                    for left, right in pairs))
            raw_row = {"batch": batch, "replica": replica, "hand": hand,
                       "residue": residue, "translation": translation,
                       "fiber": fiber, "orientation": orientation,
                       "anchor_D_field": oracle.vertex_to_field[anchor_D],
                       "anchor_J_field": oracle.vertex_to_field[anchor_J],
                       "landing_field": oracle.vertex_to_field[landing],
                       **primary["responses"], "R_plus": primary["R_plus"],
                       "R_minus": primary["R_minus"], "dual_reason": dual_reason}
            for name in SUPPORT_NAMES:
                support_raw(name, primary["supports"][name], raw_row)
                if dual is not None:
                    support_raw("dual_" + name, dual["supports"][name], raw_row)
            raw.append(raw_row)
    row["field_sha256"] = field_digest.hexdigest()
    row["mark_sha256"] = mark_digest.hexdigest()
    return row, raw


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=20_000)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--p", type=float, default=P_FIXED)
    parser.add_argument("--seed", type=int, default=25050510120263000)
    parser.add_argument("--replica-offset", type=int, default=25_050_500_000)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batches-output", type=Path, required=True)
    parser.add_argument("--defined-output", type=Path, required=True)
    parser.add_argument("--git-commit", default="unknown")
    parser.add_argument("--environment", default="unknown")
    args = parser.parse_args()
    if args.samples <= 0 or args.batches < 10 or args.samples % args.batches:
        raise ValueError("samples must be positive and divisible by batches>=10")
    per_batch = args.samples // args.batches
    tasks = [(batch, args.replica_offset + batch * per_batch, per_batch,
              args.p, args.seed) for batch in range(args.batches)]
    started = time.perf_counter()
    if args.workers == 1:
        results = [run_batch(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(run_batch, tasks))
    results.sort(key=lambda value: value[0]["batch"])
    batch_rows = [value[0] for value in results]
    raw_rows = [row for value in results for row in value[1]]
    batch_hash = write_csv(args.batches_output, batch_rows)
    raw_hash = write_csv(args.defined_output, raw_rows)
    payload = {
        "schema": "matching-one/p250-adaptive-spatial-joint/v1",
        "run": {"samples": args.samples, "batches": args.batches,
                "workers": args.workers, "p": args.p, "seed": args.seed,
                "replica_offset": args.replica_offset,
                "elapsed_seconds": time.perf_counter() - started,
                "python": sys.version.split()[0], "platform": platform.platform(),
                "git_commit": args.git_commit, "environment": args.environment},
        "spatial_group": "Z/101, j(a,b)=a-10b mod 101; nonzero residues sampled uniformly",
        "hands": list(HANDS), "defined_rows": len(raw_rows),
        "batch_sha256": batch_hash, "defined_sha256": raw_hash,
        "typed_controls": {
            "defined_mismatch": sum(row["typed_defined_mismatch"] for row in batch_rows),
            "support_mismatch": sum(row["typed_support_mismatch"] for row in batch_rows),
            "max_Rminus_residual": max(row["typed_Rminus_residual_max"] for row in batch_rows),
            "max_Rplus_sum": max(row["typed_Rplus_sum_max"] for row in batch_rows)},
        "claim_boundary": "adaptive intervention response on two N505 children; not a spontaneous path-memory or field-identity claim"
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
