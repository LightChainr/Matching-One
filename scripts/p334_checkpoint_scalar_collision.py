#!/usr/bin/env python3
"""Zero-new-sample real-lattice (H2,b2) branching nonclosure witnesses."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MASK = (1 << 64) - 1
STEP = 0x9e3779b97f4a7c15


def mix(value):
    value = (value + STEP) & MASK
    value = ((value ^ (value >> 30))*0xbf58476d1ce4e5b9) & MASK
    value = ((value ^ (value >> 27))*0x94d049bb133111eb) & MASK
    return value ^ (value >> 31)


def archived_permutation(n, seed, replica):
    """Decode the archived runner's exact uint64 counter; no new sampling."""
    state = mix(seed ^ mix((replica + 0xd1b54a32d192ed03) & MASK))
    permutation = list(range(n))
    for stop in range(n-1, 0, -1):
        bound = stop+1
        remainder = ((MASK % bound)+1) % bound
        while True:
            value = mix(state)
            state = (state+STEP) & MASK
            if remainder == 0 or value <= MASK-remainder:
                other = value % bound
                break
        permutation[stop], permutation[other] = permutation[other], permutation[stop]
    return permutation


def fraction(value):
    return {"numerator": value.numerator, "denominator": value.denominator,
            "text": str(value), "decimal": float(value)}


def grouping(rows, fields):
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[field] for field in fields)].append(row)
    collisions = [(key, value) for key, value in groups.items()
                  if len({r["checkpoint_sum_child_b1_sq"] for r in value}) > 1]
    return groups, collisions


def witness_row(row, metadata, matrix):
    n, k0 = row["n"], row["k0"]
    d = n-k0
    b1, b2 = row["checkpoint_b1_safe_count"], row["checkpoint_b2_safe_pairs"]
    squares = row["checkpoint_sum_child_b1_sq"]
    numerator = b1*squares-(2*b2)**2
    permutation = archived_permutation(n, metadata["seed"], row["replica"])
    if permutation[k0] != row["next_site"]:
        raise ValueError("decoded archive counter disagrees with saved next site")
    return {
        "original_checkpoint_row": row,
        "replay": {"seed": metadata["seed"], "replica_counter": row["replica"],
                   "runner_commit": metadata["git_commit"], "k0": k0,
                   "period_matrix": matrix,
                   "coordinate_rule": metadata["quotient_coordinates"],
                   "occupied_prefix_labels": permutation[:k0],
                   "occupied_mask_hex": hex(sum(1 << v for v in permutation[:k0])),
                   "next_site_decode_matches": True},
        "safe_insertion_graph": {
            "safe_singletons_b1": b1, "edges_b2": b2,
            "degree_sum": 2*b2, "degree_square_sum": squares,
            "unordered_2stars": (squares-2*b2)//2,
        },
        "branch_success": fraction(Fraction(squares, d*(d-1)**2)),
        "delta_variance_integer_numerator": numerator,
        "delta_coop": fraction(Fraction(numerator, d*b1*(d-1)**2)) if b1 else fraction(Fraction(0)),
    }


def build():
    lock = json.loads((ROOT/"analysis/p334_cooperative_closure_raw_lock.json").read_text())
    coarse_fields = ("n", "orientation", "k0", "H2", "checkpoint_b2_safe_pairs")
    strong_fields = coarse_fields + ("age_steps", "ell_u", "ell_v")
    environments = {}
    for size, source in lock["runs"].items():
        prefix = ROOT/source["prefix"]
        csv_path = Path(str(prefix)+".geometry_pilot.csv")
        csv_hash = hashlib.sha256(csv_path.read_bytes()).hexdigest()
        if csv_hash != source["csv_sha256"]:
            raise ValueError("raw archive differs from recorded production")
        metadata = json.loads(Path(str(prefix)+".metadata.json").read_text())
        with csv_path.open(newline="") as stream:
            rows = [{k: (v if k == "orientation" else float(v) if k in ("q_after", "q_after2") else int(v))
                     for k, v in row.items()} for row in csv.DictReader(stream)]
        for orientation in ("first", "second"):
            selected = [r for r in rows if r["orientation"] == orientation]
            groups, collisions = grouping(selected, coarse_fields)
            strong_groups, strong_collisions = grouping(selected, strong_fields)
            item = {"source_csv_sha256": csv_hash, "archived_rows": len(selected),
                    "coarse_groups": len(groups), "coarse_collision_groups": len(collisions),
                    "coarse_collision_rows": sum(len(v) for _, v in collisions),
                    "age_line_matched_groups": len(strong_groups),
                    "age_line_matched_collision_groups": len(strong_collisions),
                    "age_line_matched_collision_rows": sum(len(v) for _, v in strong_collisions)}
            if strong_collisions:
                # Post-hoc exact witness choice, not an effect-size hypothesis test.
                key, group = max(strong_collisions, key=lambda kv: (
                    max(r["checkpoint_sum_child_b1_sq"] for r in kv[1])-
                    min(r["checkpoint_sum_child_b1_sq"] for r in kv[1]), kv[0]))
                low = min(group, key=lambda r: (r["checkpoint_sum_child_b1_sq"], r["replica"]))
                high = max(group, key=lambda r: (r["checkpoint_sum_child_b1_sq"], -r["replica"]))
                matrix = metadata["designs"][0][f"{orientation}_period_matrix"]
                first, second = (witness_row(r, metadata, matrix) for r in (low, high))
                difference = Fraction(second["branch_success"]["numerator"], second["branch_success"]["denominator"])-Fraction(first["branch_success"]["numerator"], first["branch_success"]["denominator"])
                item["witness"] = {
                    "identical_fields": dict(zip(strong_fields, key)),
                    "A": first, "B": second,
                    "degree_square_difference": high["checkpoint_sum_child_b1_sq"]-low["checkpoint_sum_child_b1_sq"],
                    "branch_success_difference_B_minus_A": fraction(difference),
                    "delta_coop_difference_B_minus_A": fraction(difference),
                }
            environments[f"{size}_{orientation}"] = item
    return {
        "schema": "matching-one/p334-real-checkpoint-scalar-nonclosure/v1",
        "production_result_commit": "e81dd59ff6be69056e504e0e81cfeccf73dc5e97",
        "new_samples": 0,
        "group_fields": list(coarse_fields), "strong_group_fields": list(strong_fields),
        "statement": "In each of the four fixed N/orientation geometries, real archived checkpoints with identical H2,b2,k0,age and primitive line have different exact one-common-update/two-clone survival probabilities. Thus this scalar state, even augmented by age and line, is not sufficient for that branching prediction.",
        "graph_identity": "Safe singleton sites are vertices; safe unordered two-site sets are edges. Their degrees are c_v. Sum c_v=2*b2, sum c_v^2=2*b2+2*unordered_2stars. The branching observable is the degree second moment divided by d*(d-1)^2.",
        "selection": "Within each environment, choose the age-and-line-matched collision with largest integer degree-square range; deterministic tie-break. Selection is for an exact counterexample, not significance inference.",
        "environments": environments,
        "claim_boundary": "No geometries or orientations are pooled. This is insufficiency of the specified scalar state for the specified branching observable; it does not prove hidden temporal memory of the full lattice configuration, complete single-chain trace equivalence, or a continuum field.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT/"results/local-20260831/P334-cooperative-closure/scalar_state_collisions.json")
    args = parser.parse_args()
    payload = build()
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
    for key, value in payload["environments"].items():
        witness = value.get("witness")
        print(key, "coarse/strong_collision_groups", value["coarse_collision_groups"],
              value["age_line_matched_collision_groups"],
              "branch_difference", witness["branch_success_difference_B_minus_A"] if witness else None)


if __name__ == "__main__":
    main()
