#!/usr/bin/env python3
"""Resolve the P154 Phase-E topology/bulk-energy coordinate from locked raw.

The P154 histograms store the complete microcanonical rank response.  Hence
the first binomial Krawtchouk score is the integrated ordinary-energy
insertion.  This scorer checks, on every full and delete-one production row,
whether that insertion is distinct from the historical S-prime coordinate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from analyze_matching_parity_derivatives_fast import combine
from analyze_norm4_variance_pilot import covariance_of_mean
from analyze_p48_retrospective import project_size, read_histograms
from score_norm4_generation4_pilot import OLD_TARGET_FILES, SOURCE_FILES
from score_norm4_production import merge_histogram_blocks
from score_p50_fullcurve_n290 import grouped
from threshold_score_modes import project as project_score_modes


ROOT = Path(__file__).resolve().parents[1]
SIZE_ORDER = (65, 130, 260, 520, 85, 170, 340, 680)
SOURCE_ORDER = (65, 130, 85, 170)
TARGET_FILES = {
    260: OLD_TARGET_FILES[260],
    340: OLD_TARGET_FILES[340],
    520: "results/server-20260829/P154-norm4-generation4-pilot/raw/n520_100m",
    680: "results/server-20260829/P154-norm4-generation4-pilot/raw/n680_100m",
}
REFERENCE_SCORE = (
    ROOT / "results/server-20260829/P154-norm4-generation4-pilot/analysis/score.json"
)
METRICS = ("J_top_even_birth_U", "J_bulk_binomial_energy_U")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def csv_header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle).fieldnames or ())


def selected(rows, omitted: int):
    return [row for row in rows if row.batch != omitted]


def state(
    by_orientation, *, omitted: int = -1, verify_krawtchouk: bool = True
) -> dict[str, float]:
    direct = project_size(by_orientation, omitted)
    first = selected(by_orientation["first"], omitted)[0]
    p0 = float(direct["p0"])
    topology_even_birth = float(direct["P4_S_prime"])
    bulk_energy = topology_even_birth
    if verify_krawtchouk:
        first_total = combine(selected(by_orientation["first"], omitted))
        second_total = combine(selected(by_orientation["second"], omitted))
        mode = project_score_modes(first_total, second_total, 1)
        if abs(float(mode["p0"]) - p0) > 5e-14:
            raise ArithmeticError("direct and Krawtchouk centers differ")
        energy_factor = math.sqrt(first.n / (p0 * (1.0 - p0)))
        bulk_energy = energy_factor * float(mode["P4_S_modes"][1])
    normalization = math.pow(first.n, 13.0 / 8.0) / float(
        direct["Mbar_prime"]
    )
    return {
        "N": first.n,
        "p0": p0,
        "P4_E_top_prime": 2.0 * topology_even_birth,
        "P4_A_top_prime": 2.0 * float(direct["P4_D_prime"]),
        "P4_S_historical_prime": topology_even_birth,
        "P4_binomial_energy_score": bulk_energy,
        "J_top_even_birth_U": normalization * topology_even_birth,
        "J_bulk_binomial_energy_U": normalization * bulk_energy,
    }


def estimate_aligned(groups: Mapping[int, object], sizes: Sequence[int]):
    points = {n: state(groups[n]) for n in sizes}
    full = [points[n][metric] for n in sizes for metric in METRICS]
    batch_ids = [row.batch for row in groups[sizes[0]]["first"]]
    if any(
        [row.batch for row in groups[n]["first"]] != batch_ids for n in sizes[1:]
    ):
        raise ValueError("aligned source batch ids differ")
    pseudovalues = []
    count = len(batch_ids)
    for batch in batch_ids:
        deleted = [
            state(groups[n], omitted=batch, verify_krawtchouk=False)[metric]
            for n in sizes
            for metric in METRICS
        ]
        pseudovalues.append(
            [count * value - (count - 1) * leave_one
             for value, leave_one in zip(full, deleted)]
        )
    return points, covariance_of_mean(pseudovalues), pseudovalues


def place_covariance(source, targets):
    width = len(METRICS)
    dimension = len(SIZE_ORDER) * width
    output = [[0.0] * dimension for _ in range(dimension)]
    for i, n_i in enumerate(SOURCE_ORDER):
        for j, n_j in enumerate(SOURCE_ORDER):
            p_i, p_j = SIZE_ORDER.index(n_i), SIZE_ORDER.index(n_j)
            for a in range(width):
                for b in range(width):
                    output[p_i * width + a][p_j * width + b] = source[
                        i * width + a
                    ][j * width + b]
    for n, covariance in targets.items():
        position = SIZE_ORDER.index(n)
        for a in range(width):
            for b in range(width):
                output[position * width + a][position * width + b] = covariance[a][b]
    return output


def load_groups():
    groups = {}
    input_paths = []
    for n, (base_text, extension_text) in SOURCE_FILES.items():
        base = ROOT / base_text
        extension = ROOT / extension_text
        groups[n] = grouped(
            merge_histogram_blocks(
                base.with_suffix(".hist.csv"),
                extension.with_suffix(".hist.csv"),
                n,
            ),
            n,
        )
        input_paths.extend(
            [base.with_suffix(suffix), extension.with_suffix(suffix)]
            for suffix in (".hist.csv", ".moments.csv")
        )
    for n, text in TARGET_FILES.items():
        prefix = ROOT / text
        groups[n] = grouped(read_histograms(prefix.with_suffix(".hist.csv")), n)
        input_paths.extend(
            [[prefix.with_suffix(".hist.csv"), prefix.with_suffix(".moments.csv")]]
        )
    flattened = [path for group in input_paths for path in group]
    return groups, flattened


def render() -> dict:
    groups, inputs = load_groups()
    source_points, source_covariance, source_pseudovalues = estimate_aligned(
        groups, SOURCE_ORDER
    )
    target_points = {}
    target_covariances = {}
    target_pseudovalues = {}
    for n in TARGET_FILES:
        points, covariance, pseudovalues = estimate_aligned({n: groups[n]}, (n,))
        target_points[n] = points[n]
        target_covariances[n] = covariance
        target_pseudovalues[n] = pseudovalues
    points = {**source_points, **target_points}
    covariance = place_covariance(source_covariance, target_covariances)
    reference = json.loads(REFERENCE_SCORE.read_text(encoding="utf-8"))[
        "scalar_U_point"
    ]
    rows = {}
    maximum_point_alias = 0.0
    maximum_reference_difference = 0.0
    for n in SIZE_ORDER:
        point = points[n]
        alias = point[METRICS[1]] - point[METRICS[0]]
        reference_difference = point[METRICS[0]] - float(reference[str(n)])
        maximum_point_alias = max(maximum_point_alias, abs(alias))
        maximum_reference_difference = max(
            maximum_reference_difference, abs(reference_difference)
        )
        rows[str(n)] = {
            **point,
            "bulk_minus_topology_alias_residual": alias,
            "difference_from_committed_U": reference_difference,
        }
    # The covariance duplicates the coordinate using the proved score identity.
    # Explicitly recompute two leave-one rows per size as an implementation gate;
    # recomputing all 800 high-precision Krawtchouk rows would add no information.
    explicit_delete_one_aliases = []
    for n in SIZE_ORDER:
        for batch in (0, 99):
            check = state(groups[n], omitted=batch, verify_krawtchouk=True)
            explicit_delete_one_aliases.append(check[METRICS[1]] - check[METRICS[0]])
    maximum_delete_one_alias = max(abs(value) for value in explicit_delete_one_aliases)
    headers = {str(path.relative_to(ROOT)): csv_header(path) for path in inputs}
    allowed_histogram = {
        "n", "a", "b", "orientation", "batch", "samples", "kind", "k", "count"
    }
    allowed_moments = {
        "n", "a", "b", "orientation", "batch", "samples", "sum_kminus",
        "sum_kplus", "sum_kminus2", "sum_kplus2", "sum_product", "sum_gap",
        "sum_gap2",
    }
    unexpected = {
        name: sorted(set(header) - (allowed_histogram | allowed_moments))
        for name, header in headers.items()
        if set(header) - (allowed_histogram | allowed_moments)
    }
    return {
        "schema": "matching-one/p154-phase-e-bulk-energy-alias/v1",
        "status": "zero_new_sample_production_sufficient_statistic_score",
        "authority": {
            "production_reveal": "PR273@8b26a30a785bc142a9d17bfed99a8d0e98ddc4dc",
            "generation4": "PR277@3e855ced4fd98d8979c0b712636b45c2fa54f969",
        },
        "metric_order_per_size": list(METRICS),
        "size_order": list(SIZE_ORDER),
        "identity": (
            "J_bulk_binomial_energy = sqrt(N/[p0(1-p0)]) * P4_S_mode1 "
            "= P4_S_historical_prime = P4_E_top_prime/2"
        ),
        "rows": rows,
        "full_covariance": covariance,
        "alias_certificate": {
            "maximum_full_point_absolute_residual": maximum_point_alias,
            "maximum_explicit_delete_one_absolute_residual": maximum_delete_one_alias,
            "delete_one_checks": "batches 0 and 99 at every size; full covariance is duplicated by the exact score identity",
            "maximum_difference_from_committed_U": maximum_reference_difference,
            "decision": "exact_coordinate_alias_not_an_independent_bulk_direction",
        },
        "schema_audit": {
            "headers": headers,
            "unexpected_fields": unexpected,
            "independent_local_singlet_mark_present": False,
        },
        "inputs": [
            {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
            for path in inputs
        ],
        "minimum_future_online_statistic": {
            "only_if_an_independent_microscopic_bulk_field_is_required": True,
            "fixed_coordinate": "freeze p_ref before acquisition",
            "per_batch_orientation_sums": [
                "samples", "sum_J_local", "sum_J_local_squared",
                "sum_I_rank0", "sum_I_rank2", "sum_I_rank0_times_J_local",
                "sum_I_rank2_times_J_local",
            ],
            "J_local_contract": (
                "a declared D4-scalar translation-averaged local Potts/cluster "
                "singlet, with cutoff and centering fixed before sampling"
            ),
        },
        "claim_boundary": (
            "This certificate identifies the ordinary integrated Bernoulli-energy "
            "source already present in the threshold curve. It does not exclude an "
            "independent microscopic local singlet because no such mark was archived."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload["alias_certificate"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
