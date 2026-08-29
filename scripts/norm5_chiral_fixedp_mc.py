#!/usr/bin/env python3
"""Deployable fixed-p N325 chiral Z5 response stream.

This is intentionally one design, not a general cover framework.  A common
Bernoulli field is indexed by the exact ``(parent site, Z5 fiber)`` labels for
the two children of ``8+i`` under ``2+i`` and ``2-i``.  One uniformly cycled
marked root per replica gives an unbiased estimator of the opposite-character
sum while retaining the fixed-root pivotal cost of #215/#225.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import json
from math import atan2, cos, pi, sin, sqrt
from pathlib import Path
import sys
from typing import Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))

from integer_period_torus import (  # noqa: E402
    IntegerTorusGeometry,
    gaussian_integer_torus,
    integer_torus_geometry,
    matrix_product,
    matrix_vector,
)
from marked_pivotal_h4_reference import marked_pair  # noqa: E402
from norm5_chiral_hecke_phase import gaussian_ratio_power  # noqa: E402


MASK64 = (1 << 64) - 1
P_FIXED = 0.592746050790
SMOKE_CAP = 20_000
PRODUCTION_ID = "P226-N325-norm5-chiral-fixedp-v1"
PRODUCTION_RUN = {
    "samples": 200_000,
    "batches": 100,
    "workers": 16,
    "p": P_FIXED,
    "seed": 2265325000829,
    "radius": 1,
    "output": "/workspace/Matching-One/results/huawei-20260829/P226-norm5-chiral-fixedp/chiral_response.json",
    "batches_output": "/workspace/Matching-One/results/huawei-20260829/P226-norm5-chiral-fixedp/chiral_response.batches.csv",
}

PARENT = (8, 1)
PARENT_MATRIX = ((8, -1), (1, 8))
M_PLUS = ((2, -1), (1, 2))
M_MINUS = ((2, 1), (-1, 2))
CHILD_PLUS = matrix_product(PARENT_MATRIX, M_PLUS)  # 15+10i
CHILD_MINUS = matrix_product(PARENT_MATRIX, M_MINUS)  # 17-6i
REFLECTION = ((1, 0), (0, -1))
PARENT_MIRROR = ((8, 1), (-1, 8))
CHILD_PLUS_MIRROR = matrix_product(PARENT_MIRROR, M_MINUS)  # 15-10i
PHASES = tuple((cos(2 * pi * k / 5), sin(2 * pi * k / 5)) for k in range(5))
PRIMARY_ORDER = ("plus_re", "plus_im", "minus_re", "minus_im")


def splitmix64(value: int) -> int:
    value = (value + 0x9E3779B97F4A7C15) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return (value ^ (value >> 31)) & MASK64


def counter_uniform(seed: int, replica: int, site: int) -> float:
    key = (
        seed
        ^ splitmix64(325)
        ^ splitmix64(replica + 0xD1B54A32D192ED03)
        ^ splitmix64(site + 0x94D049BB133111EB)
    )
    return (splitmix64(key) >> 11) * (2.0**-53)


def reflected(point: tuple[int, int]) -> tuple[int, int]:
    return point[0], -point[1]


class CoverContext:
    def __init__(
        self,
        periods,
        parent_matrix,
        parent_representatives: Sequence[tuple[int, int]],
        *,
        name: str,
    ) -> None:
        self.name = name
        self.geometry = integer_torus_geometry(periods, name=name)
        if self.geometry.n != 325 or len(parent_representatives) != 65:
            raise ValueError("frozen design requires a 65-to-325 cover")
        self.field_to_vertex: list[int] = []
        self.field_coordinates: list[tuple[int, int]] = []
        for representative in parent_representatives:
            for fiber in range(5):
                offset = matrix_vector(parent_matrix, (fiber, 0))
                point = representative[0] + offset[0], representative[1] + offset[1]
                self.field_coordinates.append(point)
                self.field_to_vertex.append(self.geometry.vertex(point))
        if len(set(self.field_to_vertex)) != 325:
            raise AssertionError(f"{name}: parent/fiber labels do not cover child")
        self.shift_permutations = [
            self._shift_to_origin(field) for field in range(325)
        ]

    def _shift_to_origin(self, root_field: int) -> list[int]:
        root_vertex = self.field_to_vertex[root_field]
        root = self.geometry.coordinates[root_vertex]
        permutation = []
        for point in self.geometry.coordinates:
            permutation.append(self.geometry.vertex((point[0] - root[0], point[1] - root[1])))
        if len(set(permutation)) != 325 or permutation[root_vertex] != self.geometry.vertex((0, 0)):
            raise AssertionError(f"{self.name}: invalid root translation")
        return permutation

    def active_from_field(self, field: Sequence[bool]) -> list[bool]:
        active = [False] * 325
        for index, value in enumerate(field):
            active[self.field_to_vertex[index]] = bool(value)
        return active

    def shifted_to_root(self, active: Sequence[bool], root_field: int) -> list[bool]:
        shifted = [False] * 325
        for old_vertex, new_vertex in enumerate(self.shift_permutations[root_field]):
            shifted[new_vertex] = bool(active[old_vertex])
        return shifted


_CONTEXTS: tuple[CoverContext, CoverContext, CoverContext] | None = None


def contexts() -> tuple[CoverContext, CoverContext, CoverContext]:
    global _CONTEXTS
    if _CONTEXTS is None:
        parent = gaussian_integer_torus(*PARENT)
        parent_representatives = parent.coordinates
        mirror_representatives = tuple(reflected(point) for point in parent_representatives)
        _CONTEXTS = (
            CoverContext(
                CHILD_PLUS,
                PARENT_MATRIX,
                parent_representatives,
                name="N325-parent8+i-times2+i",
            ),
            CoverContext(
                CHILD_MINUS,
                PARENT_MATRIX,
                parent_representatives,
                name="N325-parent8+i-times2-i",
            ),
            CoverContext(
                CHILD_PLUS_MIRROR,
                PARENT_MIRROR,
                mirror_representatives,
                name="N325-reflection-parent8-i-times2-i",
            ),
        )
    return _CONTEXTS


def mapping_gate() -> dict:
    plus, minus, mirror = contexts()
    parent = gaussian_integer_torus(*PARENT)
    for context in (plus, minus):
        for field_index, point in enumerate(context.field_coordinates):
            parent_index = field_index // 5
            if parent.vertex(point) != parent_index:
                raise AssertionError(f"{context.name}: child-to-parent projection failed")
    for field_index in range(325):
        if reflected(plus.field_coordinates[field_index]) != mirror.field_coordinates[field_index]:
            raise AssertionError("reflection did not preserve the common (parent,fiber) label")
    # The first parent period translates k -> k+1 for both multiplier quotients.
    deck_step = matrix_vector(PARENT_MATRIX, (1, 0))
    for context in (plus, minus):
        for parent_index in range(65):
            for fiber in range(5):
                field = 5 * parent_index + fiber
                translated = (
                    context.field_coordinates[field][0] + deck_step[0],
                    context.field_coordinates[field][1] + deck_step[1],
                )
                target = 5 * parent_index + (fiber + 1) % 5
                if context.geometry.vertex(translated) != context.field_to_vertex[target]:
                    raise AssertionError(f"{context.name}: deck step failed")
    return {
        "parent": [8, 1],
        "parent_order": 65,
        "plus_child_periods": CHILD_PLUS,
        "minus_child_periods": CHILD_MINUS,
        "child_order": 325,
        "common_field": "(deterministic parent representative j, fiber k in Z/5)",
        "fiber_representative": "u_k=(k,0) modulo M_(2+/-i)",
        "unique_labels_per_child": 325,
        "deck_step": list(deck_step),
        "deck_action": "k -> k+1 mod 5",
        "true_reflection_pair": "(8+i)(2+i) <-> (8-i)(2-i)",
        "same_parent_hands_are_reflections": False,
        "passed": True,
    }


def local_twice(context: CoverContext, field: Sequence[bool], root_field: int, radius: int) -> int:
    active = context.active_from_field(field)
    shifted = context.shifted_to_root(active, root_field)
    pair = marked_pair(context.geometry, shifted, radius)
    return pair["primal"]["h4"] - pair["matching"]["h4"]


def score(field: Sequence[bool], p: float) -> complex:
    value = 0j
    for index, state in enumerate(field):
        phase = PHASES[index % 5]
        centered = int(state) - p
        value += complex(phase[0], phase[1]) * centered
    return value


def marked_product(local: int, root_field: int, score_value: complex) -> complex:
    phase = PHASES[root_field % 5]
    opposite = complex(phase[0], -phase[1])
    marked_estimator = (325.0 / 2.0) * local * opposite
    return marked_estimator * score_value


def _run_batch(task: tuple[int, int, int, float, int, int]) -> dict:
    batch, start, samples, p, seed, radius = task
    plus, minus, mirror = contexts()
    sums = {
        "batch": batch,
        "start_replica": start,
        "samples": samples,
        "plus_re": 0.0,
        "plus_im": 0.0,
        "minus_re": 0.0,
        "minus_im": 0.0,
        "mirror_plus_re": 0.0,
        "mirror_plus_im": 0.0,
        "reflection_null_re": 0.0,
        "reflection_null_im": 0.0,
    }
    for replica in range(start, start + samples):
        field = [counter_uniform(seed, replica, site) < p for site in range(325)]
        root = replica % 325
        score_plus = score(field, p)
        local_plus = local_twice(plus, field, root, radius)
        local_minus = local_twice(minus, field, root, radius)
        local_mirror = local_twice(mirror, field, root, radius)
        plus_value = marked_product(local_plus, root, score_plus)
        minus_value = marked_product(local_minus, root, score_plus)
        # True reflection transports chi to chibar, so both score and marked
        # phase conjugate.  The local H4 mark is reflection even.
        mirror_value = marked_product(local_mirror, root, score_plus).conjugate()
        null = plus_value - mirror_value.conjugate()
        for prefix, value in (
            ("plus", plus_value),
            ("minus", minus_value),
            ("mirror_plus", mirror_value),
            ("reflection_null", null),
        ):
            sums[f"{prefix}_re"] += value.real
            sums[f"{prefix}_im"] += value.imag
    return sums


def covariance_of_mean(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(rows)
    width = len(rows[0])
    means = [sum(row[j] for row in rows) / count for j in range(width)]
    return [
        [
            sum((row[i] - means[i]) * (row[j] - means[j]) for row in rows)
            / (count * (count - 1))
            for j in range(width)
        ]
        for i in range(width)
    ]


def summarize(batches: Sequence[dict], p: float) -> dict:
    pq = p * (1 - p)
    rows = [
        [batch[name] / (batch["samples"] * pq) for name in PRIMARY_ORDER]
        for batch in batches
    ]
    point = [sum(row[j] for row in rows) / len(rows) for j in range(4)]
    covariance = covariance_of_mean(rows)
    plus = complex(point[0], point[1])
    minus = complex(point[2], point[3])
    ratio = plus / minus if minus else None

    total_samples = sum(batch["samples"] for batch in batches)
    total = [sum(batch[name] for batch in batches) for name in PRIMARY_ORDER]
    ratio_delete_one = []
    for batch in batches:
        kept = total_samples - batch["samples"]
        row = [(total[j] - batch[PRIMARY_ORDER[j]]) / (kept * pq) for j in range(4)]
        denominator = complex(row[2], row[3])
        if denominator:
            ratio_delete_one.append(complex(row[0], row[1]) / denominator)
    ratio_covariance = None
    phase_standard_error_degrees = None
    target_phase_differences = None
    if len(ratio_delete_one) == len(batches) and len(ratio_delete_one) >= 2:
        ratio_mean = sum(ratio_delete_one) / len(ratio_delete_one)
        factor = (len(ratio_delete_one) - 1) / len(ratio_delete_one)
        ratio_covariance = [
            [
                factor
                * sum(
                    ((value.real, value.imag)[i] - (ratio_mean.real, ratio_mean.imag)[i])
                    * ((value.real, value.imag)[j] - (ratio_mean.real, ratio_mean.imag)[j])
                    for value in ratio_delete_one
                )
                for j in range(2)
            ]
            for i in range(2)
        ]
        if ratio is not None and abs(ratio) > 0:
            radius_squared = ratio.real * ratio.real + ratio.imag * ratio.imag
            gradient = (-ratio.imag / radius_squared, ratio.real / radius_squared)
            phase_variance = sum(
                gradient[i] * ratio_covariance[i][j] * gradient[j]
                for i in range(2)
                for j in range(2)
            )
            phase_standard_error_degrees = sqrt(max(0.0, phase_variance)) * 180 / pi
            measured_phase = atan2(ratio.imag, ratio.real) * 180 / pi
            target_phase_differences = {}
            for spin in (4, 8, 12):
                real, imag, denominator = gaussian_ratio_power(spin)
                target = atan2(imag / denominator, real / denominator) * 180 / pi
                difference = (measured_phase - target + 180) % 360 - 180
                target_phase_differences[f"H{spin}"] = {
                    "target_phase_degrees": target,
                    "wrapped_difference_degrees": difference,
                }
    null_rows = [
        [
            batch["reflection_null_re"] / (batch["samples"] * pq),
            batch["reflection_null_im"] / (batch["samples"] * pq),
        ]
        for batch in batches
    ]
    null_point = [sum(row[j] for row in null_rows) / len(null_rows) for j in range(2)]
    null_covariance = covariance_of_mean(null_rows)
    return {
        "primary_order": list(PRIMARY_ORDER),
        "primary_point": point,
        "primary_standard_error": [sqrt(max(0.0, covariance[j][j])) for j in range(4)],
        "primary_covariance_of_mean": covariance,
        "same_parent_handed_ratio": {
            "point_re_im": [ratio.real, ratio.imag] if ratio is not None else None,
            "phase_radians": atan2(ratio.imag, ratio.real) if ratio is not None else None,
            "phase_degrees": atan2(ratio.imag, ratio.real) * 180 / pi if ratio is not None else None,
            "delete_one_covariance_re_im": ratio_covariance,
            "delete_one_valid_replicates": len(ratio_delete_one),
            "delta_method_phase_standard_error_degrees": phase_standard_error_degrees,
            "frozen_target_phase_differences": target_phase_differences,
        },
        "true_reflection_conjugacy_null": {
            "point_re_im": null_point,
            "covariance_of_mean": null_covariance,
            "relation": "R_((8+i)(2+i))-conj(R_((8-i)(2-i)))",
        },
        "same_parent_plus_minus_conjugacy_is_not_a_null": True,
    }


def validate_gate(args, batches_path: Path) -> str:
    if args.samples <= SMOKE_CAP:
        if args.production_manifest is not None:
            raise ValueError("production manifest cannot label a smoke run")
        return "engineering_smoke"
    if args.production_manifest is None:
        raise ValueError("runs above 20k require the frozen production manifest")
    manifest = json.loads(args.production_manifest.read_text())
    if manifest.get("schema") != "matching-one.norm5-chiral-fixedp-production.v1":
        raise ValueError("wrong production manifest schema")
    if manifest.get("production_id") != PRODUCTION_ID or manifest.get("run") != PRODUCTION_RUN:
        raise ValueError("manifest differs from the code-frozen production run")
    actual = {
        "samples": args.samples,
        "batches": args.batches,
        "workers": args.workers,
        "p": args.p,
        "seed": args.seed,
        "radius": args.radius,
        "output": str(args.output),
        "batches_output": str(batches_path),
    }
    if actual != PRODUCTION_RUN:
        raise ValueError("CLI differs from the frozen production run")
    return "production_under_frozen_manifest"


def run(samples: int, batches: int, workers: int, p: float, seed: int, radius: int):
    if samples <= 0 or batches < 2 or samples % batches:
        raise ValueError("samples must be divisible by batches>=2")
    per_batch = samples // batches
    tasks = [(batch, batch * per_batch, per_batch, p, seed, radius) for batch in range(batches)]
    if workers == 1:
        output = [_run_batch(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            output = list(pool.map(_run_batch, tasks))
    return output, summarize(output, p)


def write_batches(path: Path, batches: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(batches[0]))
        writer.writeheader()
        writer.writerows(batches)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=20_000)
    parser.add_argument("--batches", type=int, default=20)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--p", type=float, default=P_FIXED)
    parser.add_argument("--seed", type=int, default=2265325020000)
    parser.add_argument("--radius", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batches-output", type=Path)
    parser.add_argument("--production-manifest", type=Path)
    args = parser.parse_args()
    batches_path = args.batches_output or args.output.with_suffix(".batches.csv")
    status = validate_gate(args, batches_path)
    gate = mapping_gate()
    batches, analysis = run(
        args.samples, args.batches, args.workers, args.p, args.seed, args.radius
    )
    write_batches(batches_path, batches)
    payload = {
        "schema": "matching-one.norm5-chiral-fixedp-response.v1",
        "issues": [226, 244],
        "status": status,
        "mapping_gate": gate,
        "run": {
            "samples": args.samples,
            "batches": args.batches,
            "workers": args.workers,
            "p": args.p,
            "seed": args.seed,
            "radius": args.radius,
            "batches_output": str(batches_path),
        },
        "observable": {
            "score": "S_chi=sum_(j,k) zeta5^k (X_jk-p)",
            "marked_row_estimator": "(325/2) zeta5^(-k_root) [pivotalH4_primal-pivotalH4_matching]",
            "root_schedule": "root_field=replica mod 325",
            "shared_randomness": "same (parent,fiber) Bernoulli field and marked root for both hands",
        },
        "analysis": analysis,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(analysis, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
