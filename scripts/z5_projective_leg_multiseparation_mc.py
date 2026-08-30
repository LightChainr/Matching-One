#!/usr/bin/env python3
"""N325 charged projective-leg insertion smoke for Issue 250."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from math import gcd
from pathlib import Path
from typing import Mapping, Sequence

from integer_period_torus import (
    IntegerHomologyUnionFind,
    matrix_product,
    matrix_vector,
    unimodular_inverse,
)
from norm5_chiral_fixedp_mc import PHASES, P_FIXED, contexts, counter_uniform, mapping_gate
from z5_charged_multiseparation_mc import (
    FIELD_ORDER,
    HANDS,
    SEPARATIONS,
    hermitian_pair,
    observables_for_hand,
    summarize,
    translated_anchors,
)
from z5_charged_threepoint_mc import dft_charges


SCHEMA = "matching-one/z5-projective-leg-multiseparation-response/v1"
SMOKE_CAP = 2_000


class ProjectiveLegIndex:
    """Root membership in black-NN and white-matching rank-one components."""

    def __init__(self, geometry, active: Sequence[bool]) -> None:
        self.geometry = geometry
        self.active = [bool(value) for value in active]
        self.black = self._build(self.active, geometry.primal_edges)
        white = [not value for value in self.active]
        self.white = self._build(white, geometry.matching_edges)

    def _build(self, enabled: Sequence[bool], edges) -> IntegerHomologyUnionFind:
        union = IntegerHomologyUnionFind(self.geometry.n, self.geometry.periods)
        for edge in edges:
            if enabled[edge.i] and enabled[edge.j]:
                union.add_edge(edge.i, edge.j, edge.dx, edge.dy)
        return union

    def _component(self, vertex: int):
        return self.black.component(vertex) if self.active[vertex] else self.white.component(vertex)

    def scalar(self, vertex: int) -> int:
        component = self._component(vertex)
        if component.rank != 1:
            return 0
        return 1 if self.active[vertex] else -1

    def chi4(self, vertex: int) -> float:
        component = self._component(vertex)
        if component.rank != 1:
            return 0.0
        ell = component.basis[0]
        if gcd(abs(ell[0]), abs(ell[1])) != 1:
            raise AssertionError("rank-one component line is not primitive")
        x, y = self.geometry.periods.period_vector(ell)
        z = complex(x, y)
        character = (z**4 / (x * x + y * y) ** 2).real
        return character if self.active[vertex] else -character


def charged_rows(context, field: Sequence[bool], parent_indices: Sequence[int]):
    active = context.active_from_field(field)
    index = ProjectiveLegIndex(context.geometry, active)
    rows = {}
    for parent_index in parent_indices:
        local = [
            index.scalar(context.field_to_vertex[5 * parent_index + fiber])
            for fiber in range(5)
        ]
        rows[parent_index] = dft_charges(local)
    return rows


def _primitive_x_orbit_gate(context) -> dict:
    geometry = context.geometry
    active = [False] * geometry.n
    (a, _), (c, _) = geometry.periods.matrix
    steps = [(1 if a > 0 else -1, 0)] * abs(a)
    steps += [(0, 1 if c > 0 else -1)] * abs(c)
    point = (0, 0)
    orbit = [geometry.vertex(point)]
    for dx, dy in steps:
        point = point[0] + dx, point[1] + dy
        orbit.append(geometry.vertex(point))
    if orbit[-1] != orbit[0]:
        raise AssertionError("period-column Manhattan loop did not close")
    orbit = orbit[:-1]
    for item in set(orbit):
        active[item] = True
    index = ProjectiveLegIndex(geometry, active)
    white_active = [True] * geometry.n
    for item in set(orbit):
        white_active[item] = False
    white_index = ProjectiveLegIndex(geometry, white_active)
    roots = [geometry.vertex((distance, 0)) for distance in (0, 1, 2, 3)]
    component = index.black.component(roots[0])
    white_component = white_index.white.component(roots[0])
    ell = component.basis[0]
    physical = geometry.periods.period_vector(ell)
    chi = complex(*physical) ** 4 / (physical[0] ** 2 + physical[1] ** 2) ** 2
    shear = ((1, 1), (0, 1))
    changed_periods = matrix_product(geometry.periods.matrix, shear)
    changed_ell = matrix_vector(unimodular_inverse(shear), ell)
    changed_physical = matrix_vector(changed_periods, changed_ell)
    return {
        "context": context.name,
        "orbit_vertices": len(set(orbit)),
        "component_rank": component.rank,
        "primitive_line": list(ell),
        "physical_line": list(physical),
        "line_gcd": gcd(abs(ell[0]), abs(ell[1])),
        "black_separated_scalar_values": [index.scalar(root) for root in roots],
        "white_separated_scalar_values": [white_index.scalar(root) for root in roots],
        "white_component_rank": white_component.rank,
        "chi4_abs_residual": abs(abs(chi) - 1.0),
        "chi4_sign_residual": abs(
            chi - (complex(-physical[0], -physical[1]) ** 4 /
                   (physical[0] ** 2 + physical[1] ** 2) ** 2)
        ),
        "basis_change_physical_residual": [
            changed_physical[0] - physical[0], changed_physical[1] - physical[1]
        ],
        "passed": (
            component.rank == 1
            and gcd(abs(ell[0]), abs(ell[1])) == 1
            and [index.scalar(root) for root in roots] == [1, 1, 1, 1]
            and white_component.rank == 1
            and [white_index.scalar(root) for root in roots] == [-1, -1, -1, -1]
            and abs(abs(chi) - 1.0) < 1e-14
            and changed_physical == physical
        ),
    }


def exact_gate() -> dict:
    plus, minus, _ = contexts()
    base = mapping_gate()
    propagation = [_primitive_x_orbit_gate(context) for context in (plus, minus)]
    synthetic = (0, 1, -1, 2, -2)
    shifted = (synthetic[-1],) + synthetic[:-1]
    original_dft = dft_charges(synthetic)
    shifted_dft = dft_charges(shifted)
    phase_residuals = {}
    for charge in range(1, 5):
        phase = complex(*PHASES[(-charge) % 5])
        phase_residuals[str(charge)] = abs(shifted_dft[charge] - phase * original_dft[charge])
    passed = (
        base["passed"]
        and all(row["passed"] for row in propagation)
        and max(phase_residuals.values()) < 1e-14
    )
    if not passed:
        raise AssertionError("projective-leg exact gate failed")
    return {
        "base_cover_gate": base,
        "insertion": (
            "black NN rank-one root membership minus white matching rank-one root membership"
        ),
        "integral_saturation_use": "rank-one component basis is reduced to a primitive line up to sign",
        "constructed_same_component_propagation": propagation,
        "deck_shift": "v_k -> v_(k-1)",
        "deck_character_residuals": phase_residuals,
        "secondary_character": "Re[(P ell)^4/|P ell|^4], sign-invariant and not used for candidate selection",
        "passed": True,
    }


def _run_batch(task: tuple[int, int, int, float, int]) -> dict:
    batch, start, samples, p, seed = task
    plus, minus, _ = contexts()
    sums = {name: 0.0 for name in FIELD_ORDER}
    field_digest = hashlib.sha256()
    translation_digest = hashlib.sha256()
    conjugacy_max = 0.0
    pair_imaginary_max = 0.0
    for replica in range(start, start + samples):
        field = [counter_uniform(seed, replica, site) < p for site in range(325)]
        field_digest.update(bytes(field))
        translation_index, triangles = translated_anchors(seed, replica)
        translation_digest.update(translation_index.to_bytes(2, "little"))
        unique = sorted({index for triangle in triangles.values() for index in triangle})
        for hand, context in (("plus", plus), ("minus", minus)):
            rows = charged_rows(context, field, unique)
            for separation, anchors in triangles.items():
                channels, pairs = observables_for_hand(rows, anchors)
                conjugacy_max = max(
                    conjugacy_max,
                    abs(channels["C244"] - channels["C113"].conjugate()),
                    abs(channels["C334"] - channels["C122"].conjugate()),
                )
                for channel, value in channels.items():
                    sums[f"d{separation}_{channel}_{hand}_re"] += value.real
                    sums[f"d{separation}_{channel}_{hand}_im"] += value.imag
                for name, value in pairs.items():
                    if isinstance(value, complex):
                        pair_imaginary_max = max(pair_imaginary_max, abs(value.imag))
                        value = value.real
                    sums[f"d{separation}_{name}_{hand}"] += value
    return {
        "batch": batch,
        "replica_first": start,
        "samples": samples,
        "field_sha256": field_digest.hexdigest(),
        "translation_sha256": translation_digest.hexdigest(),
        "pair_imaginary_max": pair_imaginary_max,
        "conjugacy_max_abs": conjugacy_max,
        **sums,
    }


def run(
    samples: int,
    batches: int,
    workers: int,
    p: float,
    seed: int,
    replica_offset: int,
    *,
    sample_cap: int = SMOKE_CAP,
):
    if samples <= 0 or batches < 2 or samples % batches:
        raise ValueError("samples must be divisible by batches>=2")
    if samples > sample_cap:
        raise ValueError(f"projective-leg run exceeds authorized cap {sample_cap}")
    per_batch = samples // batches
    tasks = [
        (batch, replica_offset + batch * per_batch, per_batch, p, seed)
        for batch in range(batches)
    ]
    if workers == 1:
        rows = [_run_batch(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(_run_batch, tasks))
    return rows, summarize(rows)


def write_batches(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=2_000)
    parser.add_argument("--batches", type=int, default=40)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--p", type=float, default=P_FIXED)
    parser.add_argument("--seed", type=int, default=25033433720260830)
    parser.add_argument("--replica-offset", type=int, default=0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batches-output", type=Path)
    parser.add_argument("--exact-gate", action="store_true")
    parser.add_argument("--production-manifest", type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.exact_gate:
        args.output.write_text(json.dumps(exact_gate(), indent=2) + "\n")
        return 0
    manifest = None
    sample_cap = SMOKE_CAP
    if args.production_manifest is not None:
        manifest = json.loads(args.production_manifest.read_text())
        if manifest.get("status") != "authorized_fresh_production":
            raise ValueError("production manifest is not authorized")
        expected = manifest["run"]
        observed = {
            "samples": args.samples,
            "batches": args.batches,
            "workers": args.workers,
            "p": args.p,
            "seed": args.seed,
            "replica_offset": args.replica_offset,
            "replica_last_exclusive": args.replica_offset + args.samples,
        }
        for key, value in observed.items():
            if expected.get(key) != value:
                raise ValueError(f"run differs from production manifest for {key}")
        sample_cap = int(expected["samples"])
    rows, analysis = run(
        args.samples,
        args.batches,
        args.workers,
        args.p,
        args.seed,
        args.replica_offset,
        sample_cap=sample_cap,
    )
    batches_path = args.batches_output or args.output.with_suffix(".batches.csv")
    write_batches(batches_path, rows)
    payload = {
        "schema": SCHEMA,
        "status": "frozen_low_sample_operator_propagation_smoke",
        "issues": [250],
        "exact_gate": exact_gate(),
        "run": {
            "samples": args.samples,
            "batches": args.batches,
            "workers": args.workers,
            "p": args.p,
            "seed": args.seed,
            "replica_offset": args.replica_offset,
            "replica_last_exclusive": args.replica_offset + args.samples,
            "batches_output": str(batches_path),
        },
        "observable": {
            "primary_insertion": (
                "black NN rank-one root membership minus white matching rank-one root membership"
            ),
            "deck_projection": "five fiber roots followed by the exact Z5 DFT",
            "separations": list(SEPARATIONS),
            "triangles": "(0,0),(d,0),(0,d)",
            "joint_interface": ["G1", "G2", "V1", "V2", "C113", "C122"],
        },
        "analysis": analysis,
        "production_manifest": (
            str(args.production_manifest) if args.production_manifest is not None else None
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
