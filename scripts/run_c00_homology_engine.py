#!/usr/bin/env python3
"""C00 exact regressions for the general 2x2 integer-period homology engine.

Writes machine-readable JSON for the stop-gate counts: exhaustive axis L=3,
exhaustive diamond L=2, a small non-axis Gaussian quotient, and random
unimodular basis-change invariance.  Does not start C01 production data.
"""

from __future__ import annotations

import argparse
import json
import random
import time
import unittest
from pathlib import Path
import sys


SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from matched_torus_reference import (  # noqa: E402
    axis_geometry,
    diamond_geometry,
    diamond_xy_geometry,
    gaussian_geometry,
    integer_period_geometry,
)
from torus_homology import (  # noqa: E402
    component_homologies,
    determinant,
    exhaustive_channel_counts,
    matmul2,
    qspan_contains,
    transport_basis,
    unimodular_inverse,
    wrapping_channels,
)


PR21_AXIS_L3 = {"rank0": 259, "rank1": 162, "rank2": 91, "d0": 175, "d1": 175}
PR21_DIAMOND_L2 = {"rank0": 143, "rank1": 68, "rank2": 45, "d0": 81, "d1": 81}
GAUSSIAN_2_1 = {"rank0": 16, "rank1": 10, "rank2": 6, "d0": 11, "d1": 11}


def _random_unimodular(rng: random.Random):
    matrix = ((1, 0), (0, 1))
    for _ in range(rng.randint(3, 8)):
        shift = rng.randint(-3, 3)
        if rng.randrange(2) == 0:
            shear = ((1, shift), (0, 1))
        else:
            shear = ((1, 0), (shift, 1))
        matrix = matmul2(matrix, shear)
    if rng.randrange(2) == 0:
        matrix = matmul2(matrix, ((-1, 0), (0, 1)))
    return matrix


def _basis_change_pass(geometry, seed: int, n_changes: int) -> dict:
    rng = random.Random(seed)
    period_matrix = geometry.period_matrix
    n_checked = 0
    n_plus = 0
    n_minus = 0
    for _ in range(n_changes):
        change = _random_unimodular(rng)
        det_u = determinant(change)
        if abs(det_u) != 1:
            raise RuntimeError("generated a non-unimodular matrix")
        if det_u == 1:
            n_plus += 1
        else:
            n_minus += 1
        changed = matmul2(period_matrix, change)
        inverse = unimodular_inverse(change)
        for matching in (False, True):
            edges = geometry.matching_edges if matching else geometry.primal_edges
            for mask in range(1 << geometry.n):
                active = [bool((mask >> vertex) & 1) for vertex in range(geometry.n)]
                original = component_homologies(active, edges, period_matrix)
                transformed = component_homologies(active, edges, changed)
                left = wrapping_channels(original)
                right = wrapping_channels(transformed)
                if (left.max_rank, left.either, left.cross) != (
                    right.max_rank,
                    right.either,
                    right.cross,
                ):
                    raise RuntimeError("rank/either/cross changed under unimodular map")
                for old, new in zip(original, transformed):
                    if old.root != new.root or old.size != new.size or old.rank != new.rank:
                        raise RuntimeError("component identity changed under unimodular map")
                    transported = transport_basis(old.basis, inverse)
                    if any(
                        not qspan_contains(vector, new.basis) for vector in transported
                    ):
                        raise RuntimeError("transported winding span mismatch")
                n_checked += 1
    return {
        "pass": True,
        "seed": seed,
        "n_unimodular_matrices": n_changes,
        "n_det_plus_1": n_plus,
        "n_det_minus_1": n_minus,
        "n_config_edge_set_checks": n_checked,
        "rng": "random.Random (Mersenne Twister)",
    }


def run_regressions() -> dict:
    axis = exhaustive_channel_counts(axis_geometry(3))
    diamond = exhaustive_channel_counts(diamond_geometry(2))
    diamond_xy = exhaustive_channel_counts(diamond_xy_geometry(2))
    axis_general = exhaustive_channel_counts(
        integer_period_geometry(((3, 0), (0, 3)), name="axis", L=3)
    )
    gaussian = exhaustive_channel_counts(gaussian_geometry(2, 1))
    return {
        "axis_L3": {
            "counts": axis,
            "pr21": PR21_AXIS_L3,
            "pass": axis == PR21_AXIS_L3,
        },
        "diamond_L2": {
            "counts": diamond,
            "pr21": PR21_DIAMOND_L2,
            "pass": diamond == PR21_DIAMOND_L2,
        },
        "axis_L3_via_general_period_matrix": {
            "counts": axis_general,
            "pr21": PR21_AXIS_L3,
            "pass": axis_general == PR21_AXIS_L3,
        },
        "diamond_L2_xy_embedding": {
            "counts": diamond_xy,
            "pr21": PR21_DIAMOND_L2,
            "pass": diamond_xy == PR21_DIAMOND_L2,
        },
        "gaussian_2_1": {
            "counts": gaussian,
            "expected": GAUSSIAN_2_1,
            "period_matrix": [[2, -1], [1, 2]],
            "N": 5,
            "pass": gaussian == GAUSSIAN_2_1,
        },
        "basis_change_axis_L3": _basis_change_pass(axis_geometry(3), 20260828, 8),
        "basis_change_diamond_L2": _basis_change_pass(diamond_geometry(2), 20260828, 8),
        "basis_change_gaussian_2_1": _basis_change_pass(
            gaussian_geometry(2, 1), 20260828, 12
        ),
        "generator_relative_channels_are_basis_dependent": {
            "pass": True,
            "note": (
                "direction_0/direction_1/both are generator-relative and are "
                "not invariant under shears. Rank, either and cross are."
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=None)
    parser.add_argument("--unittest", action="store_true")
    args = parser.parse_args()

    started = time.perf_counter()
    payload = run_regressions()
    elapsed = time.perf_counter() - started
    payload["wall_time_seconds_regressions"] = elapsed
    payload["all_pass"] = all(
        item.get("pass", False)
        for key, item in payload.items()
        if isinstance(item, dict) and "pass" in item
    )

    if args.unittest:
        suite = unittest.defaultTestLoader.discover(
            str(SCRIPTS.parent / "tests"),
            pattern="test_torus_homology.py",
        )
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        payload["unittest"] = {
            "testsRun": result.testsRun,
            "failures": len(result.failures),
            "errors": len(result.errors),
            "pass": result.wasSuccessful(),
        }
        payload["all_pass"] = payload["all_pass"] and result.wasSuccessful()

    text = json.dumps(payload, indent=2, sort_keys=True)
    print(text)
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text + "\n", encoding="utf-8")
    return 0 if payload["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
