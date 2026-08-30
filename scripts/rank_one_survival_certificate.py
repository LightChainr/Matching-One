#!/usr/bin/env python3
"""Exact finite survival certificate for rank-one Gaussian-torus states."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
import json
from math import comb
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

try:
    from scripts.integer_period_torus import (
        IntegerTorusGeometry,
        classify_configuration,
        gaussian_integer_torus,
    )
except ModuleNotFoundError:
    from integer_period_torus import (  # type: ignore
        IntegerTorusGeometry,
        classify_configuration,
        gaussian_integer_torus,
    )


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "rank_one_survival_certificate.json"
SCHEMA = "matching-one/rank-one-survival-certificate/v1"


def mask_from_labels(geometry: IntegerTorusGeometry, labels: Iterable[int]) -> int:
    mask = 0
    seen = set()
    for label in labels:
        if not isinstance(label, int) or not 0 <= label < geometry.n:
            raise ValueError("label outside the finite quotient")
        if label in seen:
            raise ValueError("duplicate occupied label")
        seen.add(label)
        mask |= 1 << geometry.vertex((0, label))
    return mask


def active_from_mask(n: int, mask: int) -> tuple[bool, ...]:
    if not isinstance(mask, int) or mask < 0 or mask >= 1 << n:
        raise ValueError("mask outside the finite quotient")
    return tuple(bool(mask & (1 << vertex)) for vertex in range(n))


def rank_and_line(
    geometry: IntegerTorusGeometry, mask: int
) -> tuple[int, Optional[tuple[int, int]]]:
    channels, components = classify_configuration(
        geometry, active_from_mask(geometry.n, mask)
    )
    if channels.max_rank != 1:
        return channels.max_rank, None
    lines = {
        component.basis[0]
        for component in components
        if component.rank == 1
    }
    if len(lines) != 1:
        raise ArithmeticError("rank-one configuration has inconsistent winding lines")
    return 1, next(iter(lines))


class RankCache:
    """Memoized exact rank/line queries on one finite geometry."""

    def __init__(self, geometry: IntegerTorusGeometry) -> None:
        self.geometry = geometry

    @lru_cache(maxsize=None)
    def rank_and_line(self, mask: int) -> tuple[int, Optional[tuple[int, int]]]:
        return rank_and_line(self.geometry, mask)

    def rank(self, mask: int) -> int:
        return self.rank_and_line(mask)[0]


def vacant_vertices(n: int, mask: int) -> tuple[int, ...]:
    active_from_mask(n, mask)
    return tuple(vertex for vertex in range(n) if not mask & (1 << vertex))


def trigger_layers(cache: RankCache, mask: int) -> tuple[tuple[int, ...], tuple[tuple[int, int], ...]]:
    if cache.rank(mask) != 1:
        raise ValueError("trigger layers require a rank-one state")
    vacant = vacant_vertices(cache.geometry.n, mask)
    singleton = tuple(v for v in vacant if cache.rank(mask | (1 << v)) == 2)
    singleton_set = set(singleton)
    pairs = tuple(
        (u, v)
        for u, v in combinations(vacant, 2)
        if u not in singleton_set
        and v not in singleton_set
        and cache.rank(mask | (1 << u) | (1 << v)) == 2
    )
    return singleton, pairs


def subset_survival(cache: RankCache, mask: int, horizon: int) -> Fraction:
    if cache.rank(mask) != 1:
        raise ValueError("survival requires a rank-one state")
    vacant = vacant_vertices(cache.geometry.n, mask)
    if not isinstance(horizon, int) or not 0 <= horizon <= len(vacant):
        raise ValueError("horizon outside the remaining occupation layers")
    surviving = sum(
        cache.rank(mask | sum(1 << vertex for vertex in added)) == 1
        for added in combinations(vacant, horizon)
    )
    return Fraction(surviving, comb(len(vacant), horizon))


def killed_kernel_survival(cache: RankCache, mask: int, horizon: int) -> Fraction:
    if cache.rank(mask) != 1:
        raise ValueError("killed kernel requires a rank-one state")
    vacant = vacant_vertices(cache.geometry.n, mask)
    if not isinstance(horizon, int) or not 0 <= horizon <= len(vacant):
        raise ValueError("horizon outside the remaining occupation layers")

    @lru_cache(maxsize=None)
    def recurse(state: int, remaining: int) -> Fraction:
        if remaining == 0:
            return Fraction(1)
        choices = vacant_vertices(cache.geometry.n, state)
        return sum(
            (
                recurse(state | (1 << vertex), remaining - 1)
                if cache.rank(state | (1 << vertex)) == 1
                else Fraction(0)
            )
            for vertex in choices
        ) / len(choices)

    return recurse(mask, horizon)


def permutation_exit_counts(cache: RankCache, mask: int) -> dict[int, int]:
    if cache.rank(mask) != 1:
        raise ValueError("exit counts require a rank-one state")
    counts: Counter[int] = Counter()
    vacant = vacant_vertices(cache.geometry.n, mask)
    for order in permutations(vacant):
        state = mask
        for step, vertex in enumerate(order, start=1):
            state |= 1 << vertex
            if cache.rank(state) == 2:
                counts[step] += 1
                break
        else:
            counts[len(order) + 1] += 1
    return dict(sorted(counts.items()))


def two_step_identity(cache: RankCache, mask: int) -> Fraction:
    if cache.rank(mask) != 1:
        raise ValueError("two-step identity requires a rank-one state")
    vacant = vacant_vertices(cache.geometry.n, mask)
    if len(vacant) < 2:
        raise ValueError("two-step identity requires at least two vacancies")
    singleton, pairs = trigger_layers(cache, mask)
    safe_pairs = comb(len(vacant) - len(singleton), 2) - len(pairs)
    return Fraction(safe_pairs, comb(len(vacant), 2))


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}" if value.denominator != 1 else str(value.numerator)


def label_for_vertex(geometry: IntegerTorusGeometry, vertex: int) -> int:
    matches = [label for label in range(geometry.n) if geometry.vertex((0, label)) == vertex]
    if len(matches) != 1:
        raise ArithmeticError("declared Gaussian labels are not bijective")
    return matches[0]


def witness_record(cache: RankCache, labels: Sequence[int]) -> dict[str, Any]:
    geometry = cache.geometry
    mask = mask_from_labels(geometry, labels)
    rank, line = cache.rank_and_line(mask)
    singleton, pairs = trigger_layers(cache, mask)
    subset = [subset_survival(cache, mask, horizon) for horizon in range(4)]
    kernel = [killed_kernel_survival(cache, mask, horizon) for horizon in range(4)]
    if subset != kernel:
        raise ArithmeticError("subset and killed-kernel survival disagree")
    return {
        "occupied_labels": list(labels),
        "k": len(labels),
        "rank": rank,
        "primitive_line": list(line) if line is not None else None,
        "singleton_trigger_labels": [label_for_vertex(geometry, v) for v in singleton],
        "minimal_trigger_pair_labels": [
            [label_for_vertex(geometry, u), label_for_vertex(geometry, v)]
            for u, v in pairs
        ],
        "b1": len(singleton),
        "b2": len(pairs),
        "survival_s0_to_s3": [fraction_text(value) for value in subset],
        "future_permutation_exit_counts": {
            str(step): count for step, count in permutation_exit_counts(cache, mask).items()
        },
    }


def build_artifact() -> dict[str, Any]:
    geometry = gaussian_integer_torus(3, 1)
    cache = RankCache(geometry)
    rank_one_states = 0
    pair_identity_states = 0
    horizon_checks = 0
    for mask in range(1 << geometry.n):
        if cache.rank(mask) != 1:
            continue
        rank_one_states += 1
        q = len(vacant_vertices(geometry.n, mask))
        for horizon in range(q + 1):
            if subset_survival(cache, mask, horizon) != killed_kernel_survival(cache, mask, horizon):
                raise ArithmeticError("subset and killed-kernel survival disagree")
            horizon_checks += 1
        if q >= 2:
            if two_step_identity(cache, mask) != subset_survival(cache, mask, 2):
                raise ArithmeticError("minimal-trigger pair identity failed")
            pair_identity_states += 1

    first = witness_record(cache, (0, 1, 2, 3, 4))
    second = witness_record(cache, (0, 1, 2, 3, 5))
    return {
        "schema": SCHEMA,
        "issue": 403,
        "status": "exact_rank_one_survival_certificate",
        "geometry": {
            "gaussian_period": "3+i",
            "period_matrix": [[3, -1], [1, 3]],
            "site_count": geometry.n,
            "label_rule": "label j is the coset represented by (0,j)",
        },
        "exhaustive_checks": {
            "rank_one_states": rank_one_states,
            "two_step_identity_states": pair_identity_states,
            "subset_vs_killed_kernel_horizon_checks": horizon_checks,
        },
        "witnesses": {"A": first, "B": second},
        "claim_boundary": {
            "included": "finite N=10 rank-one survival identities and counterexample",
            "excluded": "production pilot changes, sampling, scaled horizons, predictive-state learning, continuum or field identification",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    regenerated = build_artifact()
    if artifact != regenerated:
        raise ValueError("rank-one survival artifact does not exactly reproduce")
    witnesses = regenerated["witnesses"]
    return {
        "schema": SCHEMA,
        "status": "valid_exact_rank_one_survival_certificate",
        "rank_one_states": regenerated["exhaustive_checks"]["rank_one_states"],
        "A_s2": witnesses["A"]["survival_s0_to_s3"][2],
        "B_s2": witnesses["B"]["survival_s0_to_s3"][2],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        artifact = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_artifact(artifact), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
