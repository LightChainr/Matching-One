#!/usr/bin/env python3
"""Rank Gaussian-isogeny experiments by a frozen information-per-cost proxy.

The calculation is design-only.  It enumerates primitive same-N parent pairs
and primitive Gaussian multipliers, computes exact H4/H8/H12 leverage, and
compares no-fit parent-to-child fingerprints.  The default catalogue uses the
two high-signal parent pairs N=65,85 and returns a concrete four-edge campaign
that separates angular and radial alternatives under one site-update budget.

The statistical proxy is intentionally simple and inspectable: fixed-p paired
DeltaM has a common planning SE at 100M replicas, while wall cost is
proportional to child_N * replicas.  It is a ranking device, not a likelihood
for future observations.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

from gaussian_semigroup_design import Gaussian, fraction_payload


REFERENCE_ALPHA = Fraction(13, 8)
REFERENCE_A4 = 0.7885
SE_PER_100M = 9.2e-5
SITE_UPDATES_PER_SECOND_16CPU = 1.0e8


@dataclass(frozen=True)
class Parent:
    name: str
    first: Gaussian
    second: Gaussian

    @property
    def n(self) -> int:
        if self.first.norm != self.second.norm:
            raise ValueError("parent orientations must have equal norm")
        return self.first.norm


@dataclass(frozen=True)
class Model:
    name: str
    harmonic_m: int
    alpha: Fraction


PARENTS = (
    Parent("N65", Gaussian(8, 1), Gaussian(7, 4)),
    Parent("N85", Gaussian(9, 2), Gaussian(7, 6)),
)

MODELS = (
    Model("H4_alpha_13_8", 1, Fraction(13, 8)),
    Model("H12_alpha_13_8", 3, Fraction(13, 8)),
    Model("H4_alpha_4_3", 1, Fraction(4, 3)),
    Model("H4_alpha_9_8", 1, Fraction(9, 8)),
)


def _delta(parent_or_child: Tuple[Gaussian, Gaussian], harmonic_m: int) -> Fraction:
    first, second = parent_or_child
    return first.cos4m(harmonic_m) - second.cos4m(harmonic_m)


def _canonical_pair(first: Gaussian, second: Gaussian) -> Tuple[Gaussian, Gaussian]:
    return first.canonical_d4(), second.canonical_d4()


def _model_target(
    parent: Parent,
    child: Tuple[Gaussian, Gaussian],
    multiplier_norm: int,
    model: Model,
) -> float:
    """Child target after calibrating every pure model on the parent effect."""
    parent_h4 = _delta((parent.first, parent.second), 1)
    reference_parent_mean = (
        REFERENCE_A4 * float(parent_h4) * parent.n ** (-float(REFERENCE_ALPHA))
    )
    parent_delta = _delta((parent.first, parent.second), model.harmonic_m)
    child_delta = _delta(child, model.harmonic_m)
    if parent_delta == 0:
        raise ValueError("pure-model parent leverage is zero")
    return (
        reference_parent_mean
        * float(child_delta / parent_delta)
        * multiplier_norm ** (-float(model.alpha))
    )


def candidate_payload(
    parent: Parent,
    multiplier: Gaussian,
    *,
    allow_nonprimitive_child: bool = False,
) -> Dict[str, object]:
    child_raw = (
        parent.first.multiply(multiplier),
        parent.second.multiply(multiplier),
    )
    if (
        not allow_nonprimitive_child
        and (child_raw[0].content != 1 or child_raw[1].content != 1)
    ):
        raise ValueError("child pair is not primitive/cyclic")
    if child_raw[0].smith_invariants() != child_raw[1].smith_invariants():
        raise ValueError("child orientations do not have the same translation group")
    child = _canonical_pair(*child_raw)
    if child[0] == child[1]:
        raise ValueError("multiplier collapses the orientation contrast")
    q = multiplier.norm
    targets = {
        model.name: _model_target(parent, child_raw, q, model) for model in MODELS
    }
    exact: Dict[str, object] = {}
    for harmonic_m in (1, 2, 3):
        parent_delta = _delta((parent.first, parent.second), harmonic_m)
        child_delta = _delta(child_raw, harmonic_m)
        ratio = child_delta / parent_delta if parent_delta else None
        exact[f"H{4 * harmonic_m}"] = {
            "parent_delta": fraction_payload(parent_delta),
            "child_delta": fraction_payload(child_delta),
            "angular_ratio": fraction_payload(ratio) if ratio is not None else None,
        }

    angular_gap = abs(targets["H4_alpha_13_8"] - targets["H12_alpha_13_8"])
    radial_values = [
        targets["H4_alpha_13_8"],
        targets["H4_alpha_4_3"],
        targets["H4_alpha_9_8"],
    ]
    radial_gap = min(
        abs(left - right) for left, right in itertools.combinations(radial_values, 2)
    )
    child_n = parent.n * q
    return {
        "id": f"{parent.name}_q{q}_{multiplier.a}_{multiplier.b}",
        "parent": {
            "name": parent.name,
            "N": parent.n,
            "first": parent.first.as_pair(),
            "second": parent.second.as_pair(),
        },
        "multiplier": multiplier.as_pair(),
        "multiplier_norm": q,
        "child": {
            "N": child_n,
            "first_lineage_canonical": child[0].as_pair(),
            "second_lineage_canonical": child[1].as_pair(),
            "smith_invariants": list(child_raw[0].smith_invariants()),
            "cyclic": child_raw[0].cyclic_translation_group,
        },
        "exact_harmonics": exact,
        "targets_from_reference_A4": targets,
        "proxy": {
            "H4_H12_target_gap": angular_gap,
            "minimum_H4_radial_target_gap": radial_gap,
            "balanced_gap_per_sqrt_child_N": min(angular_gap, radial_gap)
            / math.sqrt(child_n),
        },
    }


def primitive_multipliers(max_norm: int) -> Iterable[Gaussian]:
    bound = math.isqrt(max_norm)
    for a in range(1, bound + 1):
        for b in range(-bound, bound + 1):
            value = Gaussian(a, b)
            if value.norm < 2 or value.norm > max_norm or value.content != 1:
                continue
            yield value


def enumerate_candidates(max_norm: int = 65) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    seen = set()
    for parent in PARENTS:
        for multiplier in primitive_multipliers(max_norm):
            try:
                row = candidate_payload(parent, multiplier)
            except ValueError:
                continue
            child = row["child"]
            assert isinstance(child, dict)
            key = (
                parent.name,
                row["multiplier_norm"],
                tuple(child["first_lineage_canonical"]),
                tuple(child["second_lineage_canonical"]),
            )
            if key in seen:
                continue
            seen.add(key)
            rows.append(row)
    rows.sort(
        key=lambda row: float(row["proxy"]["balanced_gap_per_sqrt_child_N"]),  # type: ignore[index]
        reverse=True,
    )
    return rows


def _selected_candidate(
    parent: Parent, multiplier: Tuple[int, int]
) -> Dict[str, object]:
    return candidate_payload(parent, Gaussian(*multiplier))


def selected_campaign() -> List[Dict[str, object]]:
    # q=2 supplies cheap radial leverage; q=5 supplies the first strong H4/H12
    # sign split and the best frozen q2/Jordan S-prime multiplier transfer.
    return [
        _selected_candidate(PARENTS[0], (1, -1)),
        _selected_candidate(PARENTS[1], (1, -1)),
        _selected_candidate(PARENTS[0], (1, 2)),
        _selected_candidate(PARENTS[1], (1, -2)),
    ]


def idealized_norm4_benchmark() -> List[Dict[str, object]]:
    """Return 2x dilation rows excluded from the cyclic production search.

    Multiplication by 2 has norm four and gives the exact T4=T2^2 Gaussian
    semigroup edge.  Its children have Smith invariants (2,2N), so the current
    cyclic C++ engines cannot execute them even though the two orientations at
    each N have the same finite translation group.
    """
    return [
        candidate_payload(parent, Gaussian(2, 0), allow_nonprimitive_child=True)
        for parent in PARENTS
    ]


def _pairwise_chi_square(
    campaign: Sequence[Dict[str, object]], allocations_100m: Sequence[int]
) -> Dict[str, float]:
    result: Dict[str, float] = {}
    for left_index, right_index in itertools.combinations(range(len(MODELS)), 2):
        left, right = MODELS[left_index], MODELS[right_index]
        chi_square = 0.0
        for row, allocation in zip(campaign, allocations_100m):
            targets = row["targets_from_reference_A4"]
            assert isinstance(targets, dict)
            difference = float(targets[left.name]) - float(targets[right.name])
            chi_square += difference * difference * allocation / (SE_PER_100M**2)
        result[f"{left.name}__vs__{right.name}"] = chi_square
    return result


def optimize_allocations(
    campaign: Sequence[Dict[str, object]],
    budget_billion_site_replicas: float = 750.0,
) -> Dict[str, object]:
    """Maximin allocation on a declared 100M-replica grid.

    At least 500M is assigned to each q=5 child, because those two rows carry
    all H4/H12 information.  The upper bound of 2B per edge keeps the search
    finite and matches a practical single-host campaign.
    """
    best = None
    for allocation in itertools.product(range(21), repeat=len(campaign)):
        if allocation[2] < 5 or allocation[3] < 5:
            continue
        cost = math.fsum(
            units * 0.1 * int(row["child"]["N"])  # type: ignore[index]
            for units, row in zip(allocation, campaign)
        )
        if cost > budget_billion_site_replicas:
            continue
        pairwise = _pairwise_chi_square(campaign, allocation)
        objective = min(pairwise.values())
        key = (objective, math.fsum(pairwise.values()), cost)
        if best is None or key > best[0]:
            best = (key, allocation, pairwise)
    if best is None:
        raise ValueError("budget cannot fund the mandatory q=5 minimum")
    (_objective, _sum_chi, cost), allocation, pairwise = best
    total_site_updates = cost * 1e9
    return {
        "allocation_grid_replicas": 100_000_000,
        "allocations": [
            {
                "candidate_id": row["id"],
                "child_N": row["child"]["N"],  # type: ignore[index]
                "replicas": units * 100_000_000,
            }
            for units, row in zip(allocation, campaign)
        ],
        "budget_billion_site_replicas": budget_billion_site_replicas,
        "used_billion_site_replicas": cost,
        "planning_se_at_100m": SE_PER_100M,
        "pairwise_expected_chi_square": pairwise,
        "maximin_expected_chi_square": min(pairwise.values()),
        "maximin_expected_sigma": math.sqrt(min(pairwise.values())),
        "estimated_16cpu_seconds": total_site_updates
        / SITE_UPDATES_PER_SECOND_16CPU,
    }


def design_payload(max_norm: int = 65) -> Dict[str, object]:
    candidates = enumerate_candidates(max_norm)
    campaign = selected_campaign()
    norm4_rows = idealized_norm4_benchmark()
    norm4_plus_norm5 = norm4_rows + campaign[2:]
    return {
        "schema_version": 1,
        "purpose": "exact Gaussian-isogeny maximin design before production",
        "reference": {
            "A4": REFERENCE_A4,
            "alpha_in_N": str(REFERENCE_ALPHA),
            "planning_SE_per_100m_child": SE_PER_100M,
            "site_updates_per_second_16cpu": SITE_UPDATES_PER_SECOND_16CPU,
        },
        "pure_models": [
            {
                "name": model.name,
                "harmonic": 4 * model.harmonic_m,
                "alpha_in_N": str(model.alpha),
            }
            for model in MODELS
        ],
        "enumeration": {
            "max_multiplier_norm": max_norm,
            "parents": [parent.name for parent in PARENTS],
            "candidate_count": len(candidates),
            "top_balanced_candidates": candidates[:12],
        },
        "selected_campaign": campaign,
        "maximin_budget": optimize_allocations(campaign),
        "nonprimitive_norm4_benchmark": {
            "excluded_from_primitive_ranking": True,
            "reason": (
                "Multiplication by 2 produces content-2 children with quotient "
                "groups Z/2 x Z/(2N), while current production uses cyclic Z/N labels."
            ),
            "semigroup_identity": "T_4 = T_2^2 up to the Gaussian unit i",
            "rows": norm4_rows,
            "same_budget_q4_plus_q5": optimize_allocations(norm4_plus_norm5),
            "decision": (
                "Useful as a second-tier radial-curvature and backend-universality "
                "test, but not a replacement for q=2+q=5: under the same fixed-p "
                "SE/cost proxy the best q=4+q=5 maximin design is weaker."
            ),
        },
        "scope_note": (
            "The ranking proxy ignores future parent-child covariance and is "
            "used only to freeze geometries and sample-count scale."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-multiplier-norm", type=int, default=65)
    parser.add_argument("--output")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = design_payload(args.max_multiplier_norm)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        path = Path(args.output)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
