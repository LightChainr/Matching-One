#!/usr/bin/env python3
"""Exact necessary-condition gates for the projective-line ULC conjecture."""

from __future__ import annotations

import argparse
from collections import deque
from fractions import Fraction
import json
from math import comb
from pathlib import Path
from typing import Iterable, Sequence

from digital_alexander_quotient_frontier import (
    has_four_distinct_face_corners,
    hnf_matrices,
)
from integer_period_torus import integer_torus_geometry
from p334_two_orbit_exact_atlas import _graph_connected
from projective_essential_birth_oracle import subset_marks


def _honest_geometries(maximum_order: int):
    rows = []
    for matrix in hnf_matrices(4, maximum_order):
        geometry = integer_torus_geometry(matrix, name="p334-lorentzian-gate")
        if has_four_distinct_face_corners(geometry) and _graph_connected(geometry):
            rows.append((geometry.n, matrix, geometry))
    return sorted(rows, key=lambda row: (row[0], row[1]))


def _fixed_line_families(geometry):
    marks = subset_marks(geometry, matching=False)
    lines = sorted({line for rank, line, _ in marks if rank == 1})
    for line in lines:
        family = {
            mask
            for mask, (rank, marked_line, _) in enumerate(marks)
            if rank == 1 and marked_line == line
        }
        yield line, family, marks


def homogeneous_exponent(mask: int, n: int) -> tuple[int, ...]:
    """Exponent of x0^(N-|S|) prod_(i in S) x_(i+1)."""

    return (n - mask.bit_count(),) + tuple((mask >> index) & 1 for index in range(n))


def _exponent_to_mask(exponent: Sequence[int], n: int) -> int | None:
    if exponent[0] < 0 or any(value not in (0, 1) for value in exponent[1:]):
        return None
    if exponent[0] != n - sum(exponent[1:]):
        return None
    return sum(value << index for index, value in enumerate(exponent[1:]))


def m_convex_exchange_witness(family: set[int], n: int) -> dict[str, object] | None:
    """Return the first symmetric-exchange failure for a homogeneous support."""

    exponents = {mask: homogeneous_exponent(mask, n) for mask in family}
    for left_mask in sorted(family):
        left = exponents[left_mask]
        for right_mask in sorted(family):
            right = exponents[right_mask]
            for source in range(n + 1):
                if left[source] <= right[source]:
                    continue
                targets = [
                    target
                    for target in range(n + 1)
                    if left[target] < right[target]
                ]
                exchanges = []
                passed = False
                for target in targets:
                    new_left = list(left)
                    new_right = list(right)
                    new_left[source] -= 1
                    new_left[target] += 1
                    new_right[source] += 1
                    new_right[target] -= 1
                    new_left_mask = _exponent_to_mask(new_left, n)
                    new_right_mask = _exponent_to_mask(new_right, n)
                    exchange_pass = (
                        new_left_mask in family and new_right_mask in family
                    )
                    exchanges.append(
                        {
                            "target_coordinate": target,
                            "left_mask": new_left_mask,
                            "right_mask": new_right_mask,
                            "left_in_support": new_left_mask in family,
                            "right_in_support": new_right_mask in family,
                        }
                    )
                    passed |= exchange_pass
                if not passed:
                    return {
                        "left_mask": left_mask,
                        "right_mask": right_mask,
                        "left_sites": [i for i in range(n) if left_mask >> i & 1],
                        "right_sites": [i for i in range(n) if right_mask >> i & 1],
                        "left_exponent": list(left),
                        "right_exponent": list(right),
                        "source_coordinate": source,
                        "candidate_exchanges": exchanges,
                    }
    return None


def first_m_convex_failure(maximum_order: int = 4) -> dict[str, object]:
    checked = 0
    for n, matrix, geometry in _honest_geometries(maximum_order):
        for line, family, marks in _fixed_line_families(geometry):
            checked += 1
            witness = m_convex_exchange_witness(family, n)
            if witness is not None:
                for exchange in witness["candidate_exchanges"]:
                    exchange["left_mark"] = list(marks[exchange["left_mask"]])
                    exchange["right_mark"] = list(marks[exchange["right_mask"]])
                return {
                    "minimality_gate": "all honest-face connected HNFs in increasing N and lexicographic matrix/line order",
                    "fixed_lines_checked": checked,
                    "N": n,
                    "matrix": [list(row) for row in matrix],
                    "coordinates": [list(point) for point in geometry.coordinates],
                    "line": list(line),
                    "family_masks": sorted(family),
                    "witness": witness,
                }
    raise AssertionError("no M-convex support failure found")


def _layer_counts(family: Iterable[int], n: int) -> list[int]:
    counts = [0] * (n + 1)
    for mask in family:
        counts[mask.bit_count()] += 1
    return counts


def _add_edge(graph, source: int, target: int, capacity: int) -> None:
    graph[source].append([target, capacity, len(graph[target])])
    graph[target].append([source, 0, len(graph[source]) - 1])


def _max_flow_with_cut(
    left: Sequence[int], right: Sequence[int], n: int
) -> tuple[int, list[int], list[int]]:
    """Integer max flow for the normalized-matching Hall inequalities."""

    left = list(left)
    right = list(right)
    left_count = len(left)
    right_count = len(right)
    source = 0
    left_offset = 1
    right_offset = left_offset + left_count
    sink = right_offset + right_count
    graph = [[] for _ in range(sink + 1)]
    total = left_count * right_count
    right_index = {mask: index for index, mask in enumerate(right)}
    for index in range(left_count):
        _add_edge(graph, source, left_offset + index, right_count)
    for index in range(right_count):
        _add_edge(graph, right_offset + index, sink, left_count)
    for left_index, mask in enumerate(left):
        for site in range(n):
            upper = mask | (1 << site)
            if upper != mask and upper in right_index:
                _add_edge(
                    graph,
                    left_offset + left_index,
                    right_offset + right_index[upper],
                    total,
                )

    flow = 0
    while True:
        level = [-1] * len(graph)
        level[source] = 0
        pending = deque([source])
        while pending:
            vertex = pending.popleft()
            for target, capacity, _ in graph[vertex]:
                if capacity and level[target] < 0:
                    level[target] = level[vertex] + 1
                    pending.append(target)
        if level[sink] < 0:
            break
        cursor = [0] * len(graph)

        def send(vertex: int, amount: int) -> int:
            if vertex == sink:
                return amount
            while cursor[vertex] < len(graph[vertex]):
                edge = graph[vertex][cursor[vertex]]
                target, capacity, reverse = edge
                if capacity and level[target] == level[vertex] + 1:
                    pushed = send(target, min(amount, capacity))
                    if pushed:
                        edge[1] -= pushed
                        graph[target][reverse][1] += pushed
                        return pushed
                cursor[vertex] += 1
            return 0

        while True:
            pushed = send(source, total)
            if not pushed:
                break
            flow += pushed

    reachable = {source}
    pending = [source]
    while pending:
        vertex = pending.pop()
        for target, capacity, _ in graph[vertex]:
            if capacity and target not in reachable:
                reachable.add(target)
                pending.append(target)
    cut_left = [
        left[index]
        for index in range(left_count)
        if left_offset + index in reachable
    ]
    cut_right = [
        right[index]
        for index in range(right_count)
        if right_offset + index in reachable
    ]
    return flow, cut_left, cut_right


def first_normalized_matching_failure(maximum_order: int = 12) -> dict[str, object]:
    fixed_lines_checked = 0
    layer_pairs_checked = 0
    for n, matrix, geometry in _honest_geometries(maximum_order):
        for line, family, marks in _fixed_line_families(geometry):
            fixed_lines_checked += 1
            layers = [
                sorted(mask for mask in family if mask.bit_count() == k)
                for k in range(n + 1)
            ]
            for k, (lower, upper) in enumerate(zip(layers, layers[1:])):
                if not lower or not upper:
                    continue
                layer_pairs_checked += 1
                flow, cut_lower, cut_upper = _max_flow_with_cut(lower, upper, n)
                total = len(lower) * len(upper)
                if flow == total:
                    continue
                neighbors = sorted(
                    {
                        mask | (1 << site)
                        for mask in cut_lower
                        for site in range(n)
                        if mask | (1 << site) in set(upper)
                    }
                )
                dead_mask = next(
                    (
                        mask
                        for mask in cut_lower
                        if not any(
                            (mask | (1 << site)) in set(upper)
                            for site in range(n)
                        )
                    ),
                    None,
                )
                dead_additions = []
                if dead_mask is not None:
                    for site in range(n):
                        if dead_mask >> site & 1:
                            continue
                        new_mask = dead_mask | (1 << site)
                        dead_additions.append(
                            {
                                "site": site,
                                "coordinate": list(geometry.coordinates[site]),
                                "new_mask": new_mask,
                                "new_mark": list(marks[new_mask]),
                            }
                        )
                return {
                    "minimality_gate": "full exact normalized-matching max-flow test over honest HNFs in increasing order",
                    "fixed_lines_checked": fixed_lines_checked,
                    "layer_pairs_checked": layer_pairs_checked,
                    "N": n,
                    "matrix": [list(row) for row in matrix],
                    "line": list(line),
                    "lower_layer": k,
                    "lower_size": len(lower),
                    "upper_size": len(upper),
                    "required_flow": total,
                    "maximum_flow": flow,
                    "violating_lower_masks": cut_lower,
                    "neighbor_masks": neighbors,
                    "violating_ratio": f"{len(cut_lower)}/{len(lower)} > {len(neighbors)}/{len(upper)}",
                    "dead_end_mask": dead_mask,
                    "dead_end_sites": [
                        index for index in range(n) if dead_mask >> index & 1
                    ],
                    "dead_end_additions": dead_additions,
                    "cut_upper_masks": cut_upper,
                }
    raise AssertionError("no normalized-matching failure found")


def _trim(poly: Sequence[Fraction]) -> list[Fraction]:
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def _derivative(poly: Sequence[Fraction]) -> list[Fraction]:
    return [index * poly[index] for index in range(1, len(poly))]


def _divrem(
    numerator: Sequence[Fraction], denominator: Sequence[Fraction]
) -> tuple[list[Fraction], list[Fraction]]:
    remainder = _trim(numerator)
    denominator = _trim(denominator)
    quotient = [Fraction(0)] * max(1, len(remainder) - len(denominator) + 1)
    while len(remainder) >= len(denominator) and any(remainder):
        degree = len(remainder) - len(denominator)
        scale = remainder[-1] / denominator[-1]
        quotient[degree] = scale
        for index, value in enumerate(denominator):
            remainder[index + degree] -= scale * value
        remainder = _trim(remainder)
    return _trim(quotient), _trim(remainder)


def _gcd(left: Sequence[Fraction], right: Sequence[Fraction]) -> list[Fraction]:
    left = _trim(left)
    right = _trim(right)
    while any(right):
        _, remainder = _divrem(left, right)
        left, right = right, remainder
    scale = left[-1]
    return [value / scale for value in left]


def _sturm(poly: Sequence[Fraction]) -> list[list[Fraction]]:
    squarefree_gcd = _gcd(poly, _derivative(poly))
    squarefree, remainder = _divrem(poly, squarefree_gcd)
    assert not any(remainder)
    sequence = [squarefree, _derivative(squarefree)]
    while any(sequence[-1]):
        _, remainder = _divrem(sequence[-2], sequence[-1])
        if not any(remainder):
            break
        sequence.append([-value for value in remainder])
    return sequence


def _infinity_sign(poly: Sequence[Fraction], positive: bool) -> int:
    sign = 1 if poly[-1] > 0 else -1
    if not positive and (len(poly) - 1) % 2:
        sign *= -1
    return sign


def _variations(signs: Sequence[int]) -> int:
    return sum(left != right for left, right in zip(signs, signs[1:]))


def real_root_audit(poly: Sequence[Fraction]) -> dict[str, object]:
    poly = _trim(poly)
    sequence = _sturm(poly)
    negative = [_infinity_sign(row, False) for row in sequence]
    positive = [_infinity_sign(row, True) for row in sequence]
    squarefree_degree = len(sequence[0]) - 1
    distinct_real_roots = _variations(negative) - _variations(positive)
    return {
        "coefficients_low_to_high": [str(value) for value in poly],
        "degree": len(poly) - 1,
        "squarefree_degree": squarefree_degree,
        "distinct_real_roots": distinct_real_roots,
        "all_roots_real": distinct_real_roots == squarefree_degree,
        "sturm_sequence": [[str(value) for value in row] for row in sequence],
        "signs_at_negative_infinity": negative,
        "signs_at_positive_infinity": positive,
    }


def first_real_root_failure(
    *, normalized: bool, maximum_order: int = 12
) -> dict[str, object]:
    fixed_lines_checked = 0
    for n, matrix, geometry in _honest_geometries(maximum_order):
        for line, family, _ in _fixed_line_families(geometry):
            fixed_lines_checked += 1
            counts = _layer_counts(family, n)
            support = [k for k, value in enumerate(counts) if value]
            low = min(support)
            high = max(support)
            if normalized:
                poly = [Fraction(counts[k], comb(n, k)) for k in range(low, high + 1)]
            else:
                poly = [Fraction(counts[k]) for k in range(low, high + 1)]
            if len(poly) <= 2:
                continue
            audit = real_root_audit(poly)
            if not audit["all_roots_real"]:
                return {
                    "minimality_gate": "honest HNFs in increasing order; zero monomial factor removed",
                    "fixed_lines_checked": fixed_lines_checked,
                    "transform": "normalized_q_generating" if normalized else "raw_count_generating",
                    "N": n,
                    "matrix": [list(row) for row in matrix],
                    "line": list(line),
                    "layer_counts": counts,
                    "removed_zero_root_multiplicity": low,
                    "audit": audit,
                }
    raise AssertionError("no real-rootedness failure found")


def build_result() -> dict[str, object]:
    result = {
        "schema_version": "p334-lorentzian-support-gate-v1",
        "question": "Can standard Lorentzian, normalized-matching, or real-rooted sufficient structures prove the observed fixed-line ULC?",
        "lorentzian_support_gate": first_m_convex_failure(),
        "normalized_q_real_root_gate": first_real_root_failure(normalized=True),
        "raw_count_real_root_gate": first_real_root_failure(normalized=False),
        "normalized_matching_gate": first_normalized_matching_failure(),
        "verdict": {
            "homogenized_multiaffine_lorentzian": "closed_by_minimal_N4_M_convex_exchange_failure",
            "hessian_signature_test": "not_run_support_necessary_condition_already_fails",
            "normalized_q_generating_real_rootedness": "closed_by_exact_N6_quadratic",
            "raw_count_generating_real_rootedness": "closed_by_exact_N8_Sturm_count",
            "fixed_line_normalized_matching_property": "closed_by_exact_N11_layer_cut",
            "remaining_claim": "ULC_through_N12_remains_exact_bounded_evidence_but_requires_a_weaker_direct_rank_sequence_mechanism",
        },
    }
    # Canonicalize tuples from the exact geometry oracle to the checked JSON schema.
    return json.loads(json.dumps(result))


def render_markdown(result: dict[str, object]) -> str:
    m = result["lorentzian_support_gate"]
    q = result["normalized_q_real_root_gate"]
    a = result["raw_count_real_root_gate"]
    nmp = result["normalized_matching_gate"]
    return "\n".join(
        [
            "# Exact gates on the projective-line ULC proof routes",
            "",
            "## 1. The multiaffine Lorentzian route closes at the smallest quotient",
            "",
            f"At N={m['N']}, matrix `{m['matrix']}`, fixed line `{m['line']}`, the family is `{m['family_masks']}`.",
            f"The two support exponents `{m['witness']['left_exponent']}` and `{m['witness']['right_exponent']}` fail symmetric exchange at coordinate {m['witness']['source_coordinate']}; every allowed target leaves at least one exchanged exponent outside the fixed-line support.",
            "Because M-convex support is necessary for a homogeneous Lorentzian polynomial with nonnegative coefficients, the proposed homogenized multiaffine polynomial is not Lorentzian in general. No Hessian test is needed after this support failure.",
            "",
            "## 2. Natural rank-polynomial real-rooted strengthenings also fail",
            "",
            f"The normalized-q generating polynomial first fails at N={q['N']}, `{q['matrix']}`, line `{q['line']}`. After its zero factor is removed its coefficients are `{q['audit']['coefficients_low_to_high']}` and it has {q['audit']['distinct_real_roots']} real roots out of squarefree degree {q['audit']['squarefree_degree']}.",
            f"The raw count generating polynomial first fails at N={a['N']}, `{a['matrix']}`, line `{a['line']}`. Exact Sturm variations leave {a['audit']['distinct_real_roots']} real roots out of squarefree degree {a['audit']['squarefree_degree']}.",
            "Thus the observed ULC cannot be promoted to either of these standard real-rooted statements.",
            "",
            "## 3. Even normalized matching is too strong",
            "",
            f"The first exact normalized-matching failure occurs at N={nmp['N']}, `{nmp['matrix']}`, line `{nmp['line']}`, between layers {nmp['lower_layer']} and {nmp['lower_layer'] + 1}.",
            f"The maximum flow is {nmp['maximum_flow']}/{nmp['required_flow']}; the violating cut is `{nmp['violating_ratio']}`. In particular mask {nmp['dead_end_mask']} is a same-line rank-one state with no same-line one-site extension even though the next layer contains {nmp['upper_size']} states; each missing-site insertion jumps to rank two.",
            "",
            "## Mechanism classification",
            "",
            "- **Exact closure:** multiaffine Lorentzian support, normalized-q real-rootedness, raw-count real-rootedness, and full normalized matching are all false in this family.",
            "- **Still exact bounded evidence:** the normalized layer sequence remains ULC for every line checked through N=12 in the preceding atlas.",
            "- **Revised conjecture:** fixed-line ULC is a rank-sequence phenomenon weaker than the standard support and layer-expansion certificates above.",
            "- **Next proof target:** derive a direct two-step injection or coefficient inequality for `A_k^2 binom(N,k-1) binom(N,k+1) >= A_(k-1) A_(k+1) binom(N,k)^2`; do not seek a global exchange axiom.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = build_result()
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    args.markdown.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()
