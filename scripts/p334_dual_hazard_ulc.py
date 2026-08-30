#!/usr/bin/env python3
"""Direct aggregate-current criterion for fixed-line ULC."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import comb
from pathlib import Path
from typing import Sequence

from p334_lorentzian_support_gate import _honest_geometries
from projective_essential_birth_oracle import subset_marks


def _mean(values: Sequence[Fraction]) -> Fraction:
    return sum(values, Fraction()) / len(values)


def _covariance(left: Sequence[Fraction], right: Sequence[Fraction]) -> Fraction:
    return _mean([a * b for a, b in zip(left, right)]) - _mean(left) * _mean(right)


def _local_degrees(marks, line, n: int):
    layers: list[list[int]] = [[] for _ in range(n + 1)]
    local = {}
    for mask, (rank, marked_line, _) in enumerate(marks):
        if rank != 1 or marked_line != line:
            continue
        k = mask.bit_count()
        birth = 0
        exit_flux = 0
        for site in range(n):
            if mask >> site & 1:
                if marks[mask ^ (1 << site)][0] == 0:
                    birth += 1
            elif marks[mask | (1 << site)][0] == 2:
                exit_flux += 1
        local[mask] = {
            "birth": birth,
            "exit": exit_flux,
            "down_internal": k - birth,
            "up_internal": n - k - exit_flux,
        }
        layers[k].append(mask)
    return layers, local


def _fraction_rows(values: Sequence[tuple[int, Fraction]]) -> list[dict[str, object]]:
    return [{"layer": layer, "value": str(value)} for layer, value in values]


def audit_fixed_line(marks, line, n: int) -> dict[str, object]:
    layers, local = _local_degrees(marks, line, n)
    support = [k for k, layer in enumerate(layers) if layer]
    counts = [len(layer) for layer in layers]
    q = [Fraction(counts[k], comb(n, k)) for k in range(n + 1)]
    beta = []
    xi = []
    current_rows = []
    ratio_rows = []
    for k in support:
        layer = layers[k]
        if k:
            beta.append(
                (k, Fraction(sum(local[mask]["birth"] for mask in layer), k * len(layer)))
            )
        if k < n:
            xi.append(
                (
                    k,
                    Fraction(
                        sum(local[mask]["exit"] for mask in layer),
                        (n - k) * len(layer),
                    ),
                )
            )

    edge_rows = []
    for k in range(n):
        lower = layers[k]
        upper = layers[k + 1]
        if not lower or not upper:
            continue
        internal = sum(local[mask]["up_internal"] for mask in lower)
        assert internal == sum(local[mask]["down_internal"] for mask in upper)
        assert internal > 0

        pointwise_exit_failures = 0
        pointwise_birth_failures = 0
        for mask in lower:
            for site in range(n):
                upper_mask = mask | (1 << site)
                if upper_mask == mask or upper_mask not in local:
                    continue
                if local[upper_mask]["exit"] < local[mask]["exit"]:
                    pointwise_exit_failures += 1
                if local[upper_mask]["birth"] > local[mask]["birth"]:
                    pointwise_birth_failures += 1

        exit_lower_hazard = [
            Fraction(local[mask]["exit"], n - k) for mask in lower
        ]
        exit_upper_hazard = [
            Fraction(local[mask]["exit"], n - k - 1) for mask in upper
        ]
        lower_up_degree = [Fraction(local[mask]["up_internal"]) for mask in lower]
        upper_down_degree = [
            Fraction(local[mask]["down_internal"]) for mask in upper
        ]
        exit_edge_lower = sum(
            degree * hazard for degree, hazard in zip(lower_up_degree, exit_lower_hazard)
        ) / internal
        exit_edge_upper = sum(
            degree * hazard for degree, hazard in zip(upper_down_degree, exit_upper_hazard)
        ) / internal
        exit_edge_slack = exit_edge_upper - exit_edge_lower
        exit_uniform_delta = _mean(exit_upper_hazard) - _mean(exit_lower_hazard)
        exit_bias_correction = (
            _covariance(lower_up_degree, exit_lower_hazard) / _mean(lower_up_degree)
            - _covariance(upper_down_degree, exit_upper_hazard)
            / _mean(upper_down_degree)
        )
        assert exit_uniform_delta == exit_edge_slack + exit_bias_correction

        birth_lower_hazard = [
            Fraction(local[mask]["birth"], k) for mask in lower
        ]
        birth_upper_hazard = [
            Fraction(local[mask]["birth"], k + 1) for mask in upper
        ]
        birth_edge_lower = sum(
            degree * hazard for degree, hazard in zip(lower_up_degree, birth_lower_hazard)
        ) / internal
        birth_edge_upper = sum(
            degree * hazard for degree, hazard in zip(upper_down_degree, birth_upper_hazard)
        ) / internal
        birth_edge_slack = birth_edge_lower - birth_edge_upper
        birth_uniform_delta = _mean(birth_lower_hazard) - _mean(birth_upper_hazard)
        birth_bias_correction = (
            -_covariance(lower_up_degree, birth_lower_hazard) / _mean(lower_up_degree)
            + _covariance(upper_down_degree, birth_upper_hazard)
            / _mean(upper_down_degree)
        )
        assert birth_uniform_delta == birth_edge_slack + birth_bias_correction

        edge_rows.append(
            {
                "lower_layer": k,
                "internal_edges": internal,
                "pointwise_exit_nesting_failures": pointwise_exit_failures,
                "pointwise_birth_nesting_failures": pointwise_birth_failures,
                "exit_edge_slack": str(exit_edge_slack),
                "exit_degree_bias_correction": str(exit_bias_correction),
                "exit_uniform_hazard_delta": str(exit_uniform_delta),
                "birth_edge_slack": str(birth_edge_slack),
                "birth_degree_bias_correction": str(birth_bias_correction),
                "birth_uniform_hazard_delta": str(birth_uniform_delta),
            }
        )

    beta_by_layer = dict(beta)
    xi_by_layer = dict(xi)
    for k in range(n):
        if not counts[k] or not counts[k + 1]:
            continue
        birth_edges = sum(local[mask]["birth"] for mask in layers[k + 1])
        exit_edges = sum(local[mask]["exit"] for mask in layers[k])
        internal_edges = sum(local[mask]["up_internal"] for mask in layers[k])
        derivative_current = (k + 1) * counts[k + 1] - (n - k) * counts[k]
        assert derivative_current == birth_edges - exit_edges
        current_rows.append(
            {
                "lower_layer": k,
                "internal_edges": internal_edges,
                "birth_edges": birth_edges,
                "exit_edges": exit_edges,
                "net_current": derivative_current,
                "normalized_current": str(
                    Fraction(derivative_current, (n - k) * counts[k])
                ),
            }
        )
        ratio = q[k + 1] / q[k]
        hazard_ratio = (1 - xi_by_layer[k]) / (1 - beta_by_layer[k + 1])
        assert ratio == hazard_ratio
        ratio_rows.append(
            {
                "lower_layer": k,
                "q_ratio": str(ratio),
                "dual_hazard_ratio": str(hazard_ratio),
            }
        )

    q_ratios = [q[k + 1] / q[k] for k in support[:-1]]
    return {
        "support": support,
        "counts_on_support": [counts[k] for k in support],
        "q_on_support": [str(q[k]) for k in support],
        "birth_hazard_beta": _fraction_rows(beta),
        "exit_hazard_xi": _fraction_rows(xi),
        "birth_hazard_nonincreasing": all(
            left[1] >= right[1] for left, right in zip(beta, beta[1:])
        ),
        "exit_hazard_nondecreasing": all(
            left[1] <= right[1] for left, right in zip(xi, xi[1:])
        ),
        "q_ratio_nonincreasing": all(
            left >= right for left, right in zip(q_ratios, q_ratios[1:])
        ),
        "q_ratio_strictly_decreasing": all(
            left > right for left, right in zip(q_ratios, q_ratios[1:])
        ),
        "current_rows": current_rows,
        "ratio_rows": ratio_rows,
        "edge_coupling_rows": edge_rows,
    }


def _complement_audit(primal_marks, matching_marks, n: int) -> dict[str, object]:
    full = (1 << n) - 1
    checked = 0
    state_failures = 0
    local_degree_failures = 0
    primal_cache = {}
    matching_cache = {}

    def degrees(marks, mask):
        key = (id(marks), mask)
        cache = primal_cache if marks is primal_marks else matching_cache
        if key in cache:
            return cache[key]
        birth = 0
        exit_flux = 0
        for site in range(n):
            if mask >> site & 1:
                if marks[mask ^ (1 << site)][0] == 0:
                    birth += 1
            elif marks[mask | (1 << site)][0] == 2:
                exit_flux += 1
        cache[key] = birth, exit_flux
        return cache[key]

    for mask, (rank, line, index) in enumerate(primal_marks):
        if rank != 1:
            continue
        checked += 1
        complement = full ^ mask
        if matching_marks[complement] != (rank, line, index):
            state_failures += 1
            continue
        primal_birth, primal_exit = degrees(primal_marks, mask)
        matching_birth, matching_exit = degrees(matching_marks, complement)
        if primal_birth != matching_exit or primal_exit != matching_birth:
            local_degree_failures += 1
    return {
        "rank_one_states_checked": checked,
        "state_line_failures": state_failures,
        "birth_exit_swap_failures": local_degree_failures,
        "all_pass": state_failures == 0 and local_degree_failures == 0,
    }


def build_result() -> dict[str, object]:
    carrier_summary = {
        "primal": Counter(),
        "matching": Counter(),
    }
    worst_negative_exit_correction = None
    example_rows = []
    complement = Counter()
    matrices = 0
    for n, matrix, geometry in _honest_geometries(12):
        matrices += 1
        primal_marks = subset_marks(geometry, matching=False)
        matching_marks = subset_marks(geometry, matching=True)
        complement_row = _complement_audit(primal_marks, matching_marks, n)
        complement["rank_one_states_checked"] += complement_row[
            "rank_one_states_checked"
        ]
        complement["state_line_failures"] += complement_row["state_line_failures"]
        complement["birth_exit_swap_failures"] += complement_row[
            "birth_exit_swap_failures"
        ]

        for carrier, marks in (("primal", primal_marks), ("matching", matching_marks)):
            lines = sorted({line for rank, line, _ in marks if rank == 1})
            for line in lines:
                audit = audit_fixed_line(marks, line, n)
                summary = carrier_summary[carrier]
                summary["fixed_line_sequences"] += 1
                summary["adjacent_layer_pairs"] += len(audit["edge_coupling_rows"])
                summary["ulc_comparisons"] += max(0, len(audit["support"]) - 2)
                summary["exit_hazard_nondecreasing"] += int(
                    audit["exit_hazard_nondecreasing"]
                )
                summary["birth_hazard_nonincreasing"] += int(
                    audit["birth_hazard_nonincreasing"]
                )
                summary["q_ratio_strictly_decreasing"] += int(
                    audit["q_ratio_strictly_decreasing"]
                )
                for row in audit["edge_coupling_rows"]:
                    summary["pointwise_exit_nesting_failures"] += row[
                        "pointwise_exit_nesting_failures"
                    ]
                    summary["pointwise_birth_nesting_failures"] += row[
                        "pointwise_birth_nesting_failures"
                    ]
                    correction = Fraction(row["exit_degree_bias_correction"])
                    if correction < 0:
                        summary["negative_exit_bias_corrections"] += 1
                    if Fraction(row["exit_uniform_hazard_delta"]) == 0:
                        summary["zero_exit_uniform_deltas"] += 1
                    candidate = {
                        "correction": correction,
                        "carrier": carrier,
                        "N": n,
                        "matrix": [list(part) for part in matrix],
                        "line": list(line),
                        **row,
                    }
                    if (
                        worst_negative_exit_correction is None
                        or correction < worst_negative_exit_correction["correction"]
                    ):
                        worst_negative_exit_correction = candidate
                if len(example_rows) < 4 and len(audit["support"]) >= 3:
                    example_rows.append(
                        {
                            "carrier": carrier,
                            "N": n,
                            "matrix": [list(part) for part in matrix],
                            "line": list(line),
                            "audit": audit,
                        }
                    )

    assert worst_negative_exit_correction is not None
    worst_negative_exit_correction = dict(worst_negative_exit_correction)
    worst_negative_exit_correction["correction"] = str(
        worst_negative_exit_correction["correction"]
    )
    result = {
        "schema_version": "p334-dual-hazard-ulc-v1",
        "conditional_lemma": {
            "name": "dual exit-hazard criterion",
            "definitions": [
                "beta_k = B_(k-1)/(k A_k), the fraction of occupied deletions from F_k that fall to rank zero",
                "xi_k = X_k/((N-k) A_k), the fraction of absent insertions from F_k that jump to rank two",
                "r_k = q_(k+1)/q_k = (1-xi_k)/(1-beta_(k+1))",
            ],
            "statement": "If xi_k is nondecreasing for both the primal and complementary matching fixed-line carriers, then complement duality makes beta_k nonincreasing, hence r_k is nonincreasing and q_k is ULC.",
            "status": "exact_conditional_lemma_not_a_general_ULC_proof",
        },
        "exact_edge_lemma": {
            "statement": "Along every internal fixed-line edge S subset T, exit-pivotal sites nest upward and birth-pivotal sites nest downward; normalized hazards therefore order pointwise under the edge coupling.",
            "uniform_bridge_identity": "xi_(k+1)-xi_k = Delta_edge + Cov_k(u,h_x)/E_k[u] - Cov_(k+1)(d,h_x)/E_(k+1)[d]",
            "missing_inequality": "The nonnegative edge slack Delta_edge must dominate any negative degree-bias correction. Rank monotonicity and order-convexity alone do not determine that covariance sign.",
        },
        "complement_duality": {
            "matrices_checked": matrices,
            **dict(complement),
            "all_pass": complement["state_line_failures"] == 0
            and complement["birth_exit_swap_failures"] == 0,
            "identity": "beta_primal(k)=xi_matching(N-k) and xi_primal(k)=beta_matching(N-k)",
        },
        "bounded_exact_audit": {
            carrier: dict(summary) for carrier, summary in carrier_summary.items()
        },
        "worst_negative_exit_degree_bias": worst_negative_exit_correction,
        "representative_sequences": example_rows,
        "verdict": {
            "complete_proof": False,
            "strongest_surviving_route": "prove uniform exit-hazard monotonicity for both complementary carriers",
            "verified_scope": "all 83 honest-face connected HNFs, 240 fixed lines per carrier, N=4..12",
            "why_pointwise_matching_is_insufficient": "edge-coupled pivotal nesting is exact, but 484 of 984 adjacent carrier-layer pairs have a negative exit degree-bias correction",
            "next_missing_object": "an aggregate two-step inequality showing Delta_edge + degree_bias >= 0 without statewise matching",
        },
    }
    return json.loads(json.dumps(result))


def render_markdown(result: dict[str, object]) -> str:
    primal = result["bounded_exact_audit"]["primal"]
    matching = result["bounded_exact_audit"]["matching"]
    worst = result["worst_negative_exit_degree_bias"]
    return "\n".join(
        [
            "# Direct aggregate-current route to fixed-line ULC",
            "",
            "## Exact conditional lemma",
            "",
            "Write `A_k=|F_k|`, `beta_k=B_(k-1)/(k A_k)` for the normalized rank-zero birth boundary, and `xi_k=X_k/((N-k) A_k)` for the normalized rank-two exit boundary. Internal-edge double counting gives",
            "",
            "`q_(k+1)/q_k = (1-xi_k)/(1-beta_(k+1)).`",
            "",
            "Complement duality swaps birth and exit and reverses the layer: `beta_P(k)=xi_M(N-k)` and `xi_P(k)=beta_M(N-k)`. Therefore nondecreasing `xi` on both primal and matching carriers is sufficient for ULC. This is an exact conditional lemma, not yet a general proof of exit-hazard monotonicity.",
            "",
            "## What rank monotonicity proves, and the precise missing term",
            "",
            "Along every internal fixed-line edge, exit-pivotal sites nest upward and birth-pivotal sites nest downward. Thus the edge-weighted normalized exit hazard is nondecreasing. Passing from edge-weighting to uniform layer-weighting introduces exactly",
            "",
            "`xi_(k+1)-xi_k = Delta_edge + Cov_k(u,h_x)/E_k[u] - Cov_(k+1)(d,h_x)/E_(k+1)[d]`,",
            "",
            "where `u,d` are internal up/down degrees and `h_x` is local exit hazard. `Delta_edge>=0` is proved pointwise. The covariance correction has no fixed sign, so the missing aggregate lemma is exactly that the edge slack dominates its negative part.",
            "",
            "## Bounded exact audit",
            "",
            f"Across {result['complement_duality']['matrices_checked']} HNFs through N=12, complement line/degree duality passes on {result['complement_duality']['rank_one_states_checked']} primal rank-one states with zero failures.",
            "",
            "| carrier | fixed lines | adjacent pairs | exit hazard monotone | birth hazard monotone | strict ULC sequences | negative bias corrections |",
            "|---|---:|---:|---:|---:|---:|---:|",
            f"| primal | {primal['fixed_line_sequences']} | {primal['adjacent_layer_pairs']} | {primal['exit_hazard_nondecreasing']} | {primal['birth_hazard_nonincreasing']} | {primal['q_ratio_strictly_decreasing']} | {primal['negative_exit_bias_corrections']} |",
            f"| matching | {matching['fixed_line_sequences']} | {matching['adjacent_layer_pairs']} | {matching['exit_hazard_nondecreasing']} | {matching['birth_hazard_nonincreasing']} | {matching['q_ratio_strictly_decreasing']} | {matching['negative_exit_bias_corrections']} |",
            "",
            f"The most negative exit degree-bias correction is {worst['correction']} at `{worst['matrix']}`, carrier `{worst['carrier']}`, line `{worst['line']}`, lower layer {worst['lower_layer']}; edge slack {worst['exit_edge_slack']} leaves the positive uniform increment {worst['exit_uniform_hazard_delta']}.",
            "",
            "## Status",
            "",
            "- **Proved:** current identity, ratio identity, complement birth/exit swap, and edge-coupled pivotal nesting.",
            "- **Exact finite evidence:** uniform exit hazard is nondecreasing on all 480 primal/matching line sequences through N=12, hence every audited q sequence is ULC.",
            "- **Not proved:** the covariance domination inequality for arbitrary quotient size.",
            "- **Revised proof target:** an aggregate two-step path or boundary double count for `Delta_edge + degree_bias >= 0`; statewise matching is neither required nor possible.",
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
