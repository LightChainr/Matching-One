#!/usr/bin/env python3
"""Collapse the replicated TM Hall family by translation regularity."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path

from p334_dual_hazard_ulc import _local_degrees
from p334_lorentzian_support_gate import _honest_geometries
from p334_tm_replicated_switching_oracle import (
    _marked_pair_capacities,
    _pair_from_supply,
)
from projective_essential_birth_oracle import subset_marks


def weighted_degrees(weights, n: int):
    """Weighted degrees with a loop counted twice."""

    degrees = [0] * n
    for (left, right), capacity in weights.items():
        degrees[left] += capacity
        degrees[right] += capacity
    return degrees


def inside_weight(weights, vertex_bits: int):
    return sum(
        capacity
        for (left, right), capacity in weights.items()
        if vertex_bits >> left & 1 and vertex_bits >> right & 1
    )


def incident_weight(weights, vertex_bits: int):
    return sum(
        capacity
        for (left, right), capacity in weights.items()
        if vertex_bits >> left & 1 or vertex_bits >> right & 1
    )


def cut_weight(weights, vertex_bits: int):
    return sum(
        capacity
        for (left, right), capacity in weights.items()
        if bool(vertex_bits >> left & 1) != bool(vertex_bits >> right & 1)
    )


def gap(weights_demand, weights_supply, vertex_bits: int):
    return incident_weight(weights_supply, vertex_bits) - inside_weight(
        weights_demand, vertex_bits
    )


def _edge_gap(bits: int, supply: int, demand: int):
    left = bool(bits & 1)
    right = bool(bits & 2)
    return supply * int(left or right) - demand * int(left and right)


def uncrossing_truth_table():
    """Verify the one-edge identity from which submodularity follows."""

    rows = []
    failures = 0
    for left_set in range(4):
        for right_set in range(4):
            union = left_set | right_set
            intersection = left_set & right_set
            slack = (
                _edge_gap(left_set, 2, 3)
                + _edge_gap(right_set, 2, 3)
                - _edge_gap(union, 2, 3)
                - _edge_gap(intersection, 2, 3)
            )
            first_side = left_set & ~right_set
            second_side = right_set & ~left_set
            separated = bool(
                (first_side & 1 and second_side & 2)
                or (first_side & 2 and second_side & 1)
            )
            expected = 5 * int(separated)
            failures += int(slack != expected)
            rows.append(
                {
                    "A_membership": left_set,
                    "B_membership": right_set,
                    "slack": slack,
                    "expected": expected,
                }
            )
    return rows, failures


def _descriptor(n, matrix, carrier, line, k, **extra):
    return {
        "N": n,
        "matrix": [list(row) for row in matrix],
        "carrier": carrier,
        "line": list(line),
        "lower_layer": k,
        **extra,
    }


def build_result():
    counters = Counter()
    global_ratio_rows = []
    largest_global_ratio = Fraction(-1)
    largest_proper_ratio = (Fraction(-1), None)
    near_tight_by_size = Counter()
    truth_table, truth_failures = uncrossing_truth_table()

    for n, matrix, geometry in _honest_geometries(12):
        for carrier, marks in (
            ("primal", subset_marks(geometry, matching=False)),
            ("matching", subset_marks(geometry, matching=True)),
        ):
            for line in sorted({line for rank, line, _ in marks if rank == 1}):
                layers, _ = _local_degrees(marks, line, n)
                for k in range(n):
                    if not layers[k] or not layers[k + 1]:
                        continue
                    counters["line_layer_rows"] += 1
                    demand, typed_supply, _, _, _ = _marked_pair_capacities(
                        marks, line, n, k, layers
                    )
                    supply = Counter()
                    for key, capacity in typed_supply.items():
                        supply[_pair_from_supply(key)] += capacity

                    demand_degrees = weighted_degrees(demand, n)
                    supply_degrees = weighted_degrees(supply, n)
                    demand_regular = len(set(demand_degrees)) == 1
                    supply_regular = len(set(supply_degrees)) == 1
                    counters["demand_regular_rows"] += int(demand_regular)
                    counters["supply_regular_rows"] += int(supply_regular)
                    assert demand_regular and supply_regular
                    demand_total = sum(demand.values())
                    supply_total = sum(supply.values())
                    assert sum(demand_degrees) == 2 * demand_total
                    assert sum(supply_degrees) == 2 * supply_total
                    assert demand_total <= supply_total

                    if not demand_total:
                        counters["zero_demand_rows"] += 1
                    else:
                        counters["positive_demand_rows"] += 1
                        global_ratio = Fraction(demand_total, supply_total)
                        row = _descriptor(
                            n,
                            matrix,
                            carrier,
                            line,
                            k,
                            demand=demand_total,
                            supply=supply_total,
                            ratio=str(global_ratio),
                        )
                        if global_ratio > largest_global_ratio:
                            largest_global_ratio = global_ratio
                            global_ratio_rows = [row]
                        elif global_ratio == largest_global_ratio:
                            global_ratio_rows.append(row)

                    full = (1 << n) - 1
                    row_maximum = Fraction(-1)
                    row_maximizers = []
                    for vertex_bits in range(1, 1 << n):
                        counters["nonempty_site_cuts_checked"] += 1
                        size = vertex_bits.bit_count()
                        demand_inside = inside_weight(demand, vertex_bits)
                        supply_incident = incident_weight(supply, vertex_bits)
                        demand_cut = cut_weight(demand, vertex_bits)
                        supply_cut = cut_weight(supply, vertex_bits)
                        decomposition = (
                            Fraction(size, n) * (supply_total - demand_total)
                            + Fraction(demand_cut + supply_cut, 2)
                        )
                        actual_gap = supply_incident - demand_inside
                        counters["regular_cut_decomposition_failures"] += int(
                            decomposition != actual_gap
                        )
                        assert decomposition == actual_gap
                        counters["negative_Hall_cuts"] += int(actual_gap < 0)
                        assert actual_gap >= 0
                        if not demand_inside:
                            continue
                        counters["positive_demand_cuts_checked"] += 1
                        ratio = Fraction(demand_inside, supply_incident)
                        global_ratio = Fraction(demand_total, supply_total)
                        counters["ratio_bound_failures"] += int(
                            ratio > global_ratio
                        )
                        assert ratio <= global_ratio
                        if ratio > row_maximum:
                            row_maximum = ratio
                            row_maximizers = [vertex_bits]
                        elif ratio == row_maximum:
                            row_maximizers.append(vertex_bits)
                        if vertex_bits != full and ratio > largest_proper_ratio[0]:
                            largest_proper_ratio = (
                                ratio,
                                _descriptor(
                                    n,
                                    matrix,
                                    carrier,
                                    line,
                                    k,
                                    vertices=[
                                        site
                                        for site in range(n)
                                        if vertex_bits >> site & 1
                                    ],
                                    demand=demand_inside,
                                    supply=supply_incident,
                                    ratio=str(ratio),
                                    demand_cut=demand_cut,
                                    supply_cut=supply_cut,
                                ),
                            )
                        if ratio >= Fraction(9, 10):
                            near_tight_by_size[
                                ("all_site" if vertex_bits == full else "proper", size)
                            ] += 1

                    if demand_total:
                        unique_full = row_maximizers == [full]
                        counters["unique_all_site_ratio_max_rows"] += int(unique_full)
                        counters["nonunique_or_proper_ratio_max_rows"] += int(
                            not unique_full
                        )
                        assert unique_full

    assert counters["line_layer_rows"] == 984
    assert counters["demand_regular_rows"] == 984
    assert counters["supply_regular_rows"] == 984
    assert counters["regular_cut_decomposition_failures"] == 0
    assert counters["negative_Hall_cuts"] == 0
    assert counters["ratio_bound_failures"] == 0
    assert counters["nonunique_or_proper_ratio_max_rows"] == 0
    assert truth_failures == 0
    proper_ratio, proper_row = largest_proper_ratio

    result = {
        "schema_version": "p334-tm-hall-uncrossing-theorem-v1",
        "translation_regularity": {
            "topological_input": "the finite torus translation group acts transitively on sites and preserves carrier, layer, ambient projective line ell, rank transitions, and the replicated token type",
            "consequence": "the weighted demand graph and combined synergy-plus-reservoir supply multigraph are both regular; a reservoir loop counts twice in weighted degree",
            "status": "exact symmetry consequence",
        },
        "regular_two_mark_cut_theorem": {
            "statement": "For regular nonnegative demand and supply weighted multigraphs on N sites, g(U)=S_inc(U)-D_in(U)=|U|/N*(S_tot-D_tot)+(D_cut(U)+S_cut(U))/2.",
            "corollary": "The all-site aggregate TM inequality implies every induced Hall cut. Moreover D_in(U)/S_inc(U) <= D_tot/S_tot, so all-site is ratio-worst; it is uniquely worst when the combined positive graph connects every proper cut.",
            "proof": [
                "regularity gives 2 D_in(U)+D_cut(U)=2|U|D_tot/N",
                "regularity gives S_inc(U)=|U|S_tot/N+S_cut(U)/2",
                "subtract and use nonnegative cut weights",
            ],
            "status": "exact theorem",
        },
        "uncrossing_lemma": {
            "identity": "g(A)+g(B)-g(A union B)-g(A intersection B)=(D+S)(A\\B,B\\A)>=0",
            "consequence": "g is submodular and the deficit -g is supermodular; maximally dangerous Hall sets uncross by union/intersection",
            "one_edge_truth_table": truth_table,
            "truth_table_failures": truth_failures,
            "status": "exact by summing the two-site identity",
        },
        "bounded_audit": dict(counters),
        "tight_cut_classification": {
            "largest_ratio": str(largest_global_ratio),
            "maximizer_count": len(global_ratio_rows),
            "maximizers": global_ratio_rows,
            "near_tight_ratio_at_least_9_over_10": [
                {"class": key[0], "set_size": key[1], "cuts": count}
                for key, count in sorted(near_tight_by_size.items())
            ],
            "largest_proper_ratio": str(proper_ratio),
            "largest_proper_cut": proper_row,
            "classification": "all positive-demand rows have all-site as the unique ratio maximizer; every >=9/10 cut is all-site, while the strongest proper cut is a complement-of-one-site cut",
        },
        "canonical_cut_reduction": {
            "required_cuts": "one all-site aggregate TM cut per line/carrier/layer row",
            "connected_set_restriction": "true but unnecessary: the unique worst cut is all-site on the connected quotient",
            "line_coset_restriction": "unnecessary: no proper line-coset or other proper set can beat all-site under regularity",
            "bisubmodular_fallback": "not needed; ordinary submodularity plus the regular cut decomposition closes the family",
        },
        "scientific_boundary": {
            "closed": "proper-subset Hall inequalities require no new topological injection once translation regularity and aggregate TM are known",
            "remaining": "prove aggregate TM itself for arbitrary digital Alexander quotients, equivalently the all-site two-carrier moment inequality",
        },
    }
    return json.loads(json.dumps(result))


def render_markdown(result):
    audit = result["bounded_audit"]
    tight = result["tight_cut_classification"]
    return "\n".join(
        [
            "# Translation regularity collapses every TM Hall cut",
            "",
            "Let `D` be the replicated demand graph and `S` the combined synergy-plus-independent supply multigraph. Torus translations preserve the fixed projective line, layer and transition type, so both weighted graphs are regular.",
            "",
            "For every site set `U`, weighted degree counting gives",
            "",
            "`g(U)=S_inc(U)-D_in(U)=|U|/N (S_tot-D_tot) + (D_cut(U)+S_cut(U))/2`.",
            "",
            "Every term on the right is nonnegative once the all-site aggregate TM cut passes. Therefore the entire induced Hall family follows from one canonical cut. Also",
            "",
            "`D_in(U)/S_inc(U) <= D_tot/S_tot`,",
            "",
            "so all-site is the worst ratio cut. Positive independent-reservoir edges cross every proper nonempty cut whenever demand is nonzero, making the maximizer unique.",
            "",
            "## Uncrossing",
            "",
            "The same two-site expansion gives",
            "",
            "`g(A)+g(B)-g(A union B)-g(A intersection B)=(D+S)(A\\B,B\\A)>=0`.",
            "",
            "Thus `g` is submodular and the deficit is supermodular. A bisubmodular fallback is unnecessary: translation regularity supplies a stronger exact decomposition.",
            "",
            "## Exact bounded census",
            "",
            f"All {audit['line_layer_rows']} rows have regular demand and supply degrees. The oracle checks {audit['nonempty_site_cuts_checked']} nonempty site cuts, with zero decomposition, Hall or ratio-bound failures. Among {audit['positive_demand_rows']} nonzero rows, all-site is the unique ratio maximizer every time.",
            f"All ratio-`>=9/10` cuts are all-site. The global maximum is `{tight['largest_ratio']}` in {tight['maximizer_count']} N=12 rows. The strongest proper cut is only `{tight['largest_proper_ratio']}` and is a complement-of-one-site cut.",
            "",
            "## Revised frontier",
            "",
            "The proper-subset Hall family is now a theorem, not an additional conjecture. The remaining topology problem is exactly aggregate TM on arbitrary digital Alexander quotients; connected-set and line-coset reductions add no further burden.",
            "",
        ]
    )


def main():
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
