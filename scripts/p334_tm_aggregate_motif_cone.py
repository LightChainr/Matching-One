#!/usr/bin/env python3
"""Four-motif and relative-displacement anatomy of aggregate TM."""

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


def aggregate_face_motifs(marks, line, n: int, k: int, layers):
    """Count ordered two-site faces above the fixed-line lower layer."""

    coexit = 0
    mixed = 0
    synergy = 0
    flat = 0
    for mask in layers[k]:
        exits = []
        internal = []
        for site in range(n):
            if mask >> site & 1:
                continue
            rank, marked_line, _ = marks[mask | (1 << site)]
            if rank == 2:
                exits.append(site)
            elif rank == 1 and marked_line == line:
                internal.append(site)
        coexit += len(exits) * (len(exits) - 1)
        mixed += 2 * len(exits) * len(internal)
        for left in internal:
            for right in internal:
                if left == right:
                    continue
                if marks[mask | (1 << left) | (1 << right)][0] == 2:
                    synergy += 1
                else:
                    flat += 1
    total = len(layers[k]) * (n - k) * (n - k - 1)
    assert coexit + mixed + synergy + flat == total
    return {
        "T": total,
        "D_coexit": coexit,
        "M_mixed": mixed,
        "Y_synergy": synergy,
        "F_flat": flat,
    }


def pair_face_motifs(marks, line, pair, layer):
    left, right = pair
    counts = Counter()
    for mask in layer:
        if mask >> left & 1 or mask >> right & 1:
            continue
        left_rank = marks[mask | (1 << left)][0]
        right_rank = marks[mask | (1 << right)][0]
        double_rank = marks[mask | (1 << left) | (1 << right)][0]
        if left_rank == 2 and right_rank == 2:
            counts["D"] += 1
        elif left_rank == 2 or right_rank == 2:
            counts["M"] += 1
        elif double_rank == 2:
            counts["Y"] += 1
        else:
            counts["F"] += 1
    counts["T"] = sum(counts[key] for key in "DMYF")
    return counts


def curvature_polynomial(motifs):
    total = motifs["T"]
    coexit = motifs.get("D_coexit", motifs.get("D", 0))
    mixed = motifs.get("M_mixed", motifs.get("M", 0))
    synergy = motifs.get("Y_synergy", motifs.get("Y", 0))
    flat = motifs.get("F_flat", motifs.get("F", 0))
    mixed_cover = mixed * mixed
    synergy_cover = 4 * synergy * (total - coexit)
    hard = 4 * coexit * flat
    return mixed_cover, synergy_cover, hard


def quotient_order(geometry, left: int, right: int):
    first = geometry.coordinates[left]
    second = geometry.coordinates[right]
    displacement = second[0] - first[0], second[1] - first[1]
    for order in range(1, geometry.n + 1):
        if geometry.periods.quotient_key(
            (order * displacement[0], order * displacement[1])
        ) == (0, 0):
            return order
    raise AssertionError("a finite quotient displacement had no finite order")


def _descriptor(n, matrix, carrier, line, k, **extra):
    return {
        "N": n,
        "matrix": [list(row) for row in matrix],
        "carrier": carrier,
        "line": list(line),
        "lower_layer": k,
        **extra,
    }


def lower_convex_hull(points):
    """Exact lower hull of Pareto cover-ratio points."""

    def cross(origin, left, right):
        return (left[0] - origin[0]) * (right[1] - origin[1]) - (
            left[1] - origin[1]
        ) * (right[0] - origin[0])

    hull = []
    for point in sorted(points):
        while len(hull) >= 2 and cross(hull[-2], hull[-1], point) <= 0:
            hull.pop()
        hull.append(point)
    return hull


def build_result():
    counters = Counter()
    pareto_points = Counter()
    first_rayleigh_failure = None
    first_synergy_only_failure = None
    first_pair_failure = None
    first_order_failure = None
    first_simple_contrast_failure = None
    first_reflected_margin_mismatch = None
    minimum_corrected_cover = (None, None)

    for n, matrix, geometry in _honest_geometries(12):
        primal_marks = subset_marks(geometry, matching=False)
        matching_marks = subset_marks(geometry, matching=True)
        margin_by_row = {}
        lines = sorted(
            {line for rank, line, _ in primal_marks if rank == 1}
            | {line for rank, line, _ in matching_marks if rank == 1}
        )
        for carrier, marks in (
            ("primal", primal_marks),
            ("matching", matching_marks),
        ):
            for line in lines:
                layers, _ = _local_degrees(marks, line, n)
                for k in range(n):
                    if not layers[k] or not layers[k + 1]:
                        continue
                    counters["line_layer_rows"] += 1
                    motifs = aggregate_face_motifs(marks, line, n, k, layers)
                    demand, typed_supply, synergy_supply, _, _ = (
                        _marked_pair_capacities(marks, line, n, k, layers)
                    )
                    tm_margin = sum(typed_supply.values()) - sum(demand.values())
                    margin_by_row[carrier, line, k] = tm_margin
                    mixed_cover, synergy_cover, hard = curvature_polynomial(motifs)
                    polynomial = mixed_cover + synergy_cover - hard
                    exact_identity = polynomial == 4 * (n - k - 1) * tm_margin
                    counters["four_motif_identity_pass"] += int(exact_identity)
                    assert exact_identity
                    corrected_pass = polynomial >= 0
                    rayleigh_pass = mixed_cover >= hard
                    synergy_only_pass = synergy_cover >= hard
                    counters["corrected_Rayleigh_pass"] += int(corrected_pass)
                    counters["ordinary_Rayleigh_pass"] += int(rayleigh_pass)
                    counters["ordinary_Rayleigh_fail"] += int(not rayleigh_pass)
                    counters["synergy_only_pass"] += int(synergy_only_pass)
                    counters["synergy_only_fail"] += int(not synergy_only_pass)
                    assert corrected_pass
                    descriptor = _descriptor(
                        n,
                        matrix,
                        carrier,
                        line,
                        k,
                        motifs=motifs,
                        mixed_square_cover=mixed_cover,
                        synergy_cover=synergy_cover,
                        hard_coexit_flat=hard,
                        TM_margin=tm_margin,
                    )
                    if not rayleigh_pass and first_rayleigh_failure is None:
                        first_rayleigh_failure = descriptor
                    if not synergy_only_pass and first_synergy_only_failure is None:
                        first_synergy_only_failure = descriptor
                    if hard:
                        point = Fraction(mixed_cover, hard), Fraction(
                            synergy_cover, hard
                        )
                        pareto_points[point] += 1
                        ratio = Fraction(mixed_cover + synergy_cover, hard)
                        if (
                            minimum_corrected_cover[0] is None
                            or ratio < minimum_corrected_cover[0]
                        ):
                            minimum_corrected_cover = ratio, descriptor

                    combined_supply = Counter()
                    for key, capacity in typed_supply.items():
                        combined_supply[_pair_from_supply(key)] += capacity
                    signed = combined_supply.copy()
                    signed.subtract(demand)
                    contrast_minimum = None
                    for left in range(n):
                        for right in range(left + 1, n):
                            diagonal_left = 2 * signed[left, left]
                            diagonal_right = 2 * signed[right, right]
                            contrast = (
                                diagonal_left
                                + diagonal_right
                                - 2 * signed[left, right]
                            )
                            candidate = contrast, left, right
                            if contrast_minimum is None or candidate < contrast_minimum:
                                contrast_minimum = candidate
                    assert contrast_minimum is not None
                    if contrast_minimum[0] < 0:
                        counters["signed_kernel_simple_contrast_negative_rows"] += 1
                        if first_simple_contrast_failure is None:
                            first_simple_contrast_failure = _descriptor(
                                n,
                                matrix,
                                carrier,
                                line,
                                k,
                                contrast_vector=[
                                    contrast_minimum[1],
                                    contrast_minimum[2],
                                ],
                                quadratic_form=contrast_minimum[0],
                            )

                    order_tables = Counter()
                    for left in range(n):
                        for right in range(left + 1, n):
                            pair_motifs = pair_face_motifs(
                                marks, line, (left, right), layers[k]
                            )
                            if not pair_motifs["T"]:
                                continue
                            counters["relative_pair_tables"] += 1
                            pair_polynomial = (
                                sum(curvature_polynomial(pair_motifs)[:2])
                                - curvature_polynomial(pair_motifs)[2]
                            )
                            pair_pass = pair_polynomial >= 0
                            counters["relative_pair_corrected_pass"] += int(pair_pass)
                            counters["relative_pair_corrected_fail"] += int(
                                not pair_pass
                            )
                            if not pair_pass and first_pair_failure is None:
                                first_pair_failure = _descriptor(
                                    n,
                                    matrix,
                                    carrier,
                                    line,
                                    k,
                                    pair=[left, right],
                                    quotient_order=quotient_order(
                                        geometry, left, right
                                    ),
                                    motifs=dict(pair_motifs),
                                    polynomial=pair_polynomial,
                                )
                            order = quotient_order(geometry, left, right)
                            for key in "DMYF":
                                order_tables[order, key] += pair_motifs[key]
                    for order in sorted({order for order, _ in order_tables}):
                        order_motifs = {
                            "T": sum(order_tables[order, key] for key in "DMYF"),
                            **{key: order_tables[order, key] for key in "DMYF"},
                        }
                        order_polynomial = (
                            sum(curvature_polynomial(order_motifs)[:2])
                            - curvature_polynomial(order_motifs)[2]
                        )
                        order_pass = order_polynomial >= 0
                        counters["displacement_order_tables"] += 1
                        counters["displacement_order_corrected_pass"] += int(
                            order_pass
                        )
                        counters["displacement_order_corrected_fail"] += int(
                            not order_pass
                        )
                        if not order_pass and first_order_failure is None:
                            first_order_failure = _descriptor(
                                n,
                                matrix,
                                carrier,
                                line,
                                k,
                                quotient_order=order,
                                motifs=order_motifs,
                                polynomial=order_polynomial,
                            )

        for (carrier, line, k), primal_margin in sorted(margin_by_row.items()):
            if carrier != "primal":
                continue
            reflected_key = "matching", line, n - k - 1
            if reflected_key not in margin_by_row:
                continue
            matching_margin = margin_by_row[reflected_key]
            counters["Alexander_reflected_row_pairs"] += 1
            equality = primal_margin == matching_margin
            counters["Alexander_reflected_margin_equal"] += int(equality)
            counters["Alexander_reflected_margin_unequal"] += int(not equality)
            if not equality and first_reflected_margin_mismatch is None:
                first_reflected_margin_mismatch = _descriptor(
                    n,
                    matrix,
                    "primal",
                    line,
                    k,
                    primal_margin=primal_margin,
                    reflected_matching_layer=n - k - 1,
                    reflected_matching_margin=matching_margin,
                )

    pareto_frontier = []
    pareto_point_values = []
    for point, count in sorted(pareto_points.items()):
        if any(
            other[0] <= point[0]
            and other[1] <= point[1]
            and other != point
            for other in pareto_points
        ):
            continue
        pareto_point_values.append(point)
        pareto_frontier.append(
            {
                "mixed_square_over_hard": str(point[0]),
                "synergy_over_hard": str(point[1]),
                "total_cover_over_hard": str(point[0] + point[1]),
                "rows": count,
            }
        )
    convex_lower_hull = [
        {
            "mixed_square_over_hard": str(point[0]),
            "synergy_over_hard": str(point[1]),
            "total_cover_over_hard": str(point[0] + point[1]),
        }
        for point in lower_convex_hull(pareto_point_values)
    ]

    assert counters["line_layer_rows"] == 984
    assert counters["four_motif_identity_pass"] == 984
    assert counters["corrected_Rayleigh_pass"] == 984
    assert counters["ordinary_Rayleigh_fail"] > 0
    assert counters["synergy_only_fail"] > 0
    assert counters["relative_pair_corrected_fail"] > 0
    assert counters["displacement_order_corrected_fail"] > 0
    minimum_ratio, minimum_row = minimum_corrected_cover
    result = {
        "schema_version": "p334-tm-aggregate-motif-cone-v1",
        "random_two_mark_experiment": {
            "sampling": "choose S uniformly from F_k(ell), then an ordered pair of distinct absent sites",
            "motifs": {
                "D": "both single insertions exit to rank two",
                "M": "exactly one single insertion exits",
                "Y": "both singles stay on ell but the double insertion exits",
                "F": "both singles and the double insertion stay on ell",
            },
            "probability_form": "with p=P(a single mark exits), aggregate TM is P(D) <= p^2 + P(Y)",
        },
        "curvature_corrected_Rayleigh_theorem": {
            "integer_identity": "M^2 + 4 Y (T-D) - 4 D F = 4 (m-1) TM_margin",
            "equivalent_inequality": "4 D F <= M^2 + 4 Y (T-D)",
            "interpretation": "hard coexit-times-flat four-face tokens must inject into mixed-times-mixed tokens or a synergy face paired with any non-coexit face",
            "status": "exact equivalence; this is the unique remaining local motif inequality",
        },
        "bounded_audit": dict(counters),
        "mechanism_independence": {
            "ordinary_Rayleigh_without_synergy": {
                "status": "false",
                "first_failure": first_rayleigh_failure,
            },
            "synergy_without_mixed_square": {
                "status": "false",
                "first_failure": first_synergy_only_failure,
            },
            "minimum_corrected_cover_over_hard": str(minimum_ratio),
            "minimum_row": minimum_row,
            "Pareto_frontier": pareto_frontier,
            "exact_convex_lower_hull": convex_lower_hull,
        },
        "failed_stronger_routes": {
            "relative_displacement_locality": {
                "status": "false",
                "first_failure": first_pair_failure,
            },
            "quotient_order_orbit_pairing": {
                "status": "false",
                "first_failure": first_order_failure,
            },
            "Fourier_sum_of_squares": {
                "status": "false",
                "reason": "the signed supply-minus-demand kernel has an exact negative e_v-e_w quadratic contrast",
                "first_failure": first_simple_contrast_failure,
            },
            "Alexander_primal_matching_identity": {
                "status": "false",
                "reason": "complement reverses an exit face into a birth face rather than the reflected exit motif table",
                "first_margin_mismatch": first_reflected_margin_mismatch,
            },
        },
        "conditional_topology_theorem": {
            "statement": "If every fixed-line two-mark table satisfies the aggregate curvature-corrected Rayleigh inequality after summing all relative displacements, then aggregate TM holds; with translation regularity this also supplies the explicit one-mark Hall injection and every proper cut.",
            "single_missing_motif_inequality": "4 D F <= M^2 + 4 Y (T-D)",
            "proof_target": "construct a topology-respecting injection from coexit x flat face pairs into mixed x mixed or synergy x non-coexit face pairs",
        },
    }
    return json.loads(json.dumps(result))


def render_markdown(result):
    audit = result["bounded_audit"]
    mechanism = result["mechanism_independence"]
    return "\n".join(
        [
            "# Aggregate TM is a curvature-corrected Rayleigh inequality",
            "",
            "Choose a uniform fixed-line state and an ordered pair of distinct absent sites. The two-site face has four possible motifs: coexit `D`, mixed `M`, synergy `Y`, or flat `F`; write `T=D+M+Y+F`.",
            "",
            "Direct integer expansion gives",
            "",
            "`M^2 + 4Y(T-D) - 4DF = 4(m-1) TM_margin`.",
            "",
            "Equivalently, with `p` the one-mark exit probability, TM is `P(D) <= p^2+P(Y)`. The ordinary Rayleigh determinant is corrected by exactly the concave two-site faces where neither single exits but the double insertion does.",
            "",
            "## The two covering mechanisms are both essential",
            "",
            f"Ordinary Rayleigh `4DF<=M^2` fails {audit['ordinary_Rayleigh_fail']} rows. Synergy-only coverage fails {audit['synergy_only_fail']} rows. Their corrected sum passes all {audit['corrected_Rayleigh_pass']} rows and satisfies the integer identity in all {audit['four_motif_identity_pass']} rows. The finite motif cone has {len(mechanism['Pareto_frontier'])} Pareto-minimal cover rays; no single cover term spans it.",
            "",
            "## Why aggregation over every displacement is essential",
            "",
            f"The corrected inequality fails on {audit['relative_pair_corrected_fail']} of {audit['relative_pair_tables']} individual site-pair tables, beginning at N=6. Even grouping by quotient order leaves {audit['displacement_order_corrected_fail']} failures in {audit['displacement_order_tables']} order tables. Thus delta-by-delta and inversion/order pairing are too strong.",
            "",
            "A Fourier sum-of-squares route is also false: an exact two-site contrast is negative in "
            f"{audit['signed_kernel_simple_contrast_negative_rows']} rows. Alexander reflection does not identify the two margins: only {audit['Alexander_reflected_margin_equal']} of {audit['Alexander_reflected_row_pairs']} reflected primal/matching pairs are equal, because an exit face reflects to a birth face.",
            "",
            "## Single remaining motif inequality",
            "",
            "The general topology target is now one explicit four-face injection:",
            "",
            "`coexit x flat -> mixed x mixed  OR  synergy x non-coexit`.",
            "",
            "Proving its aggregate cardinality `4DF<=M^2+4Y(T-D)` after summing all relative displacements is exactly aggregate TM. By the regular-cut theorem, it then supplies the one-mark Hall injection and every proper Hall cut automatically.",
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
