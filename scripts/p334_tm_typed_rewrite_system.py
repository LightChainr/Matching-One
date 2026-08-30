#!/usr/bin/env python3
"""Finite typed rewrite system for the last aggregate-TM critical pair."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path

from p334_dual_hazard_ulc import _local_degrees
from p334_lorentzian_support_gate import _honest_geometries
from p334_tm_aggregate_motif_cone import (
    aggregate_face_motifs,
    curvature_polynomial,
    lower_convex_hull,
)
from projective_essential_birth_oracle import subset_marks


EXIT_PATTERNS = {
    "D": "1222",
    "M_left": "1212",
    "M_right": "1122",
    "Y": "1112",
    "F": "1111",
}


def alexander_dual_pattern(pattern: str):
    """Reverse the square and apply r*=2-r."""

    bottom, left, right, top = (int(value) for value in pattern)
    return "".join(
        str(value)
        for value in (2 - top, 2 - right, 2 - left, 2 - bottom)
    )


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
    motif_totals = Counter()
    point_counts = Counter()
    point_witness = {}
    rescue_signatures = Counter()
    rescue_witnesses = {}
    maximum_rescue_fraction = (Fraction(-1), None)

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
                    motifs = aggregate_face_motifs(marks, line, n, k, layers)
                    for key, value in motifs.items():
                        motif_totals[key] += value
                    assert motifs["M_mixed"] % 2 == 0
                    counters["pattern_1222_D"] += motifs["D_coexit"]
                    counters["pattern_1212_M_left"] += motifs["M_mixed"] // 2
                    counters["pattern_1122_M_right"] += motifs["M_mixed"] // 2
                    counters["pattern_1112_Y"] += motifs["Y_synergy"]
                    counters["pattern_1111_F"] += motifs["F_flat"]

                    mixed_cover, synergy_cover, hard = curvature_polynomial(motifs)
                    mixed_used = min(hard, mixed_cover)
                    residual = hard - mixed_used
                    synergy_used = residual
                    closed = synergy_used <= synergy_cover
                    counters["rewrite_closed_rows"] += int(closed)
                    counters["rewrite_open_rows"] += int(not closed)
                    assert closed
                    if residual:
                        counters["mixed_then_synergy_rows"] += 1
                        signature = tuple(
                            motifs[key]
                            for key in (
                                "T",
                                "D_coexit",
                                "M_mixed",
                                "Y_synergy",
                                "F_flat",
                            )
                        )
                        rescue_signatures[signature] += 1
                        rescue_witnesses.setdefault(
                            signature,
                            _descriptor(
                                n,
                                matrix,
                                carrier,
                                line,
                                k,
                                motifs=motifs,
                                hard_tokens=hard,
                                mixed_tokens_used=mixed_used,
                                synergy_tokens_used=synergy_used,
                                synergy_tokens_available=synergy_cover,
                            ),
                        )
                        fraction = Fraction(synergy_used, synergy_cover)
                        if fraction > maximum_rescue_fraction[0]:
                            maximum_rescue_fraction = (
                                fraction,
                                rescue_witnesses[signature],
                            )
                    else:
                        counters["mixed_only_rows"] += 1
                    counters["hard_tokens"] += hard
                    counters["mixed_tokens_used"] += mixed_used
                    counters["synergy_tokens_used"] += synergy_used
                    counters["unmatched_hard_tokens"] += max(
                        0, residual - synergy_cover
                    )

                    if hard:
                        point = Fraction(mixed_cover, hard), Fraction(
                            synergy_cover, hard
                        )
                        point_counts[point] += 1
                        point_witness.setdefault(
                            point,
                            _descriptor(
                                n,
                                matrix,
                                carrier,
                                line,
                                k,
                                motifs=motifs,
                                mixed_square_over_hard=str(point[0]),
                                synergy_over_hard=str(point[1]),
                            ),
                        )

    pareto_points = [
        point
        for point in sorted(point_counts)
        if not any(
            other[0] <= point[0]
            and other[1] <= point[1]
            and other != point
            for other in point_counts
        )
    ]
    lower_hull = set(lower_convex_hull(pareto_points))
    ray_rows = []
    for index, point in enumerate(pareto_points):
        mixed_ratio, synergy_ratio = point
        if mixed_ratio >= 1:
            rule = "R_M: lexicographically inject every hard token into M x M"
            rescue_fraction = Fraction()
        else:
            rule = "R_M then R_Y: exhaust M x M, then inject the residual into 4 Y x nonD"
            rescue_fraction = (1 - mixed_ratio) / synergy_ratio
        ray_rows.append(
            {
                "ray": f"R{index}",
                "mixed_square_over_hard": str(mixed_ratio),
                "synergy_over_hard": str(synergy_ratio),
                "total_cover_over_hard": str(mixed_ratio + synergy_ratio),
                "bounded_rows": point_counts[point],
                "on_exact_lower_hull": point in lower_hull,
                "canonical_rule": rule,
                "synergy_pool_fraction_used_after_mixed": str(rescue_fraction),
                "minimal_quotient_witness": point_witness[point],
            }
        )

    dual_table = {
        name: {
            "exit_pattern": pattern,
            "Alexander_dual_pattern": alexander_dual_pattern(pattern),
        }
        for name, pattern in EXIT_PATTERNS.items()
    }
    assert counters["line_layer_rows"] == 984
    assert counters["rewrite_closed_rows"] == 984
    assert counters["rewrite_open_rows"] == 0
    assert counters["mixed_only_rows"] == 968
    assert counters["mixed_then_synergy_rows"] == 16
    assert counters["unmatched_hard_tokens"] == 0
    assert len(pareto_points) == 9
    assert len(lower_hull) == 4
    assert len(rescue_signatures) == 4
    rescue_fraction, rescue_row = maximum_rescue_fraction

    result = {
        "schema_version": "p334-tm-typed-rewrite-system-v1",
        "finite_local_type_theorem": {
            "ordered_square_vertices": ["00", "10", "01", "11"],
            "statement": "Above a rank-one fixed-line base, monotonicity permits exactly D, two oriented M types, Y, and F; Alexander complement reverses the square and sends r to 2-r.",
            "patterns": dual_table,
            "observed_ordered_pattern_counts": {
                key: counters[key]
                for key in (
                    "pattern_1222_D",
                    "pattern_1212_M_left",
                    "pattern_1122_M_right",
                    "pattern_1112_Y",
                    "pattern_1111_F",
                )
            },
            "unique_negative_product_type": "D x F",
            "status": "exact exhaustive local rank-pattern classification",
        },
        "canonical_finite_rewrite_system": {
            "hard_tokens": "four orientation replicas of D x F",
            "R_M": "order hard and M x M tokens lexicographically and match as many as possible",
            "R_Y": "match the remaining hard tokens to the lexicographically ordered four replicas of Y x nonD",
            "priority": ["R_M", "R_Y"],
            "termination_measure": "number of unmatched D x F hard tokens",
            "bounded_status": "terminates with zero unmatched tokens on all 984 rows",
            "scientific_boundary": "this is a machine-verifiable aggregate token rewrite; a general topology proof must realize R_M/R_Y by canonical operations on the underlying squares rather than by global lexicographic labels",
        },
        "bounded_audit": dict(counters),
        "bounded_motif_totals": dict(motif_totals),
        "cover_cone_rays": ray_rows,
        "synergy_rescue_classification": {
            "signature_count": len(rescue_signatures),
            "signatures": [
                {
                    "T_D_M_Y_F": list(signature),
                    "rows": count,
                    "minimal_quotient_witness": rescue_witnesses[signature],
                }
                for signature, count in sorted(rescue_signatures.items())
            ],
            "maximum_synergy_pool_fraction_used": str(rescue_fraction),
            "maximum_row": rescue_row,
        },
        "unique_unclosed_general_critical_pair": {
            "type": "D x F",
            "critical_residual": "K=max(0,4DF-M^2)-4Y(T-D)",
            "bounded_value": "K<=0 on every audited row",
            "why_local_patterns_do_not_close_it": "the N=6 opposite-pair counterexample has D and F but no M or Y at that displacement/order, so any realization of R_M/R_Y must exchange mass between relative-displacement classes",
            "missing_topology_rule": "given a coexit square and a flat square, canonically cross-switch their ordered missing marks or invoke the Alexander-dual birth square so that two mixed squares or a synergy-plus-noncoexit pair is produced without collisions",
        },
        "conditional_classification_theorem": {
            "statement": "If every D x F token admits one of the two typed rewrites R_M or R_Y with globally injective images, then aggregate TM holds; translation regularity then yields the one-mark Hall injection and all proper cuts.",
            "status": "exact conditional theorem with one named critical pair",
        },
    }
    return json.loads(json.dumps(result))


def render_markdown(result):
    audit = result["bounded_audit"]
    rescue = result["synergy_rescue_classification"]
    return "\n".join(
        [
            "# A finite typed rewrite system for aggregate TM",
            "",
            "An ordered two-site square above a fixed-line base has only five oriented rank patterns: `D=1222`, `M_left=1212`, `M_right=1122`, `Y=1112`, and `F=1111`. Alexander complement reverses the square and gives `0001`, `0101`, `0011`, `0111`, and `1111` respectively.",
            "",
            "In the curvature-corrected Rayleigh polynomial, the unique negative product type is `D x F`. There are exactly two positive rewrite reservoirs:",
            "",
            "1. `R_M`: two mixed squares (`M x M`);",
            "2. `R_Y`: a synergy square paired with any non-coexit square (`4Y x nonD`).",
            "",
            "Fixing lexicographic token labels and priority `R_M` then `R_Y` gives a terminating, machine-verifiable aggregate rewrite system.",
            "",
            "## Extreme-ray witnesses and rewrite regimes",
            "",
            f"All nine bounded Pareto rays have explicit minimal quotient witnesses; four lie on the exact lower convex hull. Six rays close by `R_M` alone. Three mixed-deficient rays use `R_M` then `R_Y`. Across the full atlas, {audit['mixed_only_rows']} rows are mixed-only and {audit['mixed_then_synergy_rows']} require synergy rescue.",
            f"The 16 rescue rows collapse to {rescue['signature_count']} exact `(T,D,M,Y,F)` signatures. Even in the most demanding signature, only `{rescue['maximum_synergy_pool_fraction_used']}` of the available synergy pool is used after exhausting mixed tokens.",
            "",
            "## The unique unclosed general critical pair",
            "",
            "The bounded rewrite closes all 984 rows with zero unmatched hard tokens. It is not yet a general topology injection: global lexicographic labels do not tell us how to transform the underlying configurations.",
            "",
            "The sole unresolved critical pair is `D x F`. A general rule must cross-switch the two ordered missing-site pairs, or pass through the Alexander-dual birth square, to create `M x M` or `Y x nonD` without image collisions. The exact residual is",
            "",
            "`K=max(0,4DF-M^2)-4Y(T-D)`.",
            "",
            "All audited rows have `K<=0`. The known N=6 displacement counterexample explains why the rule cannot be delta-local: its hard pair has no mixed or synergy cover in the same displacement class.",
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
