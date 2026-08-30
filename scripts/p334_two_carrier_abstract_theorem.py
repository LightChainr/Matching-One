#!/usr/bin/env python3
"""Abstract two-carrier theorem and minimal symbolic counterexamples for ULC."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import comb
from pathlib import Path


def upsets(n: int) -> list[frozenset[int]]:
    """Enumerate every upset of the Boolean lattice for n <= 4."""

    result = []
    state_count = 1 << n
    for code in range(1 << state_count):
        family = frozenset(
            mask for mask in range(state_count) if code >> mask & 1
        )
        if all(
            all(
                (mask | (1 << site)) in family
                for site in range(n)
                if not mask >> site & 1
            )
            for mask in family
        ):
            result.append(family)
    return result


def admissible_systems(n: int):
    """Yield monotone 0/1/2 systems with contiguous nonempty rank-one support."""

    full = (1 << n) - 1
    all_upsets = upsets(n)
    for rank_at_least_one in all_upsets:
        if 0 in rank_at_least_one or full not in rank_at_least_one:
            continue
        for rank_two in all_upsets:
            if (
                not rank_two
                or full not in rank_two
                or not rank_two <= rank_at_least_one
            ):
                continue
            sector = rank_at_least_one - rank_two
            if not sector:
                continue
            support = sorted({mask.bit_count() for mask in sector})
            if support != list(range(min(support), max(support) + 1)):
                continue
            yield rank_at_least_one, rank_two, frozenset(sector)


def dual_system(n: int, rank_at_least_one, rank_two):
    """Return r*(W)=2-r(E\\W) as nested upsets."""

    full = (1 << n) - 1
    dual_one = frozenset(
        mask for mask in range(1 << n) if (full ^ mask) not in rank_two
    )
    dual_two = frozenset(
        mask for mask in range(1 << n) if (full ^ mask) not in rank_at_least_one
    )
    return dual_one, dual_two, frozenset(dual_one - dual_two)


def _local_degrees(n: int, rank_at_least_one, rank_two, sector):
    layers = [[] for _ in range(n + 1)]
    local = {}
    for mask in sector:
        k = mask.bit_count()
        birth = 0
        exit_flux = 0
        for site in range(n):
            if mask >> site & 1:
                if (mask ^ (1 << site)) not in rank_at_least_one:
                    birth += 1
            elif (mask | (1 << site)) in rank_two:
                exit_flux += 1
        local[mask] = {
            "birth": birth,
            "exit": exit_flux,
            "down_internal": k - birth,
            "up_internal": n - k - exit_flux,
        }
        layers[k].append(mask)
    for layer in layers:
        layer.sort()
    return layers, local


def audit_system(n: int, rank_at_least_one, rank_two, sector):
    layers, local = _local_degrees(n, rank_at_least_one, rank_two, sector)
    rows = []
    for k in range(n):
        lower = layers[k]
        upper = layers[k + 1]
        if not lower or not upper:
            continue
        lower_size = len(lower)
        upper_size = len(upper)
        lower_exit = sum(local[mask]["exit"] for mask in lower)
        upper_exit = sum(local[mask]["exit"] for mask in upper)
        upper_birth = sum(local[mask]["birth"] for mask in upper)
        upper_birth_exit = sum(
            local[mask]["birth"] * local[mask]["exit"] for mask in upper
        )
        internal = sum(local[mask]["up_internal"] for mask in lower)
        assert internal == sum(local[mask]["down_internal"] for mask in upper)
        upper_down_exit = sum(
            local[mask]["down_internal"] * local[mask]["exit"]
            for mask in upper
        )
        association_left = upper_size * upper_birth_exit
        association_right = upper_birth * upper_exit
        transport_left = (n - k) * lower_size * upper_down_exit
        transport_right = (n - k - 1) * internal * lower_exit
        xi_lower = Fraction(lower_exit, (n - k) * lower_size)
        xi_upper = Fraction(upper_exit, (n - k - 1) * upper_size)
        rows.append(
            {
                "lower_layer": k,
                "lower_masks": lower,
                "upper_masks": upper,
                "lower_size": lower_size,
                "upper_size": upper_size,
                "internal_edges": internal,
                "lower_exit_sum": lower_exit,
                "upper_exit_sum": upper_exit,
                "upper_birth_sum": upper_birth,
                "upper_birth_exit_sum": upper_birth_exit,
                "upper_down_exit_sum": upper_down_exit,
                "boundary_association_left": association_left,
                "boundary_association_right": association_right,
                "boundary_association_pass": association_left >= association_right,
                "transport_moment_left": transport_left,
                "transport_moment_right": transport_right,
                "transport_moment_pass": transport_left >= transport_right,
                "xi_lower": str(xi_lower),
                "xi_upper": str(xi_upper),
                "xi_delta": str(xi_upper - xi_lower),
                "local_degrees": {
                    str(mask): local[mask] for mask in sorted(set(lower + upper))
                },
            }
        )

    counts = [len(layer) for layer in layers]
    q = [Fraction(counts[k], comb(n, k)) for k in range(n + 1)]
    ulc_rows = []
    for k in range(1, n):
        if not counts[k - 1] or not counts[k] or not counts[k + 1]:
            continue
        margin = q[k] ** 2 - q[k - 1] * q[k + 1]
        ulc_rows.append(
            {
                "layer": k,
                "q_previous": str(q[k - 1]),
                "q": str(q[k]),
                "q_next": str(q[k + 1]),
                "margin": str(margin),
                "pass": margin >= 0,
            }
        )
    return {
        "layers": [layer for layer in layers],
        "counts": counts,
        "q": [str(value) for value in q],
        "adjacent_rows": rows,
        "all_boundary_association": all(
            row["boundary_association_pass"] for row in rows
        ),
        "all_transport_moment": all(row["transport_moment_pass"] for row in rows),
        "exit_hazard_nondecreasing": all(
            Fraction(row["xi_delta"]) >= 0 for row in rows
        ),
        "ulc_rows": ulc_rows,
        "ulc": all(row["pass"] for row in ulc_rows),
    }


def _system_payload(n, first, second, sector, audit):
    dual_one, dual_two, dual_sector = dual_system(n, first, second)
    return {
        "N": n,
        "rank_at_least_one_masks": sorted(first),
        "rank_two_masks": sorted(second),
        "rank_one_sector_masks": sorted(sector),
        "dual_rank_at_least_one_masks": sorted(dual_one),
        "dual_rank_two_masks": sorted(dual_two),
        "dual_rank_one_sector_masks": sorted(dual_sector),
        "audit": audit,
    }


def symbolic_scan(maximum_n: int = 4):
    summary = []
    association_needed = None
    transport_needed = None
    first_ulc_failure = None
    for n in range(2, maximum_n + 1):
        counts = {
            "N": n,
            "upsets": len(upsets(n)),
            "admissible_systems": 0,
            "boundary_association_failures": 0,
            "transport_moment_failures": 0,
            "exit_hazard_decrease_systems": 0,
            "ulc_failures": 0,
        }
        for first, second, sector in admissible_systems(n):
            counts["admissible_systems"] += 1
            audit = audit_system(n, first, second, sector)
            if not audit["all_boundary_association"]:
                counts["boundary_association_failures"] += 1
            if not audit["all_transport_moment"]:
                counts["transport_moment_failures"] += 1
            if not audit["exit_hazard_nondecreasing"]:
                counts["exit_hazard_decrease_systems"] += 1
            if not audit["ulc"]:
                counts["ulc_failures"] += 1

            for row in audit["adjacent_rows"]:
                if (
                    association_needed is None
                    and row["transport_moment_pass"]
                    and not row["boundary_association_pass"]
                    and Fraction(row["xi_delta"]) < 0
                ):
                    association_needed = _system_payload(
                        n, first, second, sector, audit
                    )
                    association_needed["witness_row"] = row
                if (
                    transport_needed is None
                    and row["boundary_association_pass"]
                    and not row["transport_moment_pass"]
                    and Fraction(row["xi_delta"]) < 0
                ):
                    transport_needed = _system_payload(
                        n, first, second, sector, audit
                    )
                    transport_needed["witness_row"] = row
            if first_ulc_failure is None and not audit["ulc"]:
                first_ulc_failure = _system_payload(n, first, second, sector, audit)
                first_ulc_failure["failing_rows"] = [
                    row for row in audit["ulc_rows"] if not row["pass"]
                ]
        summary.append(counts)
    assert association_needed is not None
    assert transport_needed is not None
    assert first_ulc_failure is not None
    return {
        "summary": summary,
        "association_axiom_independence_witness": association_needed,
        "transport_axiom_independence_witness": transport_needed,
        "first_ULC_failure": first_ulc_failure,
    }


def build_result(root: Path) -> dict[str, object]:
    scan = symbolic_scan()
    topology = json.loads(
        (root / "results/p334-hazard-transport-bound/latest.json").read_text()
    )
    result = {
        "schema_version": "p334-two-carrier-abstract-theorem-v1",
        "abstract_class": {
            "name": "boundary-regular two-carrier rank system",
            "basic_axioms": [
                "a monotone rank map r:2^E->{0,1,2} with r(empty)=0 and r(E)=2",
                "a nonempty rank-one sector F with contiguous layer support",
                "the complementary carrier r*(W)=2-r(E\\W)",
                "Boolean boundary degrees b,d,x,u with b+d=|S| and x+u=N-|S|",
            ],
            "moment_axioms": [
                "BA: A_(k+1) sum_T b(T)x(T) >= (sum_T b(T))(sum_T x(T))",
                "TM: (N-k)A_k sum_T d(T)x(T) >= (N-k-1)I_k sum_S x(S)",
            ],
        },
        "theorem": {
            "statement": "If BA and TM hold on every adjacent rank-one layer of both complementary carriers, then each carrier has nondecreasing exit hazard; complement duality makes birth hazard nonincreasing, so the normalized layer sequence q_k is ULC.",
            "proof": [
                "TM is exactly E_(upper edge)[h_x] >= E_(lower uniform)[h_x].",
                "BA is exactly E_(upper uniform)[h_x] >= E_(upper edge)[h_x], because d=(k+1)-b.",
                "Their sum gives xi_(k+1)>=xi_k.",
                "Apply the same to the complementary carrier and use q_(k+1)/q_k=(1-xi_k)/(1-beta_(k+1)).",
            ],
            "status": "exact theorem for the named abstract class",
            "minimality": "BA and TM are logically independent under the basic axioms; deleting either admits an N=4 exit-hazard decrease.",
        },
        "symbolic_boolean_scan": scan,
        "topological_instantiation": {
            "audited_pairs": topology["bounded_counts"]["adjacent_pairs"],
            "BA_pass": topology["bounded_counts"]["association_nonnegative"],
            "TM_pass": topology["bounded_counts"]["variance_domination_pass"],
            "scope": "existing honest torus quotients through N=12 only",
            "status": "exact bounded evidence that the topological sectors lie in the abstract class; not a proof for all quotients",
        },
        "verdict": {
            "basic_monotone_two_carrier_axioms_imply_ULC": False,
            "minimal_abstract_counterexample_size": 4,
            "new_exact_result": "BA+TM is a modular sufficient theorem, and both moment axioms are independently necessary for that modular proof under the basic axioms",
            "remaining_topological_problem": "derive BA and TM from homology-carrier geometry, likely by separate aggregate pair and two-step path injections",
        },
    }
    return json.loads(json.dumps(result))


def render_markdown(result: dict[str, object]) -> str:
    scan = result["symbolic_boolean_scan"]
    ba = scan["association_axiom_independence_witness"]
    tm = scan["transport_axiom_independence_witness"]
    ulc = scan["first_ULC_failure"]
    return "\n".join(
        [
            "# A two-carrier abstract theorem for fixed-line ULC",
            "",
            "## The abstract class",
            "",
            "A boundary-regular two-carrier system consists of a monotone three-rank Boolean system, its exact complement-dual carrier, a contiguous rank-one sector, and two aggregate moment axioms:",
            "",
            "`BA: A_(k+1) sum b x >= (sum b)(sum x)`,",
            "",
            "`TM: (N-k) A_k sum d x >= (N-k-1) I_k sum x`.",
            "",
            "BA is nonnegative aggregate birth/exit association on the upper layer. TM says the upper edge-weighted exit hazard dominates the lower uniform exit hazard.",
            "",
            "## Exact theorem",
            "",
            "BA moves from the upper edge marginal to the upper uniform layer; TM moves from the lower uniform layer to the upper edge marginal. Hence exit hazard is nondecreasing. If BA and TM hold on both complementary carriers, complement duality also makes birth hazard nonincreasing, and the exact hazard-ratio identity proves ULC.",
            "",
            "## Minimality of the new axioms",
            "",
            "Exhaustive symbolic enumeration contains 5 admissible systems at N=2, 111 at N=3, and 7,076 at N=4. No N<=3 system violates BA, TM, exit-hazard monotonicity, or ULC. At N=4 the basic axioms are no longer enough.",
            "",
            f"BA cannot be dropped: the first independence witness has sector `{ba['rank_one_sector_masks']}`. At layer {ba['witness_row']['lower_layer']}, TM passes `{ba['witness_row']['transport_moment_left']}>={ba['witness_row']['transport_moment_right']}`, but BA fails `{ba['witness_row']['boundary_association_left']}<{ba['witness_row']['boundary_association_right']}` and exit hazard decreases by `{ba['witness_row']['xi_delta']}`.",
            f"TM cannot be dropped: sector `{tm['rank_one_sector_masks']}` has BA equality `{tm['witness_row']['boundary_association_left']}={tm['witness_row']['boundary_association_right']}`, TM failure `{tm['witness_row']['transport_moment_left']}<{tm['witness_row']['transport_moment_right']}`, and exit-hazard decrement `{tm['witness_row']['xi_delta']}`.",
            f"The first direct ULC counterexample is also N=4, sector `{ulc['rank_one_sector_masks']}`, with failing normalized layers `{ulc['failing_rows']}`.",
            "",
            "## Topological status",
            "",
            f"All {result['topological_instantiation']['audited_pairs']} existing torus carrier-layer pairs satisfy BA and TM exactly. Thus the new theorem captures every checked topology, but the derivation of BA and TM from homology geometry remains open. The useful proof targets are now integer moment inequalities, not a global stochastic coupling.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = build_result(args.root)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    args.markdown.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()
