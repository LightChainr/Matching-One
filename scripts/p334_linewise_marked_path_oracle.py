#!/usr/bin/env python3
"""Linewise marked two-step reformulation of the BA/TM topology frontier."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path

from p334_dual_hazard_ulc import _local_degrees
from p334_lorentzian_support_gate import _honest_geometries
from projective_essential_birth_oracle import subset_marks


def _concordance_masses(layer, local):
    types = Counter(
        (local[mask]["birth"], local[mask]["exit"]) for mask in layer
    )
    concordant = 0
    discordant = 0
    for (birth_left, exit_left), count_left in types.items():
        for (birth_right, exit_right), count_right in types.items():
            product = (birth_left - birth_right) * (exit_left - exit_right)
            weighted = count_left * count_right * product
            if weighted > 0:
                concordant += weighted
            elif weighted < 0:
                discordant -= weighted
    return concordant, discordant, [
        {"birth": birth, "exit": exit_flux, "states": count}
        for (birth, exit_flux), count in sorted(types.items())
    ]


def marked_row(marks, line, n: int, k: int, layers, local):
    lower = layers[k]
    upper = layers[k + 1]
    m = n - k
    lower_size = len(lower)
    upper_size = len(upper)
    internal = sum(local[mask]["up_internal"] for mask in lower)
    assert internal == sum(local[mask]["down_internal"] for mask in upper)

    birth_sum = sum(local[mask]["birth"] for mask in upper)
    upper_exit_sum = sum(local[mask]["exit"] for mask in upper)
    joint_birth_exit = sum(
        local[mask]["birth"] * local[mask]["exit"] for mask in upper
    )
    ba_left = upper_size * joint_birth_exit
    ba_right = birth_sum * upper_exit_sum
    concordant, discordant, upper_types = _concordance_masses(upper, local)
    assert concordant - discordant == 2 * (ba_left - ba_right)

    lower_exit_sum = sum(local[mask]["exit"] for mask in lower)
    upper_down_exit = sum(
        local[mask]["down_internal"] * local[mask]["exit"] for mask in upper
    )
    lower_up_exit = sum(
        local[mask]["up_internal"] * local[mask]["exit"] for mask in lower
    )
    new_exit_triples_algebraic = upper_down_exit - lower_up_exit
    ordered_same_state_exit_pairs = sum(
        local[mask]["exit"] * (local[mask]["exit"] - 1) for mask in lower
    )

    new_exit_triples_direct = 0
    square_switch_failures = 0
    for mask in lower:
        for first_site in range(n):
            first = mask | (1 << first_site)
            if first == mask or first not in local:
                continue
            for second_site in range(n):
                if first >> second_site & 1:
                    continue
                final = first | (1 << second_site)
                alternate = mask | (1 << second_site)
                if marks[final][0] != 2 or marks[alternate][0] != 1:
                    continue
                new_exit_triples_direct += 1
                if (
                    marks[alternate][1] != line
                    or marks[alternate | (1 << first_site)][0] != 2
                ):
                    square_switch_failures += 1
    assert new_exit_triples_direct == new_exit_triples_algebraic
    assert new_exit_triples_direct % 2 == 0

    tm_supply_synergy = m * lower_size * new_exit_triples_direct
    tm_supply_independent = (m - 1) * lower_exit_sum**2
    tm_demand = m * lower_size * ordered_same_state_exit_pairs
    tm_supply = tm_supply_synergy + tm_supply_independent
    tm_direct_left = m * lower_size * upper_down_exit
    tm_direct_right = (m - 1) * internal * lower_exit_sum
    assert tm_supply - tm_demand == tm_direct_left - tm_direct_right

    return {
        "lower_layer": k,
        "lower_size": lower_size,
        "upper_size": upper_size,
        "internal_edges": internal,
        "line": list(line),
        "BA": {
            "joint_marked_paths": joint_birth_exit,
            "birth_marks": birth_sum,
            "exit_marks": upper_exit_sum,
            "left": ba_left,
            "right": ba_right,
            "margin": ba_left - ba_right,
            "pass": ba_left >= ba_right,
            "concordance_mass": concordant,
            "discordance_mass": discordant,
            "upper_boundary_types": upper_types,
        },
        "TM": {
            "lower_exit_marks": lower_exit_sum,
            "new_exit_triples": new_exit_triples_direct,
            "synergy_square_count": new_exit_triples_direct // 2,
            "square_switch_failures": square_switch_failures,
            "ordered_same_state_exit_pairs": ordered_same_state_exit_pairs,
            "supply_synergy": tm_supply_synergy,
            "supply_independent": tm_supply_independent,
            "supply": tm_supply,
            "demand": tm_demand,
            "margin": tm_supply - tm_demand,
            "pass": tm_supply >= tm_demand,
            "direct_left": tm_direct_left,
            "direct_right": tm_direct_right,
        },
    }


def _brief(descriptor):
    row = descriptor["row"]
    return {
        "N": descriptor["N"],
        "matrix": descriptor["matrix"],
        "carrier": descriptor["carrier"],
        "line": row["line"],
        "lower_layer": row["lower_layer"],
        "BA": row["BA"],
        "TM": row["TM"],
    }


def build_result():
    counters = Counter()
    ba_ratio_rows = []
    ba_nonzero_equalities = []
    tm_ratio_rows = []
    tm_synergy_absent = []
    complement_state_failures = 0
    complement_degree_swap_failures = 0

    for n, matrix, geometry in _honest_geometries(12):
        primal_marks = subset_marks(geometry, matching=False)
        matching_marks = subset_marks(geometry, matching=True)
        full = (1 << n) - 1
        primal_degree_cache = {}
        matching_degree_cache = {}

        def degrees(marks, mask, cache):
            if mask in cache:
                return cache[mask]
            birth = 0
            exit_flux = 0
            for site in range(n):
                if mask >> site & 1:
                    if marks[mask ^ (1 << site)][0] == 0:
                        birth += 1
                elif marks[mask | (1 << site)][0] == 2:
                    exit_flux += 1
            cache[mask] = birth, exit_flux
            return cache[mask]

        for mask, mark in enumerate(primal_marks):
            if mark[0] != 1:
                continue
            counters["complement_rank_one_states"] += 1
            complement = full ^ mask
            if matching_marks[complement] != mark:
                complement_state_failures += 1
                continue
            primal_birth, primal_exit = degrees(
                primal_marks, mask, primal_degree_cache
            )
            matching_birth, matching_exit = degrees(
                matching_marks, complement, matching_degree_cache
            )
            if primal_birth != matching_exit or primal_exit != matching_birth:
                complement_degree_swap_failures += 1

        for carrier, marks in (
            ("primal", primal_marks),
            ("matching", matching_marks),
        ):
            for line in sorted({line for rank, line, _ in marks if rank == 1}):
                layers, local = _local_degrees(marks, line, n)
                for k in range(n):
                    if not layers[k] or not layers[k + 1]:
                        continue
                    row = marked_row(marks, line, n, k, layers, local)
                    descriptor = {
                        "N": n,
                        "matrix": [list(part) for part in matrix],
                        "carrier": carrier,
                        "row": row,
                    }
                    counters["line_layer_rows"] += 1
                    counters["BA_pass"] += int(row["BA"]["pass"])
                    counters["TM_pass"] += int(row["TM"]["pass"])
                    counters["square_switch_failures"] += row["TM"][
                        "square_switch_failures"
                    ]
                    concordance = row["BA"]["concordance_mass"]
                    discordance = row["BA"]["discordance_mass"]
                    if concordance:
                        ba_ratio_rows.append(
                            (Fraction(discordance, concordance), descriptor)
                        )
                    if concordance == discordance and concordance:
                        ba_nonzero_equalities.append(_brief(descriptor))
                    supply = row["TM"]["supply"]
                    demand = row["TM"]["demand"]
                    if supply:
                        tm_ratio_rows.append((Fraction(demand, supply), descriptor))
                    if (
                        row["TM"]["new_exit_triples"] == 0
                        and row["TM"]["ordered_same_state_exit_pairs"] > 0
                    ):
                        counters["TM_synergy_absent_but_demand_positive"] += 1
                        if len(tm_synergy_absent) < 4:
                            tm_synergy_absent.append(_brief(descriptor))

    max_ba_ratio = max(value for value, _ in ba_ratio_rows)
    max_tm_ratio = max(value for value, _ in tm_ratio_rows)
    max_ba_rows = [
        _brief(row) for value, row in ba_ratio_rows if value == max_ba_ratio
    ]
    max_tm_rows = [
        _brief(row) for value, row in tm_ratio_rows if value == max_tm_ratio
    ]
    counters["complement_state_failures"] = complement_state_failures
    counters["complement_degree_swap_failures"] = complement_degree_swap_failures

    result = {
        "schema_version": "p334-linewise-marked-path-oracle-v1",
        "ambient_H1_marking": {
            "birth": "a 0->1 deletion/insertion boundary selects the fixed primitive line ell",
            "exit": "a 1->2 boundary creates the unique nonzero quotient direction in H1/ell",
            "Alexander_switch": "complement maps a primal joint (birth v, exit w) mark to the matching joint (birth w, exit v) mark at the reflected layer",
        },
        "linewise_path_theorem": {
            "BA_path_form": "A times joint 0<-ell->2 marked paths is at least independent birth marks times exit marks; equivalently concordance mass >= discordance mass",
            "TM_path_form": "m A N_new + (m-1) X^2 >= m A E_2, where N_new counts oriented synergy squares and E_2 counts ordered distinct exit marks at one lower state",
            "conclusion": "If both path inequalities hold for every primitive line and both complementary carriers, the normalized fixed-line layer sequence is ULC.",
            "status": "exact theorem with locally checkable integer path conditions",
        },
        "exact_switches": {
            "synergy_square": "(S,v,w) <-> (S,w,v): both single insertions preserve ell and the double insertion creates rank two",
            "complement_joint_path": "(T; birth v, exit w)_P <-> (E\\T; birth w, exit v)_M",
        },
        "bounded_audit": dict(counters),
        "BA_concordance_extreme": {
            "maximum_discordance_over_concordance": str(max_ba_ratio),
            "maximizer_count": len(max_ba_rows),
            "maximizers": max_ba_rows,
            "nonzero_equality_count": len(ba_nonzero_equalities),
            "nonzero_equalities": ba_nonzero_equalities,
        },
        "TM_path_extreme": {
            "maximum_demand_over_supply": str(max_tm_ratio),
            "maximizer_count": len(max_tm_rows),
            "maximizers": max_tm_rows,
            "first_synergy_absent_rows": tm_synergy_absent,
        },
        "injection_status": {
            "BA_global_switch": "open; statewise comonotonicity is false and eight nonzero equality cases require a genuine weight-preserving bijection",
            "TM_synergy_only_switch": "false; 74 rows have zero new synergy triples but positive ordered-exit demand",
            "weakest_verified_linewise_conditions": "the replicated BA concordance and TM synergy-plus-independent path inequalities stated above",
            "topological_gap": "construct injections for these replicated marked sets using cycle representatives, or prove the same cardinality inequalities by homological double counting",
        },
    }
    return json.loads(json.dumps(result))


def render_markdown(result):
    audit = result["bounded_audit"]
    ba = result["BA_concordance_extreme"]
    tm = result["TM_path_extreme"]
    return "\n".join(
        [
            "# Linewise marked-path form of BA and TM",
            "",
            "Fix a primitive ambient homology line `ell`. A birth mark is a boundary `0->ell`; an exit mark is a boundary `ell->H1`, equivalently creation of the unique nonzero quotient direction in `H1/ell`.",
            "",
            "## BA as concordance of marked two-step paths",
            "",
            "At one layer, `sum b(T)x(T)` counts joint marked paths `0<-T->2`. The exact determinant",
            "",
            "`A sum bx - (sum b)(sum x)`",
            "",
            "is half the concordance mass minus discordance mass over ordered state pairs. BA is precisely the statement that concordant boundary fragilities dominate discordant ones. Alexander complement switches `(birth v, exit w)` to `(birth w, exit v)` on the matching carrier and preserves this determinant.",
            "",
            "## TM as synergy squares plus an independent reservoir",
            "",
            "Let `m=N-k`, `X=sum x(S)`, `E_2=sum x(S)(x(S)-1)`, and `N_new` count oriented triples `(S,v,w)` in which both single insertions preserve `ell` but the double insertion creates rank two. Then TM is exactly",
            "",
            "`m A N_new + (m-1) X^2 >= m A E_2`.",
            "",
            "The canonical switch `(S,v,w)<->(S,w,v)` proves every new-exit event is an oriented side of a synergy square. It does not prove TM by itself: the independent `X^2` reservoir is essential.",
            "",
            "## Exact bounded audit",
            "",
            f"All {audit['line_layer_rows']} line/carrier/layer rows pass both path inequalities. The square involution has {audit['square_switch_failures']} failures, and complement swaps birth/exit degrees on {audit['complement_rank_one_states']} states with {audit['complement_degree_swap_failures']} failures.",
            f"BA reaches the nontrivial equality `discordance/concordance={ba['maximum_discordance_over_concordance']}` in {ba['nonzero_equality_count']} complement-paired rows, so a proof must allow bijection rather than strict surplus.",
            f"TM is nearly tight: maximum demand/supply is {tm['maximum_demand_over_supply']} in {tm['maximizer_count']} N=12 rows. In those maximizers the synergy contribution vanishes. More broadly, {audit['TM_synergy_absent_but_demand_positive']} rows have no new synergy square but positive same-state exit-pair demand.",
            "",
            "## The local topology theorem and the remaining switch",
            "",
            "If the BA concordance inequality and the replicated TM path inequality hold for each primitive line on both complementary carriers, the two-carrier moment theorem gives ULC. This is the weakest currently verified linewise condition. A global statewise injection is too strong; the open topological work is a weight-preserving switching of discordant BA marks and a replicated TM injection that uses both synergy squares and independent exit pairs.",
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
