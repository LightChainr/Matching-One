#!/usr/bin/env python3
"""Exact physical-toggle obstruction for the P250 two-morphism rectangle.

The existing P215/P225/P334 operations are site overwrite maps.  This
certificate checks their operator algebra and exhausts the first honest square
torus (L=3) with a marked root, occupied deletion site and vacant join site.
It deliberately stops before inventing a state-dependent connector rule.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

try:
    from integer_period_torus import (
        IntegerHomologyUnionFind,
        axis_integer_torus,
    )
except ModuleNotFoundError:
    from scripts.integer_period_torus import (
        IntegerHomologyUnionFind,
        axis_integer_torus,
    )


DEFAULT_OUTPUT = Path("results/exact-p250-physical-toggle-obstruction/latest.json")


def force(mask: int, vertex: int, occupied: bool) -> int:
    """Force one site without retaining history."""

    bit = 1 << vertex
    return mask | bit if occupied else mask & ~bit


def _component_rank(
    geometry,
    active: Sequence[bool],
    root: int,
    *,
    matching: bool,
) -> int:
    union = IntegerHomologyUnionFind(geometry.n, geometry.periods)
    edges = geometry.matching_edges if matching else geometry.primal_edges
    for edge in edges:
        if active[edge.i] and active[edge.j]:
            union.add_edge(edge.i, edge.j, edge.dx, edge.dy)
    return union.component(root).rank


def typed_projective_leg(geometry, mask: int, root: int, *, matching_hand: bool) -> int:
    """Rank-one leg with a declared primary/matching graph hand.

    ``matching_hand=False`` is the production black-NN minus white-matching
    leg.  The other hand exchanges the two graph types.  Complement together
    with this hand exchange is the exact typed involution.
    """

    occupied = [bool(mask & (1 << vertex)) for vertex in range(geometry.n)]
    if occupied[root]:
        return int(
            _component_rank(
                geometry, occupied, root, matching=matching_hand
            )
            == 1
        )
    white = [not value for value in occupied]
    return -int(
        _component_rank(
            geometry, white, root, matching=not matching_hand
        )
        == 1
    )


def overwrite_algebra() -> dict:
    """Check absorption on one support and commutation on disjoint supports."""

    same_support_checks = 0
    disjoint_checks = 0
    for n in range(1, 10):
        for mask in range(1 << n):
            for vertex in range(n):
                d = force(mask, vertex, False)
                j = force(mask, vertex, True)
                # Words are applied left-to-right: DJ means first D, then J.
                dj = force(d, vertex, True)
                jd = force(j, vertex, False)
                if dj != j or jd != d:
                    raise AssertionError("same-site overwrite absorption failed")
                # D+J-DJ-JD is zero pointwise for every response L.
                if sorted((d, j)) != sorted((dj, jd)):
                    raise AssertionError("same-site connected rectangle failed")
                same_support_checks += 1
            for deleted in range(n):
                for joined in range(n):
                    if deleted == joined:
                        continue
                    dj = force(force(mask, deleted, False), joined, True)
                    jd = force(force(mask, joined, True), deleted, False)
                    if dj != jd:
                        raise AssertionError("disjoint site overwrites did not commute")
                    disjoint_checks += 1
    return {
        "word_convention": "DJ means D first, then J",
        "same_support_identity": "DJ=J and JD=D, hence R_plus=0 for every L",
        "disjoint_support_identity": "DJ=JD, hence R_minus=0 for every L",
        "same_support_state_checks": same_support_checks,
        "disjoint_support_state_checks": disjoint_checks,
    }


def typed_involution_gate(geometry) -> dict:
    """Check colour complement, graph-hand exchange and toggle conjugacy."""

    full = (1 << geometry.n) - 1
    leg_checks = 0
    morphism_checks = 0
    for mask in range(1 << geometry.n):
        complement = full ^ mask
        for root in range(geometry.n):
            left = typed_projective_leg(
                geometry, mask, root, matching_hand=False
            )
            right = typed_projective_leg(
                geometry, complement, root, matching_hand=True
            )
            if left != -right:
                raise AssertionError("typed projective-leg involution failed")
            leg_checks += 1
        for vertex in range(geometry.n):
            # C D_v = J_v C and C J_v = D_v C.
            if full ^ force(mask, vertex, False) != force(
                complement, vertex, True
            ):
                raise AssertionError("delete/add complement conjugacy failed")
            if full ^ force(mask, vertex, True) != force(
                complement, vertex, False
            ):
                raise AssertionError("add/delete complement conjugacy failed")
            morphism_checks += 2
    return {
        "involution": "(field, NN-hand) -> (complement field, matching-hand)",
        "leg_identity": "L_NN(omega,r)=-L_matching(C omega,r)",
        "morphism_conjugacy": ["C D_v = J_v C", "C J_v = D_v C"],
        "leg_checks": leg_checks,
        "morphism_checks": morphism_checks,
        "passed": True,
    }


def exhaustive_marked_triples_l3() -> dict:
    """Exhaust all physical D/J rectangles on the first honest square torus."""

    geometry = axis_integer_torus(3)
    n = geometry.n
    values = [
        [
            typed_projective_leg(geometry, mask, root, matching_hand=False)
            for root in range(n)
        ]
        for mask in range(1 << n)
    ]
    r_plus_histogram: dict[int, int] = {}
    r_minus_histogram: dict[int, int] = {}
    witness = None
    checked = 0
    for root in range(n):
        for deleted in range(n):
            for joined in range(n):
                if len({root, deleted, joined}) != 3:
                    continue
                for mask in range(1 << n):
                    if not (mask & (1 << deleted)):
                        continue
                    if mask & (1 << joined):
                        continue
                    d = force(mask, deleted, False)
                    j = force(mask, joined, True)
                    dj = force(d, joined, True)
                    jd = force(j, deleted, False)
                    responses = {
                        "L_0": values[mask][root],
                        "L_D": values[d][root],
                        "L_J": values[j][root],
                        "L_DJ": values[dj][root],
                        "L_JD": values[jd][root],
                    }
                    r_plus = (
                        responses["L_D"]
                        + responses["L_J"]
                        - responses["L_DJ"]
                        - responses["L_JD"]
                    )
                    r_minus = responses["L_DJ"] - responses["L_JD"]
                    r_plus_histogram[r_plus] = r_plus_histogram.get(r_plus, 0) + 1
                    r_minus_histogram[r_minus] = r_minus_histogram.get(r_minus, 0) + 1
                    checked += 1
                    if r_plus and witness is None:
                        witness = {
                            "root": root,
                            "delete_site": deleted,
                            "join_site": joined,
                            "coordinates": {
                                "root": list(geometry.coordinates[root]),
                                "delete_site": list(geometry.coordinates[deleted]),
                                "join_site": list(geometry.coordinates[joined]),
                            },
                            "base_mask": mask,
                            "base_occupied_vertices": [
                                vertex for vertex in range(n) if mask & (1 << vertex)
                            ],
                            "responses": responses,
                            "R_plus": r_plus,
                            "R_minus": r_minus,
                            "interpretation": (
                                "adding the third site closes a rank-one vertical loop; "
                                "deleting the middle occupied site destroys that closure"
                            ),
                        }
    if witness is None or not any(value != 0 for value in r_plus_histogram):
        raise AssertionError("no nonzero symmetric connected rectangle was found")
    if set(r_minus_histogram) != {0}:
        raise AssertionError("fixed distinct site operations unexpectedly retained order")
    return {
        "geometry": "axis L=3 square torus",
        "why_L3": "L=2 has collided nearest-neighbour images; L=3 is the first honest local square quotient",
        "observable": "black-NN rank-one root leg minus white-matching rank-one root leg",
        "marked_triple_contract": "root, occupied delete_site, vacant join_site are distinct",
        "eligible_rectangles_checked": checked,
        "R_plus_histogram": {str(key): value for key, value in sorted(r_plus_histogram.items())},
        "R_minus_histogram": {str(key): value for key, value in sorted(r_minus_histogram.items())},
        "R_plus_nonzero_exists": True,
        "R_minus_nonzero_exists": False,
        "minimal_R_plus_witness": witness,
    }


def runner_semantic_audit() -> dict:
    return {
        "P215": {
            "source_lineages": [
                "3881e88 scripts/homology_rank_birth_insertion.py",
                "dabe28e scripts/rank_birth_parity_channels.py",
            ],
            "operation": "force a fixed root absent, then present, on one frozen field",
            "typed_pair": "black NN insertion and complement-reversed white matching insertion",
            "order_state": "none; one binary counterfactual toggle",
        },
        "P225": {
            "source_lineages": [
                "f9bc73a/89b86f4 scripts/marked_pivotal_h4_reference.py",
                "28042d8 multiradius extension",
            ],
            "operation": "force the fixed pivotal root absent/present and compare wrapping",
            "typed_pair": "black NN root and complementary white matching root",
            "order_state": "none; one binary counterfactual toggle",
        },
        "P334": {
            "source_lineages": [
                "6f54935 scripts/projective_essential_birth_oracle.py",
                "9b8ec23 scripts/marked_birth_path_oracle.py",
            ],
            "operation": "monotone addition of the next fixed permutation site",
            "typed_pair": "P215 complement_pair is evaluated at the next absent site",
            "order_state": (
                "K1/K2 and the first projective line retain path timing, but no physical "
                "delete/join pair or common four-corner rectangle is emitted"
            ),
        },
        "search_result": (
            "No runner contains a state-dependent connector whose support changes after "
            "the first morphism. All available physical candidates reduce to fixed site overwrites."
        ),
    }


def build_result() -> dict:
    geometry = axis_integer_torus(3)
    return {
        "schema": "matching-one/p250-physical-toggle-obstruction/v1",
        "status": "exact_obstruction_no_pilot_interface",
        "runner_semantic_audit": runner_semantic_audit(),
        "overwrite_operator_algebra": overwrite_algebra(),
        "typed_involution_gate": typed_involution_gate(geometry),
        "smallest_honest_torus_certificate": exhaustive_marked_triples_l3(),
        "decision": {
            "pilot_runner_added": False,
            "reason": (
                "Distinct physical D/J site overwrites can have a nonzero symmetric mixed "
                "response R_plus, but R_minus is identically zero because the final fields "
                "coincide. Same-site overwrites make R_plus identically zero by absorption."
            ),
            "minimum_escape_from_no_go": [
                "a connector-selection rule whose microscopic support is a function of the intermediate typed state",
                "the selected support before each branch, so equality is auditable rather than inferred",
                "black-NN and white-matching component IDs/ranks/primitive lines after the first move",
                "L_D,L_J,L_DJ,L_JD on the same base field and marked triple",
            ],
        },
        "claim_boundary": (
            "This rules out the natural fixed-site embedding of formal P333 D/J into the "
            "existing projective-leg runners. It does not rule out a declared state-dependent "
            "cut/connector morphism, and it makes no continuum or path-memory claim."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--stdout", action="store_true")
    args = parser.parse_args()
    payload = build_result()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    if args.stdout:
        print(text, end="")


if __name__ == "__main__":
    main()
