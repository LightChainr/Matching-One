#!/usr/bin/env python3
"""Exact N=9 reservoir obstruction and minimal two-mark repair for P334."""

from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import combinations, product
import json
from pathlib import Path
from typing import Any, Iterable

from p334_tm_coarse_reservoir_hall import (
    capacitated_hall,
    coarse_row_audit,
    coarse_source_key,
)
from p334_tm_configuration_cross_switch import (
    _oriented_source,
    target_token,
    translate_mask,
    translation_permutations,
)
from p334_tm_corrected_reservoir_scan import rows_for_order
from p334_tm_translation_orbit_hall import (
    descriptor,
    inverse_to_origin,
    normalize_source,
    normalize_target,
    transverse_reservoir_targets,
)


SCHEMA = "p334-n9-reservoir-obstruction-v1"
FAILURE_MOTIFS = {"D": 216, "M": 432, "Y": 0, "F": 72}


def strict_descriptor(row) -> dict[str, Any]:
    n, matrix, _geometry, carrier, _marks, _line, lower_layer, faces = row
    payload = descriptor(n, matrix, carrier, (0, 0), lower_layer, faces)
    payload.pop("line")
    return payload


def selected_rows() -> list[tuple[int, tuple]]:
    rows = rows_for_order(9)
    gate = {
        "N": 9,
        "matrix": [[3, 0], [0, 3]],
        "Smith_invariants": [3, 3],
        "carrier": "matching",
        "lower_layer": 4,
        "motifs": FAILURE_MOTIFS,
    }
    selected = [
        (index, row)
        for index, row in enumerate(rows)
        if strict_descriptor(row) == gate
    ]
    if [index for index, _row in selected] != [1, 3]:
        raise AssertionError("the frozen N9 strict descriptor must select rows 1 and 3")
    return selected


def _coarse_classes(row):
    n, _matrix, geometry, _carrier, _marks, _line, _layer, faces = row
    permutations = translation_permutations(geometry)
    inverses = inverse_to_origin(permutations, n)
    raw_sources = [
        (replica, coexit, flat)
        for replica in range(4)
        for coexit in faces["D"]
        for flat in faces["F"]
    ]
    orbit_sources = sorted({normalize_source(source, inverses) for source in raw_sources})
    twins: dict[tuple, list[tuple]] = defaultdict(list)
    for source in orbit_sources:
        twins[coarse_source_key(source, inverses)].append(source)
    keys = sorted(twins)
    if any(len(twins[key]) != n for key in keys):
        raise AssertionError("relative-phase classes must have exact size N")
    representatives = [sorted(twins[key])[0] for key in keys]
    return keys, representatives, permutations, inverses


def two_output_mark_targets(marks, line, source, permutations, inverses, n: int) -> frozenset[tuple]:
    """Release two output slots while keeping both lower bases fixed.

    Targets are ordinary untagged MM/YN tokens.  The topology gate is applied
    after each release by ``target_token``; no provenance label adds capacity.
    """

    replica, coexit, flat = source
    first_base, first_left, first_right, flat_base, flat_left, flat_right = _oriented_source(
        replica, coexit, flat
    )
    targets = set()
    for permutation in permutations:
        second_base = translate_mask(flat_base, permutation)
        output_marks = [
            first_left,
            permutation[flat_left],
            first_right,
            permutation[flat_right],
        ]
        for slots in combinations(range(4), 2):
            for replacements in product(range(n), repeat=2):
                released = output_marks.copy()
                for slot, replacement in zip(slots, replacements):
                    released[slot] = replacement
                target = target_token(
                    marks,
                    line,
                    replica,
                    (first_base, released[0], released[2]),
                    (second_base, released[1], released[3]),
                )
                if target is not None:
                    targets.add(normalize_target(target, inverses))
    return frozenset(targets)


def row_certificate(index: int, row) -> dict[str, Any]:
    n, matrix, _geometry, carrier, marks, line, lower_layer, faces = row
    existing = coarse_row_audit(row, verify_all_twins=False)
    keys, representatives, permutations, inverses = _coarse_classes(row)
    old_neighborhoods = [
        frozenset(
            normalize_target(target, inverses)
            for target in transverse_reservoir_targets(
                marks, line, source, permutations, n, transport=True
            )
        )
        for source in representatives
    ]
    repair_neighborhoods = [
        two_output_mark_targets(marks, line, source, permutations, inverses, n)
        for source in representatives
    ]
    if any(target[0] != "MM" for row_targets in repair_neighborhoods for target in row_targets):
        raise AssertionError("Y=0 gate must leave only ordinary MM repair targets")
    repair = capacitated_hall(keys, repair_neighborhoods, n, "MM")
    old_targets = set().union(*old_neighborhoods)
    repair_targets = set().union(*repair_neighborhoods)
    all_mm_orbits = len(faces["M"]) ** 2 // n
    if len(repair_targets) != all_mm_orbits:
        raise AssertionError("two-mark release must exhaust the ordinary MM orbit reservoir")
    if not repair["saturates"]:
        raise AssertionError("the frozen N9 two-mark repair must saturate")
    return {
        "row_index": index,
        **descriptor(n, matrix, carrier, line, lower_layer, faces),
        "existing_one_carrier_one_mark": existing,
        "two_mark_fixed_base_repair": {
            **repair,
            "ordinary_untagged_targets_only": True,
            "reachable_MM_orbits": len(repair_targets),
            "total_MM_orbits": all_mm_orbits,
            "old_reachable_MM_orbits": len(old_targets),
            "new_MM_orbits_beyond_old_image": len(repair_targets - old_targets),
            "raw_new_MM_tokens": n * len(repair_targets - old_targets),
            "all_class_degrees_equal": len({len(row_targets) for row_targets in repair_neighborhoods}) == 1,
            "class_degree": len(repair_neighborhoods[0]),
        },
    }


def _map_face(face, permutation):
    return (
        translate_mask(face[0], permutation),
        permutation[face[1]],
        permutation[face[2]],
    )


def coordinate_swap_isomorphism(first_row, second_row) -> dict[str, Any]:
    n, matrix, geometry, carrier, _marks, first_line, layer, first_faces = first_row
    _n2, matrix2, geometry2, carrier2, _marks2, second_line, layer2, second_faces = second_row
    if (n, matrix, carrier, layer) != (_n2, matrix2, carrier2, layer2):
        raise AssertionError("strict descriptors differ")
    swap = tuple(geometry.vertex((y, x)) for x, y in geometry.coordinates)
    if tuple(swap[swap[site]] for site in range(n)) != tuple(range(n)):
        raise AssertionError("coordinate swap must be an involution")
    face_checks = {
        motif: {_map_face(face, swap) for face in first_faces[motif]} == set(second_faces[motif])
        for motif in "DMYF"
    }
    translations = set(translation_permutations(geometry))
    conjugates = {
        tuple(swap[permutation[swap[site]]] for site in range(n))
        for permutation in translations
    }
    translation_check = conjugates == translations
    if not all(face_checks.values()) or not translation_check:
        raise AssertionError("coordinate swap failed exact face/translation equivariance")
    return {
        "strict_descriptor_row_indices": [1, 3],
        "line_map": [list(first_line), list(second_line)],
        "site_permutation": list(swap),
        "involution": True,
        "face_family_bijections": face_checks,
        "translation_group_conjugated_to_itself": translation_check,
        "replica_map": "identity",
        "conclusion": (
            "The site coordinate swap x<->y maps every D/M/Y/F face and every translation phase "
            "bijectively, so both the old and two-mark compatibility graphs are exactly isomorphic."
        ),
    }


def build_result() -> dict[str, Any]:
    selected = selected_rows()
    certificates = [row_certificate(index, row) for index, row in selected]
    isomorphism = coordinate_swap_isomorphism(selected[0][1], selected[1][1])
    first = certificates[0]
    old = first["existing_one_carrier_one_mark"]["channel_flows"]["combined"]
    repair = first["two_mark_fixed_base_repair"]
    if (old["total_demand"], old["maximum_flow"], old["Hall_deficiency"]) != (6912, 4752, 2160):
        raise AssertionError("frozen N9 obstruction changed")
    if old["minimum_cut_certificate"]["class_count"] != 768:
        raise AssertionError("the minimum cut must contain every coarse class")
    if repair["new_MM_orbits_beyond_old_image"] != 15984 or repair["class_degree"] != 216:
        raise AssertionError("frozen two-mark reservoir count changed")
    return {
        "schema": SCHEMA,
        "parent_commit": "e2b489493e3ea064f38156b64ea7bf3f4ed4cde3",
        "scope": "only the strict N9 row1/row3 reservoir obstruction; no larger-HNF scan",
        "obstruction": {
            "carrier": "matching",
            "lower_layer": 4,
            "motifs": FAILURE_MOTIFS,
            "mechanism": (
                "Y=0 deletes the YN channel exactly. The surviving one-carrier/one-mark builder "
                "reaches only 4752 ordinary MM target orbits for 6912 units of coarse demand."
            ),
            "minimum_cut": old["minimum_cut_certificate"],
            "source_demand": 6912,
            "reachable_targets": 4752,
            "deficiency": 2160,
            "deficiency_fraction": "5/16",
            "classification": "all-site image-capacity obstruction, not a localized bad-class cut",
        },
        "minimal_legal_repair": {
            "move": (
                "keep both lower bases fixed and release two of the four ordered output-mark slots "
                "to arbitrary quotient sites; reapply the fixed-line topology gate"
            ),
            "capacity_rule": "ordinary untagged MM targets only; no slot, phase, source, or provenance decoration",
            "minimality": (
                "The existing builder strictly contains every fixed-base <=1-mark release and still fails. "
                "Exactly two released slots with fixed bases saturate, so two is minimal on the mark-release axis."
            ),
            "coarse_matching": [repair["maximum_flow"], repair["total_demand"]],
            "reachable_MM_orbits": repair["reachable_MM_orbits"],
            "all_MM_orbits_formula": "M^2/N = 432^2/9 = 20736",
            "new_MM_orbits": repair["new_MM_orbits_beyond_old_image"],
            "raw_new_MM_tokens": repair["raw_new_MM_tokens"],
            "regular_degree": repair["class_degree"],
            "proof": (
                "The exact capacitated max flow is integral and saturates 6912/6912. Translation-orbit "
                "lifting therefore gives a collision-free raw injection. The builder reaches the entire "
                "ordinary MM orbit reservoir without adding target multiplicity."
            ),
        },
        "same_descriptor_isomorphism": isomorphism,
        "rows": certificates,
        "scientific_card": {
            "question": "Why does the corrected combined reservoir first fail at N9, and what is the smallest local repair?",
            "answer": "Y=0 removes synergy and the one-mark MM image has a 5/16 all-site Hall deficit.",
            "repair": "A second output-mark release, with bases fixed and no decorated capacity, reaches every MM orbit and saturates.",
            "new_capacity": "15984 previously unreachable MM orbits (143856 raw tokens), all already present in M^2.",
            "boundary": "Exact for the coordinate-swap-isomorphic N9 rows 1 and 3; not an arbitrary-HNF theorem.",
        },
    }


def render_markdown(result: dict[str, Any]) -> str:
    obstruction = result["obstruction"]
    repair = result["minimal_legal_repair"]
    lines = [
        "# P334 N9 reservoir obstruction and two-mark repair",
        "",
        "## Exact obstruction",
        "",
        (
            "Rows `1` and `3` have `(D,M,Y,F)=(216,432,0,72)` in the matching carrier at layer 4. "
            "Because `Y=0`, the `YN` target channel is empty. The existing one-carrier plus one-output-mark "
            f"reservoir has coarse demand `{obstruction['source_demand']}`, reaches `{obstruction['reachable_targets']}` "
            f"ordinary `MM` targets, and leaves deficiency `{obstruction['deficiency']}={obstruction['deficiency_fraction']}`."
        ),
        "",
        (
            "The residual minimum cut contains all 768 coarse classes (192 from each source replica). "
            "This is an all-site image-capacity obstruction, not a small exceptional family."
        ),
        "",
        "## Minimal legal repair",
        "",
        repair["move"] + ".",
        "",
        (
            "The output remains an ordinary untagged `M x M` pair: no released-slot, phase, source, or provenance "
            "label is retained, and every output face is reclassified on the frozen projective line."
        ),
        "",
        (
            f"The repair reaches `{repair['reachable_MM_orbits']}` normalized targets, exactly "
            f"`{repair['all_MM_orbits_formula']}`. This adds `{repair['new_MM_orbits']}` target orbits "
            f"(`{repair['raw_new_MM_tokens']}` raw tokens) beyond the old image. Every coarse class has degree "
            f"`{repair['regular_degree']}`, and exact integral max flow is `{repair['coarse_matching'][0]}/"
            f"{repair['coarse_matching'][1]}`."
        ),
        "",
        (
            "This is minimal along the output-mark axis: the failed existing builder already strictly contains "
            "all fixed-base zero/one-mark releases, while two released slots with no base mutation saturate."
        ),
        "",
        "## Same-descriptor isomorphism",
        "",
        (
            "The strict descriptor (including `3x3` HNF/Smith type but omitting the line label) selects exactly rows "
            "1 and 3. The explicit site permutation induced by `(x,y)->(y,x)` is an involution, swaps lines "
            "`(0,1)` and `(1,0)`, maps every D/M/Y/F face family bijectively, and conjugates the translation group "
            "to itself. Hence both old and repaired compatibility graphs are exactly isomorphic."
        ),
        "",
        "## Scientific card",
        "",
        f"- **Question:** {result['scientific_card']['question']}",
        f"- **Answer:** {result['scientific_card']['answer']}",
        f"- **Repair:** {result['scientific_card']['repair']}",
        f"- **New capacity:** {result['scientific_card']['new_capacity']}",
        f"- **Boundary:** {result['scientific_card']['boundary']}",
        "",
    ]
    return "\n".join(lines)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args(argv)
    result = build_result()
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.markdown.write_text(render_markdown(result), encoding="utf-8")
    print(f"wrote {args.json}")
    print(f"wrote {args.markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
