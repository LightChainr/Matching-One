#!/usr/bin/env python3
"""Configuration-level obstruction and minimal transport repair for TM.

The aggregate four-face identity leaves one negative token type, ``D x F``.
This oracle asks whether a literal cross-switch of the two ordered faces can
inject the four copies of ``D x F`` into ``M x M`` or four copies of
``Y x nonD``.  It exhausts the first quotient on which the hard type exists.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from itertools import combinations
import json
from pathlib import Path

from p334_dual_hazard_ulc import _local_degrees
from p334_lorentzian_support_gate import _honest_geometries
from projective_essential_birth_oracle import subset_marks


Face = tuple[int, int, int]


def face_type(marks, line, face: Face) -> str | None:
    base, left, right = face
    if (
        marks[base][0] != 1
        or marks[base][1] != line
        or left == right
        or base >> left & 1
        or base >> right & 1
    ):
        return None
    left_rank = marks[base | (1 << left)][0]
    right_rank = marks[base | (1 << right)][0]
    top_rank = marks[base | (1 << left) | (1 << right)][0]
    if left_rank == right_rank == 2:
        return "D"
    if left_rank == 2 or right_rank == 2:
        return "M"
    if top_rank == 2:
        return "Y"
    return "F"


def ordered_faces(marks, line, layer, n: int):
    faces = defaultdict(list)
    for base in layer:
        for left in range(n):
            for right in range(n):
                face = base, left, right
                motif = face_type(marks, line, face)
                if motif is not None:
                    faces[motif].append(face)
    return faces


def translation_permutations(geometry):
    permutations = []
    origin = geometry.coordinates[0]
    for target in geometry.coordinates:
        dx, dy = target[0] - origin[0], target[1] - origin[1]
        permutations.append(
            tuple(
                geometry.vertex((x + dx, y + dy))
                for x, y in geometry.coordinates
            )
        )
    return permutations


def translate_mask(mask: int, permutation) -> int:
    return sum(
        1 << permutation[site]
        for site in range(len(permutation))
        if mask >> site & 1
    )


def balanced_crossovers(left: int, right: int, n: int):
    """All size-preserving recombinations with the same sitewise union."""

    common = left & right
    difference = left ^ right
    sites = [site for site in range(n) if difference >> site & 1]
    need = left.bit_count() - common.bit_count()
    for selected in combinations(sites, need):
        new_left = common | sum(1 << site for site in selected)
        yield new_left, (left | right) ^ new_left


def one_site_mutations(mask: int, n: int, *, include_identity: bool):
    if include_identity:
        yield mask
    for occupied in range(n):
        if not mask >> occupied & 1:
            continue
        for vacant in range(n):
            if mask >> vacant & 1:
                continue
            yield mask ^ (1 << occupied) ^ (1 << vacant)


def target_token(marks, line, replica: int, first: Face, second: Face):
    first_type = face_type(marks, line, first)
    second_type = face_type(marks, line, second)
    if first_type == second_type == "M":
        return "MM", first, second
    if first_type == "Y" and second_type not in (None, "D"):
        return "YN", replica, first, second
    if second_type == "Y" and first_type not in (None, "D"):
        return "YN", replica, second, first
    return None


def _oriented_source(replica: int, coexit: Face, flat: Face):
    source_base, left, right = coexit
    flat_base, flat_left, flat_right = flat
    if replica & 1:
        left, right = right, left
    if replica & 2:
        flat_left, flat_right = flat_right, flat_left
    return source_base, left, right, flat_base, flat_left, flat_right


def alignment_outcome(marks, line, source, permutations):
    """Force the first flat mark onto the first coexit mark."""

    replica, coexit, flat = source
    source_base, left, right, flat_base, flat_left, flat_right = (
        _oriented_source(replica, coexit, flat)
    )
    permutation = next(
        row for row in permutations if row[flat_left] == left
    )
    translated_base = translate_mask(flat_base, permutation)
    translated_left = permutation[flat_left]
    translated_right = permutation[flat_right]
    first = source_base, left, translated_right
    second = translated_base, translated_left, right
    return face_type(marks, line, first), face_type(marks, line, second)


def cross_targets(
    marks,
    line,
    source,
    permutations,
    n: int,
    operation: str,
):
    replica, coexit, flat = source
    source_base, left, right, flat_base, flat_left, flat_right = (
        _oriented_source(replica, coexit, flat)
    )
    targets = set()
    for permutation in permutations:
        translated_base = translate_mask(flat_base, permutation)
        translated_left = permutation[flat_left]
        translated_right = permutation[flat_right]
        if operation == "mark_only":
            base_pairs = ((source_base, translated_base),)
        elif operation == "union_preserving":
            base_pairs = balanced_crossovers(source_base, translated_base, n)
        elif operation == "one_carrier_transport":
            base_pairs = set(
                (new_base, translated_base)
                for new_base in one_site_mutations(
                    source_base, n, include_identity=True
                )
            ) | set(
                (source_base, new_base)
                for new_base in one_site_mutations(
                    translated_base, n, include_identity=True
                )
            )
        elif operation == "two_carrier_transport":
            base_pairs = (
                (new_left, new_right)
                for new_left in one_site_mutations(
                    source_base, n, include_identity=False
                )
                for new_right in one_site_mutations(
                    translated_base, n, include_identity=False
                )
            )
        else:
            raise ValueError(f"unknown cross-switch operation: {operation}")
        for first_base, second_base in base_pairs:
            first = first_base, left, translated_left
            second = second_base, right, translated_right
            target = target_token(marks, line, replica, first, second)
            if target is not None:
                targets.add(target)
    return sorted(targets)


def maximum_matching(adjacency):
    """Exact deterministic bipartite matching for the bounded oracle."""

    matched_target = {}

    def augment(source_index: int, seen: set) -> bool:
        for target in adjacency[source_index]:
            if target in seen:
                continue
            seen.add(target)
            previous = matched_target.get(target)
            if previous is None or augment(previous, seen):
                matched_target[target] = source_index
                return True
        return False

    flow = 0
    for source_index in sorted(
        range(len(adjacency)), key=lambda index: (len(adjacency[index]), index)
    ):
        flow += augment(source_index, set())
    return flow


def operation_audit(marks, line, sources, permutations, n: int, operation: str):
    adjacency = [
        cross_targets(marks, line, source, permutations, n, operation)
        for source in sources
    ]
    targets = {target for row in adjacency for target in row}
    flow = maximum_matching(adjacency)
    degrees = [len(row) for row in adjacency]
    return {
        "operation": operation,
        "source_tokens": len(sources),
        "reachable_cover_tokens": len(targets),
        "maximum_matching": flow,
        "Hall_deficiency": len(sources) - flow,
        "zero_degree_sources": sum(degree == 0 for degree in degrees),
        "minimum_degree": min(degrees),
        "maximum_degree": max(degrees),
        "collision_free_injection_exists": flow == len(sources),
    }


def phase_reconstruction_audit(marks, line, sources, permutations, n: int):
    bare = Counter()
    phase = Counter()
    phase_replica = Counter()
    admissible_phases = Counter()
    outcomes = Counter()
    for source in sources:
        replica, coexit, flat = source
        source_base, left, right, flat_base, flat_left, flat_right = (
            _oriented_source(replica, coexit, flat)
        )
        rows = []
        for phase_index, permutation in enumerate(permutations):
            translated_base = translate_mask(flat_base, permutation)
            first = source_base, left, permutation[flat_left]
            second = translated_base, right, permutation[flat_right]
            target = target_token(marks, line, replica, first, second)
            if target is not None:
                rows.append((phase_index, target))
        admissible_phases[len(rows)] += 1
        for phase_index, target in rows:
            outcomes[target[0]] += 1
            bare[target] += 1
            phase[phase_index, target] += 1
            phase_replica[phase_index, replica, target] += 1

    def summary(counter):
        return {
            "distinct_images": len(counter),
            "maximum_fiber": max(counter.values()),
            "collision_excess": sum(counter.values()) - len(counter),
        }

    return {
        "admissible_phase_count_histogram": dict(admissible_phases),
        "output_occurrences": dict(outcomes),
        "bare_target": summary(bare),
        "target_plus_translation_phase": summary(phase),
        "target_plus_phase_and_replica": summary(phase_replica),
    }


def build_result():
    rows = []
    lower_order_hard_rows = 0
    for n, matrix, geometry in _honest_geometries(6):
        for carrier, marks in (
            ("primal", subset_marks(geometry, matching=False)),
            ("matching", subset_marks(geometry, matching=True)),
        ):
            lines = sorted({line for rank, line, _ in marks if rank == 1})
            for line in lines:
                layers, _ = _local_degrees(marks, line, n)
                for lower_layer in range(n):
                    if not layers[lower_layer]:
                        continue
                    faces = ordered_faces(
                        marks, line, layers[lower_layer], n
                    )
                    if not faces["D"] or not faces["F"]:
                        continue
                    if n < 6:
                        lower_order_hard_rows += 1
                        continue
                    translations = translation_permutations(geometry)
                    sources = [
                        (replica, coexit, flat)
                        for replica in range(4)
                        for coexit in faces["D"]
                        for flat in faces["F"]
                    ]
                    alignment = Counter(
                        alignment_outcome(marks, line, source, translations)
                        for source in sources
                    )
                    audits = [
                        operation_audit(
                            marks,
                            line,
                            sources,
                            translations,
                            n,
                            operation,
                        )
                        for operation in (
                            "mark_only",
                            "union_preserving",
                            "one_carrier_transport",
                            "two_carrier_transport",
                        )
                    ]
                    rows.append(
                        {
                            "N": n,
                            "matrix": [list(row) for row in matrix],
                            "carrier": carrier,
                            "line": list(line),
                            "lower_layer": lower_layer,
                            "motifs": {
                                key: len(faces[key]) for key in "DMYF"
                            },
                            "hard_tokens_4DF": len(sources),
                            "forced_alignment_outcomes": {
                                f"{left}|{right}": count
                                for (left, right), count in sorted(
                                    alignment.items(),
                                    key=lambda item: str(item[0]),
                                )
                            },
                            "forced_alignment_positive_targets": sum(
                                count
                                for (left, right), count in alignment.items()
                                if (left == right == "M")
                                or (left == "Y" and right not in (None, "D"))
                                or (right == "Y" and left not in (None, "D"))
                            ),
                            "phase_information_loss": phase_reconstruction_audit(
                                marks, line, sources, translations, n
                            ),
                            "operation_audits": {
                                audit["operation"]: audit for audit in audits
                            },
                        }
                    )

    assert lower_order_hard_rows == 0
    assert len(rows) == 4
    for row in rows:
        assert row["motifs"] == {"D": 12, "M": 48, "Y": 24, "F": 24}
        assert row["hard_tokens_4DF"] == 1152
        assert row["forced_alignment_positive_targets"] == 0
        phase = row["phase_information_loss"]
        assert phase["admissible_phase_count_histogram"] == {1: 1152}
        assert phase["bare_target"] == {
            "distinct_images": 120,
            "maximum_fiber": 24,
            "collision_excess": 1032,
        }
        assert phase["target_plus_translation_phase"] == {
            "distinct_images": 720,
            "maximum_fiber": 4,
            "collision_excess": 432,
        }
        assert phase["target_plus_phase_and_replica"] == {
            "distinct_images": 1152,
            "maximum_fiber": 1,
            "collision_excess": 0,
        }
        audits = row["operation_audits"]
        assert not audits["mark_only"]["collision_free_injection_exists"]
        assert not audits["union_preserving"]["collision_free_injection_exists"]
        assert not audits["one_carrier_transport"][
            "collision_free_injection_exists"
        ]
        assert not audits["two_carrier_transport"][
            "collision_free_injection_exists"
        ]

    union_deficiencies = sorted(
        {
            row["operation_audits"]["union_preserving"]["Hall_deficiency"]
            for row in rows
        }
    )
    one_carrier_deficiencies = sorted(
        {
            row["operation_audits"]["one_carrier_transport"][
                "Hall_deficiency"
            ]
            for row in rows
        }
    )
    assert union_deficiencies == [1056]
    assert one_carrier_deficiencies == [612]

    result = {
        "schema_version": "p334-tm-configuration-cross-switch-v1",
        "universal_strict_injection_status": "refuted",
        "minimal_counterexample_gate": {
            "search": "all connected honest-face HNF quotients through N=6, both carriers, every fixed line and lower layer",
            "hard_rows_below_N6": lower_order_hard_rows,
            "minimal_N": 6,
            "minimal_rows": len(rows),
        },
        "exact_obstruction": {
            "statement": "A translation-equivariant two-face switch that only crosses ordered missing marks is not injective. Even allowing every balanced recombination of the two bases while preserving the sitewise occupation multiplicity u_i=1_S(i)+1_T(i) leaves a positive Hall deficiency.",
            "forced_alignment": "zero positive targets on all 4608 minimal hard tokens",
            "mark_only_Hall_deficiency": 1032,
            "union_preserving_Hall_deficiencies": union_deficiencies,
            "one_carrier_transport_Hall_deficiencies": one_carrier_deficiencies,
            "lost_information": "The mark-only output forgets the quotient translation phase; M x M also forgets the fourfold source replica. Bare fibers have maximum size 24. Adding phase leaves maximum fiber 4; phase plus replica reconstructs the source but is a decorated certificate, not an unmarked TM injection.",
            "obstruction_invariant": "sitewise union multiplicity u in {0,1,2}^Q of the two lower configurations",
        },
        "failed_two_carrier_candidate": {
            "move": "one occupied-to-vacant replacement independently in each of the two lower configurations, followed by the ordered-mark cross",
            "matching": "588 of 1152 hard tokens on each of four N=6 rows",
            "Hall_deficiency": 564,
            "interpretation": "transporting the bases does not release the crossed ordered missing-mark invariant; an Alexander-dual birth-square reservoir must also contribute a fresh transverse mark",
            "status": "exactly refuted after enforcing that every output base remains in the same fixed-line rank-one stratum",
        },
        "corrected_general_theorem": {
            "statement": "Any configuration compatibility graph from the 4DF hard tokens to M x M plus four replicas of Y x nonD proves aggregate TM if it satisfies Hall.",
            "proved": "the implication from Hall saturation to aggregate TM",
            "open": "identify the minimal Alexander-dual reservoir that releases one crossed missing mark while retaining fixed-line base semantics",
            "important_boundary": "mark-only, union-preserving, one-carrier, and two-carrier base transport all fail on the minimal N=6 rows",
        },
        "production_pair_covariance_crosscheck": {
            "commit": "a9f7d28",
            "assessment": "orthogonal",
            "reason": "its signed union histograms concern finite local embedding controls at fixed K; they do not retain the fixed projective line, global D/M/Y/F rank pattern, relative translation phase, or the two-base sitewise-union invariant that obstructs this switch",
        },
        "rows": rows,
        "scientific_card": {
            "question": "Can D x F be removed by a universal configuration-level cross-switch?",
            "answer": "Not by any two-face mark switch, nor by any sitewise-union-preserving base crossover; N=6 is an exact counterexample.",
            "obstruction": "The switch collapses translation phase and replica information and remains trapped in a fixed two-base union fiber.",
            "minimal_repair": "Base transport alone is not a repair: two-carrier transport reaches only 588/1152.",
            "next_theorem": "Add exactly one fresh transverse output mark through the Alexander-dual birth square and re-test Hall; do not return to base-only switching.",
        },
    }
    return json.loads(json.dumps(result))


def render_markdown(result):
    rows = result["rows"]
    matrices = ", ".join(str(row["matrix"]) for row in rows)
    union_deficiencies = result["exact_obstruction"][
        "union_preserving_Hall_deficiencies"
    ]
    one_deficiencies = result["exact_obstruction"][
        "one_carrier_transport_Hall_deficiencies"
    ]
    return "\n".join(
        [
            "# The D x F cross-switch obstruction and its minimal repair",
            "",
            "## Exact refutation",
            "",
            "The hoped-for universal mark-only injection is false already at `N=6`. The exhaustive minimal gate finds no hard row below `N=6` and four matching-carrier rows at `N=6`, with matrices " + matrices + ". Each has `(D,M,Y,F)=(12,48,24,24)` and therefore `4DF=1152` hard tokens.",
            "",
            "Forcing one flat mark to align with one coexit mark produces no positive target on any of the 4608 hard tokens. Allowing every quotient phase does not cure collisions: each source has exactly one admissible phase, but the unmarked output has only 120 distinct images and maximum fiber 24, so maximum matching is 120 and Hall deficiency is 1032.",
            "",
            "The exact lost labels are visible. Decorating by translation phase gives 720 images and leaves fourfold fibers. Decorating also by the source replica gives 1152 singleton images. That is a reconstruction certificate, not an unmarked TM injection: the decorations would add target capacity.",
            "",
            "## Obstruction invariant",
            "",
            "Let `u_i=1_S(i)+1_T(i)` be the sitewise occupation multiplicity of the two lower configurations after translating the flat face. All balanced configuration crossovers preserve `u`. Exhausting every such crossover still leaves Hall deficiencies " + str(union_deficiencies) + ". Thus the obstruction is not a bad choice of pairing; the entire union fiber is too small.",
            "",
            "Breaking the union invariant in only one carrier is still insufficient: the exact Hall deficiencies are " + str(one_deficiencies) + ". The signed-union covariance certificate at `a9f7d28` is orthogonal because it concerns local embedding overlaps, not the fixed-line global-rank union fiber.",
            "",
            "## Minimal corrected theorem",
            "",
            "Allowing one occupied-to-vacant replacement independently in both lower configurations is still insufficient once every output base is required to remain in the same fixed-line stratum: maximum matching is only `588/1152`, with Hall deficiency `564`, on every minimal row.",
            "",
            "For an arbitrary fixed-line HNF row, Hall saturation of any genuine configuration compatibility graph remains sufficient for aggregate TM. But the base-only two-carrier graph is not that theorem: the Alexander-dual birth square must also release a fresh transverse output mark. The original universal mark-only and base-transport injections should be considered refuted, not merely unfinished.",
            "",
            "## Scientific card",
            "",
            "- **Question:** Can `D x F` be removed by a universal configuration-level cross-switch?",
            "- **Answer:** Not by a mark switch or a sitewise-union-preserving base crossover; `N=6` is an exact counterexample.",
            "- **Obstruction:** Translation phase and source replica are collapsed inside a fixed two-base union fiber.",
            "- **Minimal repair:** Base transport alone fails: even two carriers reach only `588/1152`.",
            "- **Next theorem:** Release one fresh transverse mark through Alexander complement and test the resulting orbit Hall graph.",
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
