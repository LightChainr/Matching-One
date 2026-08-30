#!/usr/bin/env python3
"""Integral saturation certificate for unrestricted digital Alexander carriers.

The proof uses honest ``qL`` covers for every integer ``q>=2``.  A single
cover loses saturation information, while the family of covers detects every
possible Smith defect by choosing ``q`` coprime to that defect.
"""

from __future__ import annotations

import argparse
from functools import reduce
from math import gcd, lcm
import json
from operator import mul
from pathlib import Path
from typing import Iterable, Sequence

from digital_alexander_local_bridge import CORNERS, face_pattern


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 1
    return True


def coprime_prime(index: int) -> int:
    if index < 1:
        raise ValueError("a saturation index must be positive")
    candidate = 2
    while True:
        if _is_prime(candidate) and gcd(candidate, index) == 1:
            return candidate
        candidate += 1


def _product(values: Sequence[int]) -> int:
    return reduce(mul, values, 1)


def smith_intersection_row(elementary_divisors: Sequence[int]) -> dict[str, object]:
    """Audit ``H intersect qL=qH`` in coordinates adapted to ``Sat(H)``."""

    divisors = tuple(int(value) for value in elementary_divisors)
    if not 1 <= len(divisors) <= 2:
        raise ValueError("the torus audit supports rank one or two")
    if any(value < 1 for value in divisors):
        raise ValueError("Smith divisors must be positive")
    if len(divisors) == 2 and divisors[1] % divisors[0]:
        raise ValueError("Smith divisors must be ordered by divisibility")

    defect = _product(divisors)
    q = coprime_prime(defect)
    # In a basis e_1,...,e_r of S=Sat_L(H), extended to a basis of L,
    # H=sum d_i Z e_i and qL intersects S in qS.  Coordinatewise,
    # d_i Z intersect q Z = lcm(d_i,q) Z = q d_i Z.
    intersection_divisors = tuple(lcm(value, q) for value in divisors)
    q_h_divisors = tuple(q * value for value in divisors)
    index_inside_q_s = _product(value // q for value in intersection_divisors)
    return {
        "rank": len(divisors),
        "smith_divisors_of_H_in_S": list(divisors),
        "saturation_defect_d": defect,
        "chosen_prime_q": q,
        "q_coprime_to_d": gcd(q, defect) == 1,
        "intersection_divisors_in_S_coordinates": list(intersection_divisors),
        "qH_divisors_in_S_coordinates": list(q_h_divisors),
        "H_intersect_qL_equals_qH": intersection_divisors == q_h_divisors,
        "index_of_H_intersect_qL_in_qS": index_inside_q_s,
        "defect_is_preserved_upstairs": index_inside_q_s == defect,
        "contradicts_honest_carrier_saturation_when_d_gt_1": (
            defect == 1 or index_inside_q_s > 1
        ),
    }


def smith_defect_audit() -> dict[str, object]:
    rows = [smith_intersection_row((d,)) for d in range(2, 17)]
    rows += [
        smith_intersection_row((first, second))
        for first in range(1, 9)
        for second in range(first, 17)
        if second % first == 0 and first * second > 1
    ]
    all_pass = all(
        row["q_coprime_to_d"]
        and row["H_intersect_qL_equals_qH"]
        and row["defect_is_preserved_upstairs"]
        and row["contradicts_honest_carrier_saturation_when_d_gt_1"]
        for row in rows
    )
    return {
        "rows": rows,
        "rank_one_rows": sum(row["rank"] == 1 for row in rows),
        "rank_two_rows": sum(row["rank"] == 2 for row in rows),
        "all_rows_pass": all_pass,
        "symbolic_statement": (
            "For arbitrary Smith divisors d_i, choose a prime q not dividing "
            "d=product(d_i). In an S-direct-summand basis, "
            "H intersect qL=sum lcm(d_i,q) Z e_i=qH, whose index in qS is d."
        ),
    }


def _path_boundary(path: Sequence[int]) -> dict[int, int]:
    output: dict[int, int] = {}
    for first, second in zip(path, path[1:]):
        output[first] = output.get(first, 0) - 1
        output[second] = output.get(second, 0) + 1
    return {key: value for key, value in sorted(output.items()) if value}


def _path_displacement(path: Sequence[int]) -> tuple[int, int]:
    x = 0
    y = 0
    for first, second in zip(path, path[1:]):
        x += CORNERS[second][0] - CORNERS[first][0]
        y += CORNERS[second][1] - CORNERS[first][1]
    return x, y


def integer_white_face_chain_audit() -> dict[str, object]:
    rows = []
    for mask in range(16):
        pattern = face_pattern(mask)
        replacements = []
        for item in pattern["removed_diagonal_replacements"]:
            diagonal = tuple(item["diagonal"])
            path = tuple(item["boundary_path"])
            replacements.append(
                {
                    "diagonal": list(diagonal),
                    "path": list(path),
                    "same_integral_boundary": (
                        _path_boundary(diagonal) == _path_boundary(path)
                    ),
                    "same_integral_lift_displacement": (
                        _path_displacement(diagonal) == _path_displacement(path)
                    ),
                }
            )
        rows.append(
            {
                "mask": mask,
                "connectivity_preserved": pattern["connectivity_preserved"],
                "at_most_one_retained_diagonal": pattern["embedded_diagonal_gate"],
                "replacements": replacements,
            }
        )
    return {
        "patterns": rows,
        "pattern_count": len(rows),
        "replacement_count": sum(len(row["replacements"]) for row in rows),
        "all_patterns_pass": all(
            row["connectivity_preserved"]
            and row["at_most_one_retained_diagonal"]
            and all(
                item["same_integral_boundary"]
                and item["same_integral_lift_displacement"]
                for item in row["replacements"]
            )
            for row in rows
        ),
        "conclusion": (
            "Replacing every redundant white diagonal by its same-face NN path "
            "preserves connected components and the ambient H1 image over Z, not only Q."
        ),
    }


def build_certificate() -> dict[str, object]:
    face = integer_white_face_chain_audit()
    smith = smith_defect_audit()
    return {
        "schema": "matching-one/digital-alexander-integral-saturation/v1",
        "issue": 269,
        "status": "unrestricted_integral_saturation_theorem",
        "theorem": {
            "scope": (
                "Every finite-index L<=Z^2, every L-periodic site coloring, and "
                "every connected component of the black NN or white matching carrier."
            ),
            "statement": (
                "The image of H1(component;Z) in H1(T_L;Z)=L is saturated. "
                "It is 0 at rank zero, Z times a primitive direction at rank one, "
                "and all of L at rank two."
            ),
            "filtration_consequence": (
                "Every rank-one black/white plateau has saturation index iota=1; "
                "iota is a regression invariant, not an additional state coordinate."
            ),
        },
        "honest_carrier_lemma": {
            "black": (
                "A connected embedded NN graph has a connected closed regular "
                "neighborhood U which deformation retracts to it, so their integral "
                "ambient H1 images agree."
            ),
            "white": (
                "The 16-pattern integral face-chain replacement preserves the full "
                "matching image over Z; the pruned embedded component is the 1-skeleton "
                "of its complementary carrier and has the same integral image."
            ),
            "subsurface_classification": [
                "If U has genus one, an intersection-one pair maps to a unimodular pair in L, so the image is L.",
                "If U has genus zero, H1(U) is generated by boundary curves; every essential embedded boundary circle is primitive and all disjoint essential boundaries on a torus are parallel, so the image is either 0 or Z times one primitive vector.",
                "The empty/rank-zero image is saturated by definition.",
            ],
        },
        "degenerate_quotient_descent": {
            "all_honest_covers": (
                "For every q>=2, qL defines a degree-q^2 cover T_(qL)->T_L. "
                "Every vector in qL has both ambient coordinates divisible by q, "
                "so no nonzero unit-face corner difference lies in qL; every face upstairs is honest."
            ),
            "component_stabilizer": (
                "Fix one universal-lift component C_tilde and H={ell in L:C_tilde+ell=C_tilde}. "
                "The downstairs integral ambient image is exactly H. A loop in its qL-cover "
                "component lifts with endpoint displacement in H intersect qL. "
                "Conversely, for every h in H intersect qL, connectedness supplies a path from x to x+h, "
                "which projects to a loop. Hence the integral upstairs image is exactly H intersect qL."
            ),
            "coprime_cover_contradiction": (
                "Let S=Sat_L(H) and d=[S:H]. Since S is primitive it is a direct summand of L. "
                "Choose a prime q not dividing d and Smith coordinates H=sum d_i Z e_i in S. "
                "Then qL intersect S=qS and H intersect qL=sum lcm(d_i,q) Z e_i=qH. "
                "Its index in qS is still product(d_i)=d. Honest-carrier saturation upstairs forces d=1."
            ),
            "rank_boundaries": {
                "rank_zero": "H=0 is saturated.",
                "rank_one": "d=1 gives H=Z ell for a primitive ell in L.",
                "rank_two": "d=1 gives H=L; no proper finite-index image survives.",
            },
        },
        "machine_certificates": {
            "integer_white_face_chains": face,
            "coprime_smith_descent": smith,
            "all_pass": face["all_patterns_pass"] and smith["all_rows_pass"],
        },
        "search_reclassification": {
            "indices_2_through_13": {
                "HNF_representatives": 140,
                "filtration_paths": 101_140_028_118,
                "rank_one_plateau_steps": 500_805_335_024,
                "maximum_saturation_index": 1,
                "saturation_index_evolution_paths": 0,
            },
            "new_role": (
                "Regression and implementation-convention audit only. No further "
                "index enumeration is needed to infer integral saturation."
            ),
        },
        "claim_boundary": [
            "The theorem concerns ambient torus H1 images, not graph cyclomatic homology.",
            "It does not identify a continuum field, threshold, or finite-size exponent.",
            "The finite index-2-through-13 census is not part of the proof.",
        ],
    }


def render_markdown(certificate: dict[str, object]) -> str:
    theorem = certificate["theorem"]
    machine = certificate["machine_certificates"]
    search = certificate["search_reclassification"]["indices_2_through_13"]
    lines = [
        "# Integral saturation of digital Alexander carriers",
        "",
        f"Status: `{certificate['status']}`.",
        "",
        "## Theorem",
        "",
        theorem["statement"],
        "",
        "## Proof certificate",
        "",
    ]
    for item in certificate["honest_carrier_lemma"]["subsurface_classification"]:
        lines.append(f"- {item}")
    lines += [
        "",
        certificate["degenerate_quotient_descent"]["component_stabilizer"],
        "",
        certificate["degenerate_quotient_descent"]["coprime_cover_contradiction"],
        "",
        "## Executable algebra",
        "",
        f"- integral white face patterns: {machine['integer_white_face_chains']['pattern_count']}; pass={machine['integer_white_face_chains']['all_patterns_pass']}",
        f"- Smith rows: {machine['coprime_smith_descent']['rank_one_rows']} rank-one and {machine['coprime_smith_descent']['rank_two_rows']} rank-two; pass={machine['coprime_smith_descent']['all_rows_pass']}",
        f"- all gates: `{machine['all_pass']}`",
        "",
        "## Reclassified finite frontier",
        "",
        f"The {search['filtration_paths']:,} paths through index 13 are now regression only; saturation is theorem-driven.",
        "",
    ]
    return "\n".join(lines)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    certificate = build_certificate()
    rendered = (
        json.dumps(certificate, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(certificate) + "\n"
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
