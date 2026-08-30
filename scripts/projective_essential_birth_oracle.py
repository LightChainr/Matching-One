#!/usr/bin/env python3
"""Exact Phase-A oracle for the projective mark of the first H1 birth.

The oracle works on the Boolean subset lattice rather than enumerating all
``N!`` permutations.  A dynamic program nevertheless counts every ordering
exactly, retaining both birth times, their sites, and the primitive line born
at the first rank-one transition.  Integral saturation is a theorem here:
``iota`` is recorded as one, not treated as a sampled coordinate.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from math import comb, factorial, gcd
from pathlib import Path
from typing import Any, Iterable, NamedTuple, Optional, Sequence, Tuple

from c4_self_matching_exact import c4_self_matching_torus
from digital_alexander_filtration_oracle import rank_mark
from integer_period_torus import (
    IntegerTorusGeometry,
    axis_integer_torus,
    classify_configuration,
    gaussian_integer_torus,
    integer_torus_geometry,
    matrix_product,
    matrix_vector,
    unimodular_inverse,
)


Vector = Tuple[int, int]
Matrix = Tuple[Tuple[int, int], Tuple[int, int]]
ComplexQ = Tuple[Fraction, Fraction]


class BirthRecord(NamedTuple):
    k1: Optional[int]
    k2: Optional[int]
    line: Optional[Vector]
    site1: Optional[int]
    site2: Optional[int]


def canonical_projective(vector: Vector) -> Vector:
    """Return the primitive representative of ``{v,-v}``."""

    divisor = gcd(abs(vector[0]), abs(vector[1]))
    if divisor == 0:
        raise ValueError("zero vector does not define a projective line")
    x, y = vector[0] // divisor, vector[1] // divisor
    return (-x, -y) if x < 0 or (x == 0 and y < 0) else (x, y)


def chi4(period_matrix: Matrix, line: Vector) -> ComplexQ:
    """Return ``((x+iy)/|x+iy|)^4`` as an exact rational pair."""

    x, y = matrix_vector(period_matrix, line)
    radius2 = x * x + y * y
    if radius2 == 0:
        raise ValueError("a projective line cannot have zero physical image")
    denominator = radius2 * radius2
    return (
        Fraction(x**4 - 6 * x * x * y * y + y**4, denominator),
        Fraction(4 * x * y * (x * x - y * y), denominator),
    )


def _qadd(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] + right[0], left[1] + right[1]


def _qsub(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] - right[0], left[1] - right[1]


def _qpayload(value: ComplexQ) -> dict[str, str]:
    return {"real": str(value[0]), "imag": str(value[1])}


def _active(mask: int, n: int) -> tuple[bool, ...]:
    return tuple(bool(mask & (1 << vertex)) for vertex in range(n))


def subset_marks(
    geometry: IntegerTorusGeometry, *, matching: bool
) -> list[Tuple[int, Optional[Vector], Optional[int]]]:
    """Cache rank/line/index on the full Boolean subset lattice."""

    rows = []
    for mask in range(1 << geometry.n):
        row = rank_mark(geometry, _active(mask, geometry.n), matching=matching)
        if row[0] == 1 and row[2] != 1:
            raise AssertionError(
                "rank-one carrier violated the integral saturation theorem"
            )
        rows.append(row)
    return rows


def exact_path_histogram(
    geometry: IntegerTorusGeometry,
    marks: Sequence[Tuple[int, Optional[Vector], Optional[int]]],
) -> Counter[BirthRecord]:
    """Count all permutations by propagating marked states over subsets."""

    n = geometry.n
    tables: list[Counter[BirthRecord]] = [Counter() for _ in range(1 << n)]
    tables[0][BirthRecord(None, None, None, None, None)] = 1
    for mask in range(1 << n):
        old_rank, old_line, _ = marks[mask]
        for record, count in tables[mask].items():
            if old_rank == 1 and record.line != old_line:
                raise AssertionError("a path mark disagreed with the state line")
            for vertex in range(n):
                bit = 1 << vertex
                if mask & bit:
                    continue
                new_mask = mask | bit
                new_rank, new_line, _ = marks[new_mask]
                if new_rank < old_rank:
                    raise AssertionError("ambient rank decreased under site addition")
                k = mask.bit_count() + 1
                k1, k2, line, site1, site2 = record
                if k1 is None and new_rank >= 1:
                    k1, site1 = k, vertex
                    # A direct 0->2 transition has no projective first line.
                    line = new_line if new_rank == 1 else None
                if k2 is None and new_rank == 2:
                    k2, site2 = k, vertex
                tables[new_mask][BirthRecord(k1, k2, line, site1, site2)] += count
    terminal = tables[-1]
    if sum(terminal.values()) != factorial(n):
        raise AssertionError("subset DP did not recover N! paths")
    if any(row.k1 is None or row.k2 is None for row in terminal):
        raise AssertionError("a full path omitted an essential birth")
    return terminal


def _line_label(line: Optional[Vector]) -> str:
    return "DIRECT_RANK2" if line is None else f"{line[0]},{line[1]}"


def _aggregate_joint(
    histogram: Counter[BirthRecord],
) -> tuple[Counter[tuple[int, int, str]], Counter[int], Counter[tuple[int, str]]]:
    joint: Counter[tuple[int, int, str]] = Counter()
    k1: Counter[int] = Counter()
    marked_k1: Counter[tuple[int, str]] = Counter()
    for record, count in histogram.items():
        label = _line_label(record.line)
        assert record.k1 is not None and record.k2 is not None
        joint[(record.k1, record.k2, label)] += count
        k1[record.k1] += count
        marked_k1[(record.k1, label)] += count
    return joint, k1, marked_k1


def unmarked_k1_histogram(
    marks: Sequence[Tuple[int, Optional[Vector], Optional[int]]], n: int
) -> Counter[int]:
    """Recover K1 without using a line mark, by exact transition weights."""

    result: Counter[int] = Counter()
    for mask, (old_rank, _, _) in enumerate(marks):
        if old_rank != 0:
            continue
        k = mask.bit_count() + 1
        path_weight = factorial(k - 1) * factorial(n - k)
        for vertex in range(n):
            if mask & (1 << vertex):
                continue
            if marks[mask | (1 << vertex)][0] >= 1:
                result[k] += path_weight
    if sum(result.values()) != factorial(n):
        raise AssertionError("unmarked transition census did not recover N! paths")
    return result


def canonicality_audit(
    geometry: IntegerTorusGeometry,
    marks: Sequence[Tuple[int, Optional[Vector], Optional[int]]],
    *,
    matching: bool,
) -> dict[str, int | bool]:
    failures = 0
    rank_one_states = 0
    generator_checks = 0
    for mask, (rank, line, index) in enumerate(marks):
        if rank != 1:
            continue
        rank_one_states += 1
        assert line is not None
        if index != 1 or canonical_projective(line) != line:
            failures += 1
        if canonical_projective((-line[0], -line[1])) != line:
            failures += 1
        _, components = classify_configuration(
            geometry, _active(mask, geometry.n), matching=matching
        )
        for component in components:
            for generator in component.generators:
                if generator == (0, 0):
                    continue
                generator_checks += 1
                if canonical_projective(generator) != line:
                    failures += 1
    return {
        "rank_one_states": rank_one_states,
        "generator_checks": generator_checks,
        "up_to_sign_failures": failures,
        "all_pass": failures == 0,
    }


def complement_audit(
    geometry: IntegerTorusGeometry,
    black_marks: Sequence[Tuple[int, Optional[Vector], Optional[int]]],
    matching_marks: Sequence[Tuple[int, Optional[Vector], Optional[int]]],
    black_paths: Counter[BirthRecord],
    matching_paths: Counter[BirthRecord],
) -> dict[str, int | bool]:
    full = (1 << geometry.n) - 1
    state_failures = 0
    rank_one_pairs = 0
    for mask, (black_rank, black_line, black_index) in enumerate(black_marks):
        white_rank, white_line, white_index = matching_marks[full ^ mask]
        if black_rank + white_rank != 2:
            state_failures += 1
        if black_rank == 1:
            rank_one_pairs += 1
            if (
                white_rank != 1
                or black_line != white_line
                or black_index != 1
                or white_index != 1
            ):
                state_failures += 1

    reflected: Counter[BirthRecord] = Counter()
    n = geometry.n
    for record, count in black_paths.items():
        assert record.k1 is not None and record.k2 is not None
        reflected[
            BirthRecord(
                n + 1 - record.k2,
                n + 1 - record.k1,
                record.line,
                record.site2,
                record.site1,
            )
        ] += count
    path_failures = sum((reflected - matching_paths).values()) + sum(
        (matching_paths - reflected).values()
    )
    return {
        "rank_one_complement_pairs": rank_one_pairs,
        "configuration_failures": state_failures,
        "reflected_path_count_failures": path_failures,
        "all_pass": state_failures == 0 and path_failures == 0,
    }


def basis_covariance_audit(geometry: IntegerTorusGeometry) -> dict[str, Any]:
    changes: tuple[Matrix, ...] = (
        ((0, -1), (1, 0)),
        ((1, 1), (0, 1)),
        ((-1, 0), (0, 1)),
    )
    base_marks = subset_marks(geometry, matching=False)
    rows = []
    for change in changes:
        changed = integer_torus_geometry(
            matrix_product(geometry.periods.matrix, change),
            name=f"{geometry.name}-rebased",
        )
        changed_marks = subset_marks(changed, matching=False)
        inverse = unimodular_inverse(change)
        failures = 0
        checked = 0
        for mask, (rank, line, index) in enumerate(base_marks):
            changed_mask = 0
            for vertex, point in enumerate(geometry.coordinates):
                if mask & (1 << vertex):
                    changed_mask |= 1 << changed.vertex(point)
            new_rank, new_line, new_index = changed_marks[changed_mask]
            if new_rank != rank:
                failures += 1
                continue
            if rank == 1:
                checked += 1
                assert line is not None
                expected = canonical_projective(matrix_vector(inverse, line))
                if new_line != expected or index != 1 or new_index != 1:
                    failures += 1
        rows.append(
            {
                "change": [list(row) for row in change],
                "rank_one_states_checked": checked,
                "failures": failures,
            }
        )
    return {"changes": rows, "all_pass": all(row["failures"] == 0 for row in rows)}


def _d4_matrices() -> tuple[tuple[str, Matrix, str], ...]:
    identity: Matrix = ((1, 0), (0, 1))
    rotation: Matrix = ((0, -1), (1, 0))
    reflection: Matrix = ((1, 0), (0, -1))
    rotations = [identity]
    for _ in range(3):
        rotations.append(matrix_product(rotation, rotations[-1]))
    return tuple(
        [(f"R{90 * turn}", value, "identity") for turn, value in enumerate(rotations)]
        + [
            (f"R{90 * turn}X", matrix_product(value, reflection), "conjugation")
            for turn, value in enumerate(rotations)
        ]
    )


def d4_chi4_audit(geometry: IntegerTorusGeometry) -> dict[str, Any]:
    marks = subset_marks(geometry, matching=False)
    period = geometry.periods.matrix
    rows = []
    for name, action, chi_action in _d4_matrices():
        failures = 0
        checked = 0
        for mask, (rank, line, _) in enumerate(marks):
            if rank != 1:
                continue
            assert line is not None
            transformed_mask = 0
            for vertex, point in enumerate(geometry.coordinates):
                if mask & (1 << vertex):
                    transformed_mask |= 1 << geometry.vertex(matrix_vector(action, point))
            new_rank, new_line, new_index = marks[transformed_mask]
            physical = matrix_vector(period, line)
            transformed_physical = matrix_vector(action, physical)
            expected_line = canonical_projective(
                geometry.periods.winding(transformed_physical)
            )
            before = chi4(period, line)
            after = chi4(period, expected_line)
            expected_chi = before if chi_action == "identity" else (before[0], -before[1])
            checked += 1
            if (
                new_rank != 1
                or new_line != expected_line
                or new_index != 1
                or after != expected_chi
            ):
                failures += 1
        rows.append(
            {
                "element": name,
                "chi4_action": chi_action,
                "rank_one_states_checked": checked,
                "failures": failures,
            }
        )
    return {"elements": rows, "all_pass": all(row["failures"] == 0 for row in rows)}


def primitive_sector_crosswalk(
    geometry: IntegerTorusGeometry,
    marks: Sequence[Tuple[int, Optional[Vector], Optional[int]]],
) -> dict[str, Any]:
    """Verify A4 is #156's state character and dA4=j_birth-j_exit."""

    n = geometry.n
    state_crosswalk_failures = 0
    # A path which reaches a rank-one state must carry its unique state line.
    path_tables: list[Counter[Optional[Vector]]] = [Counter() for _ in range(1 << n)]
    path_tables[0][None] = 1
    for mask in range(1 << n):
        rank, line, _ = marks[mask]
        for path_line, count in path_tables[mask].items():
            if rank == 1 and path_line != line:
                state_crosswalk_failures += count
            for vertex in range(n):
                bit = 1 << vertex
                if mask & bit:
                    continue
                new_rank, new_line, _ = marks[mask | bit]
                carried = path_line
                if rank == 0 and new_rank == 1:
                    carried = new_line
                elif new_rank == 2:
                    carried = None
                path_tables[mask | bit][carried] += count

    derivative_rows = []
    derivative_failures = 0
    for k in range(n):
        direct: ComplexQ = (Fraction(0), Fraction(0))
        births: ComplexQ = (Fraction(0), Fraction(0))
        exits: ComplexQ = (Fraction(0), Fraction(0))
        edge_count = 0
        for mask in range(1 << n):
            if mask.bit_count() != k:
                continue
            old_rank, old_line, _ = marks[mask]
            old_value = (
                chi4(geometry.periods.matrix, old_line)
                if old_rank == 1 and old_line is not None
                else (Fraction(0), Fraction(0))
            )
            for vertex in range(n):
                if mask & (1 << vertex):
                    continue
                edge_count += 1
                new_rank, new_line, _ = marks[mask | (1 << vertex)]
                new_value = (
                    chi4(geometry.periods.matrix, new_line)
                    if new_rank == 1 and new_line is not None
                    else (Fraction(0), Fraction(0))
                )
                direct = _qadd(direct, _qsub(new_value, old_value))
                if old_rank == 0 and new_rank == 1:
                    births = _qadd(births, new_value)
                if old_rank == 1 and new_rank == 2:
                    # The exit/source at K2 keeps the line born at K1.
                    exits = _qadd(exits, old_value)
        if direct != _qsub(births, exits):
            derivative_failures += 1
        derivative_rows.append(
            {
                "lower_subset_size": k,
                "directed_edges": edge_count,
                "bernstein_degree_N_minus_1_normalizer": comb(n - 1, k),
                "dA4_raw": _qpayload(direct),
                "j4_birth1_raw": _qpayload(births),
                "j4_exit2_raw": _qpayload(exits),
            }
        )
    return {
        "identity": "A4(state)=primitive-sector chi4; dA4/dp=j4_birth1-j4_exit2",
        "state_crosswalk_path_failures": state_crosswalk_failures,
        "derivative_coefficient_failures": derivative_failures,
        "degree_N_minus_1_raw_coefficients": derivative_rows,
        "all_pass": state_crosswalk_failures == 0 and derivative_failures == 0,
    }


def summarize_geometry(name: str, geometry: IntegerTorusGeometry) -> dict[str, Any]:
    black_marks = subset_marks(geometry, matching=False)
    matching_marks = subset_marks(geometry, matching=True)
    black_paths = exact_path_histogram(geometry, black_marks)
    matching_paths = exact_path_histogram(geometry, matching_marks)
    joint, k1, marked_k1 = _aggregate_joint(black_paths)
    direct_k1 = unmarked_k1_histogram(black_marks, geometry.n)

    marginal_recovered = Counter()
    for (birth, _line), count in marked_k1.items():
        marginal_recovered[birth] += count
    direct_paths = sum(
        count for record, count in black_paths.items() if record.line is None
    )
    line_bearing_k1 = Counter()
    line_counts: Counter[str] = Counter()
    timing_counts: Counter[tuple[int, int]] = Counter()
    site_pair_counts: Counter[tuple[int, int]] = Counter()
    projective_paths = 0
    for record, count in black_paths.items():
        assert record.k1 is not None and record.k2 is not None
        assert record.site1 is not None and record.site2 is not None
        site_pair_counts[(record.site1, record.site2)] += count
        if record.line is not None:
            label = _line_label(record.line)
            projective_paths += count
            line_bearing_k1[record.k1] += count
            line_counts[label] += count
            timing_counts[(record.k1, record.k2)] += count
    timing_line_factorization_failures = 0
    distinct_chi4_by_timing: dict[str, set[ComplexQ]] = {}
    for (birth1, birth2, label), count in joint.items():
        if label == "DIRECT_RANK2":
            continue
        if count * projective_paths != timing_counts[(birth1, birth2)] * line_counts[label]:
            timing_line_factorization_failures += 1
        line = tuple(int(value) for value in label.split(","))
        distinct_chi4_by_timing.setdefault(f"{birth1},{birth2}", set()).add(
            chi4(geometry.periods.matrix, line)  # type: ignore[arg-type]
        )
    supports: dict[tuple[int, int], set[str]] = {}
    for birth1, birth2, line in joint:
        supports.setdefault((birth1, birth2), set()).add(line)
    multi_line_cells = {
        f"{birth1},{birth2}": sorted(lines)
        for (birth1, birth2), lines in supports.items()
        if len(lines) > 1
    }
    return {
        "id": name,
        "N": geometry.n,
        "period_matrix": [list(row) for row in geometry.periods.matrix],
        "all_paths": factorial(geometry.n),
        "subgroup_contract": {
            "after_line_birth": "Z ell with ell primitive",
            "after_second_birth": "Z^2",
            "iota": 1,
        },
        "joint_K1_K2_ell_counts": [
            {"K1": key[0], "K2": key[1], "ell": key[2], "count": count}
            for key, count in sorted(joint.items())
        ],
        "K1_histogram": {str(key): value for key, value in sorted(k1.items())},
        "unmarked_transition_K1_histogram": {
            str(key): value for key, value in sorted(direct_k1.items())
        },
        "sum_over_ell_plus_direct_atom_recovers_K1": (
            marginal_recovered == k1 == direct_k1
        ),
        "sum_over_projective_ell_recovers_line_bearing_K1": (
            sum(line_bearing_k1.values()) == projective_paths
        ),
        "direct_rank2_paths_without_projective_line": direct_paths,
        "terminal_marked_state_count": len(black_paths),
        "birth_site_pair_counts": {
            f"{first},{second}": count
            for (first, second), count in sorted(site_pair_counts.items())
        },
        "timing_cells_with_multiple_line_marks": multi_line_cells,
        "ell_not_determined_by_K1_K2": bool(multi_line_cells),
        "ell_independent_of_K1_K2_conditional_on_line_birth": (
            timing_line_factorization_failures == 0
        ),
        "timing_line_factorization_failures": timing_line_factorization_failures,
        "distinct_chi4_values_per_timing_cell": {
            timing: len(values) for timing, values in sorted(distinct_chi4_by_timing.items())
        },
        "chi4_adds_information_beyond_K1_K2_on_this_quotient": any(
            len(values) > 1 for values in distinct_chi4_by_timing.values()
        ),
        "canonicality": canonicality_audit(
            geometry, black_marks, matching=False
        ),
        "complement_Alexander": complement_audit(
            geometry, black_marks, matching_marks, black_paths, matching_paths
        ),
        "primitive_sector_crosswalk": primitive_sector_crosswalk(
            geometry, black_marks
        ),
    }


def build_certificate() -> dict[str, Any]:
    geometries = (
        ("axis-L2", axis_integer_torus(2)),
        ("gaussian-2-1", gaussian_integer_torus(2, 1)),
        ("c4-self-matching-3-1", c4_self_matching_torus(3, 1)),
    )
    rows = [summarize_geometry(name, geometry) for name, geometry in geometries]
    covariance = {
        name: basis_covariance_audit(geometry)
        for name, geometry in geometries[:2]
    }
    d4 = d4_chi4_audit(geometries[0][1])
    all_pass = (
        all(row["sum_over_ell_plus_direct_atom_recovers_K1"] for row in rows)
        and all(row["canonicality"]["all_pass"] for row in rows)
        and all(row["complement_Alexander"]["all_pass"] for row in rows)
        and all(row["primitive_sector_crosswalk"]["all_pass"] for row in rows)
        and all(row["all_paths"] == factorial(row["N"]) for row in rows)
        and all(value["all_pass"] for value in covariance.values())
        and d4["all_pass"]
    )
    return {
        "schema": "matching-one/projective-essential-birth-oracle/v1",
        "issue": 334,
        "status": "exact_phase_A_certificate",
        "geometries": rows,
        "SL2Z_basis_covariance": covariance,
        "D4_chi4_action": d4,
        "scientific_crosswalk": {
            "A4": "exactly the fixed-p rank-one primitive-sector chi4 of Issue 156",
            "not_new": "A4 alone is not an independent observable",
            "new_split": "dA4/dp = j4_birth1 - j4_exit2, with the first line retained at exit",
            "required_archive": "K1,K2,ell1,site1,site2; iota is fixed to one",
            "tiny_oracle_result": (
                "ell is not a function of K1,K2, but is exactly independent of "
                "them conditional on a line birth in all three controls; their "
                "supported lines share one chi4 value"
            ),
            "direct_birth_boundary": (
                "0->2 jumps have no ell; the full unmarked K1 histogram is the "
                "sum over projective ell plus a typed DIRECT_RANK2 atom"
            ),
        },
        "all_pass": all_pass,
    }


def render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# Projective essential-birth Phase-A certificate",
        "",
        f"Overall exact gates: **{'PASS' if payload['all_pass'] else 'FAIL'}**.",
        "",
        "| geometry | N | paths | direct 0->2 | ell adds direction beyond (K1,K2) |",
        "|---|---:|---:|---:|---|",
    ]
    for row in payload["geometries"]:
        lines.append(
            f"| {row['id']} | {row['N']} | {row['all_paths']} | "
            f"{row['direct_rank2_paths_without_projective_line']} | "
            f"{row['ell_not_determined_by_K1_K2']} |"
        )
    lines.extend(
        [
            "",
        "The primitive line is canonical up to sign and covariant under the tested "
            "unimodular basis changes.  On the axis quotient, D4 rotations leave "
            "`chi4` fixed and reflections conjugate it.  Complement/Alexander duality "
            "maps `(K1,K2,ell,site1,site2)` to "
            "`(N+1-K2,N+1-K1,ell,site2,site1)` exactly.",
            "",
            "Integral saturation fixes `iota=1`: after the first line-bearing birth the "
            "subgroup is `Z ell`, and after the second birth it is `Z^2`.",
            "",
            "The exact controls sharpen the word *independent*.  The line is not "
            "determined by `(K1,K2)`, yet conditional on a line-bearing birth it "
            "factorizes exactly from `(K1,K2)` in all three controls.  Moreover the "
            "two supported lines are related by a quarter-turn and have the same "
            "`chi4`, so these tiny quotients show a projective direction degree of "
            "freedom but no additional spin-4 value.  A larger quotient is needed to "
            "test directional spin-4 bias.",
            "",
            "There is also an exact boundary to the proposed mark: axis-L2 and the "
            "C4 control contain direct `0->2` births.  Such paths have no canonical "
            "projective line.  Consequently the full `K1` histogram is recovered by "
            "summing over projective lines **plus** the typed `DIRECT_RANK2` atom; "
            "summing over `ell` alone recovers only line-bearing births.",
            "",
            "## Crosswalk to Issue 156",
            "",
            "Configuration by configuration, the plateau character",
            "",
            "```text",
            "A4(p)=E[1{tau1<=p<tau2} chi4(ell1)]",
            "```",
            "",
            "is exactly the fixed-p rank-one primitive-sector character.  It is not a "
            "new observable.  The genuinely new information is the source/sink split",
            "",
            "```text",
            "dA4/dp = j4_birth1(p) - j4_exit2(p),",
            "```",
            "",
            "where the exit at the second birth retains the line born at the first.  "
            "The certificate verifies this identity coefficient by coefficient in the "
            "degree-`N-1` Bernstein basis.  Saving only `A4` repeats Issue 156; saving "
            "both births localizes timing versus direction.",
        ]
    )
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = build_certificate()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload) + "\n", encoding="utf-8")
    return 0 if payload["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
