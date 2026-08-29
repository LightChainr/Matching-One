#!/usr/bin/env python3
"""Exact complement-parity channels of the homology rank-birth insertion."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import comb
from pathlib import Path
from typing import Any, Optional, Sequence

from digital_alexander_filtration_oracle import rank_mark
from homology_rank_birth_insertion import (
    _add_polynomial,
    _bernstein_term,
    _differentiate,
    _evaluate,
    _fraction_text,
    _mask_active,
    _polynomial_text,
    _rank_polynomial,
    _rank_state_cache,
    _scale_polynomial,
    geometry_specs,
)
from integer_period_torus import IntegerTorusGeometry, Vector
from marked_pivotal_h4_reference import landing_mark


Polynomial = list[Fraction]
RankState = tuple[int, Optional[Vector], Optional[int]]


def spin4_character(vector: Vector) -> tuple[Fraction, Fraction]:
    """Return Re/Im of (x+iy)^4/|x+iy|^4 exactly."""

    x, y = vector
    radius_squared = x * x + y * y
    if radius_squared == 0:
        raise ValueError("zero vector has no angular character")
    denominator = radius_squared * radius_squared
    return (
        Fraction(x**4 - 6 * x * x * y * y + y**4, denominator),
        Fraction(4 * x * y * (x * x - y * y), denominator),
    )


def _character_text(character: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {"cos4": _fraction_text(character[0]), "sin4": _fraction_text(character[1])}


def _rank_states(
    geometry: IntegerTorusGeometry,
    *,
    matching: bool,
) -> list[RankState]:
    if not matching:
        return _rank_state_cache(geometry)
    return [
        rank_mark(geometry, _mask_active(mask, geometry.n), matching=True)
        for mask in range(1 << geometry.n)
    ]


def _gate_record(before: RankState, after: RankState) -> dict[str, Any]:
    rank_before, line_before, index_before = before
    rank_after, line_after, index_after = after
    if not 0 <= rank_before <= rank_after <= 2:
        raise AssertionError("rank insertion is not monotone")
    gate_01 = int(rank_before == 0 and rank_after >= 1)
    gate_12 = int(rank_before <= 1 and rank_after == 2)
    if rank_after - rank_before != gate_01 + gate_12:
        raise AssertionError("rank increment does not equal the two crossed gates")

    line: Optional[Vector] = None
    index: Optional[int] = None
    if (rank_before, rank_after) == (0, 1):
        line, index = line_after, index_after
    elif (rank_before, rank_after) == (1, 2):
        line, index = line_before, index_before
    return {
        "rank_before": rank_before,
        "rank_after": rank_after,
        "gate_01": gate_01,
        "gate_12": gate_12,
        "even": gate_01 + gate_12,
        "odd": gate_12 - gate_01,
        "ell": line,
        "iota": index,
    }


def complement_pair(
    geometry: IntegerTorusGeometry,
    black_without: Sequence[bool],
    vertex: int,
) -> tuple[dict[str, Any], dict[str, Any], list[bool]]:
    """Build primal insertion and its complement-reversed matching insertion."""

    black_zero = list(black_without)
    black_zero[vertex] = False
    black_one = list(black_zero)
    black_one[vertex] = True
    primal = _gate_record(
        rank_mark(geometry, black_zero, matching=False),
        rank_mark(geometry, black_one, matching=False),
    )

    # white v=0 is the complement of black v=1; white v=1 is the
    # complement of black v=0.  This is insertion reversal at the root.
    white_zero = [not value for value in black_one]
    white_zero[vertex] = False
    white_one = list(white_zero)
    white_one[vertex] = True
    matching = _gate_record(
        rank_mark(geometry, white_zero, matching=True),
        rank_mark(geometry, white_one, matching=True),
    )
    return primal, matching, white_zero


def _substitute_one_minus_p(polynomial: Sequence[Fraction]) -> Polynomial:
    result = [Fraction(0)] * len(polynomial)
    for power, coefficient in enumerate(polynomial):
        for degree in range(power + 1):
            result[degree] += coefficient * comb(power, degree) * ((-1) ** degree)
    return result


def _polynomial_from_counts(counts: Sequence[Fraction]) -> Polynomial:
    total = len(counts)
    result = [Fraction(0)] * total
    for occupied, count in enumerate(counts):
        result = _add_polynomial(
            result,
            [count * value for value in _bernstein_term(occupied, total - 1)],
        )
    return result


def _rank_one_probability_polynomial(
    geometry: IntegerTorusGeometry,
    states: Sequence[RankState],
) -> Polynomial:
    result = [Fraction(0)] * (geometry.n + 1)
    for mask, state in enumerate(states):
        if state[0] == 1:
            result = _add_polynomial(
                result,
                _bernstein_term(mask.bit_count(), geometry.n),
            )
    return result


def _empty_channel_arrays(n: int) -> dict[str, dict[str, list[Fraction]]]:
    return {
        side: {
            name: [Fraction(0)] * n
            for name in (
                "even",
                "odd",
                "line_cos4_even",
                "line_cos4_odd",
                "line_sin4_even",
                "line_sin4_odd",
                "local_h4_even",
                "local_h4_odd",
                "simultaneous_0_to_2",
            )
        }
        for side in ("primal", "matching")
    }


def _accumulate_channels(
    arrays: dict[str, list[Fraction]],
    occupied: int,
    record: dict[str, Any],
    multiplicity: int,
    *,
    physical_line: Optional[Vector],
    local_h4: Optional[int],
) -> None:
    even = record["even"]
    odd = record["odd"]
    arrays["even"][occupied] += multiplicity * even
    arrays["odd"][occupied] += multiplicity * odd
    if record["rank_before"] == 0 and record["rank_after"] == 2:
        arrays["simultaneous_0_to_2"][occupied] += multiplicity
    if physical_line is not None:
        cos4, sin4 = spin4_character(physical_line)
        arrays["line_cos4_even"][occupied] += multiplicity * even * cos4
        arrays["line_cos4_odd"][occupied] += multiplicity * odd * cos4
        arrays["line_sin4_even"][occupied] += multiplicity * even * sin4
        arrays["line_sin4_odd"][occupied] += multiplicity * odd * sin4
    if local_h4 is not None:
        arrays["local_h4_even"][occupied] += multiplicity * even * local_h4
        arrays["local_h4_odd"][occupied] += multiplicity * odd * local_h4


def _geometry_summary(
    name: str,
    geometry: IntegerTorusGeometry,
    *,
    roots: Sequence[int],
    root_multiplicity: int,
    local_radius: Optional[int],
) -> dict[str, Any]:
    primal_states = _rank_states(geometry, matching=False)
    matching_states = _rank_states(geometry, matching=True)
    arrays = _empty_channel_arrays(geometry.n)
    failures: Counter[str] = Counter()
    line_index_pairs: Counter[str] = Counter()
    paired_transition_counts: Counter[str] = Counter()

    for vertex in roots:
        bit = 1 << vertex
        for mask in range(1 << geometry.n):
            if mask & bit:
                continue
            black_zero = _mask_active(mask, geometry.n)
            black_one_mask = mask | bit
            primal = _gate_record(primal_states[mask], primal_states[black_one_mask])

            white_zero_mask = 0
            for other in range(geometry.n):
                if other != vertex and not (mask & (1 << other)):
                    white_zero_mask |= 1 << other
            white_one_mask = white_zero_mask | bit
            matching = _gate_record(
                matching_states[white_zero_mask], matching_states[white_one_mask]
            )
            white_zero = _mask_active(white_zero_mask, geometry.n)

            if (
                matching["rank_before"],
                matching["rank_after"],
            ) != (2 - primal["rank_after"], 2 - primal["rank_before"]):
                failures["rank_complement"] += root_multiplicity
            if (primal["gate_01"], primal["gate_12"]) != (
                matching["gate_12"],
                matching["gate_01"],
            ):
                failures["gate_exchange"] += root_multiplicity
            if primal["ell"] != matching["ell"]:
                failures["ell_preservation"] += root_multiplicity
            if primal["iota"] != matching["iota"]:
                failures["iota_preservation"] += root_multiplicity

            paired_transition_counts[
                f"{primal['rank_before']}->{primal['rank_after']} | "
                f"{matching['rank_before']}->{matching['rank_after']}"
            ] += root_multiplicity
            if primal["ell"] is not None:
                line_index_pairs[
                    f"ell={primal['ell']};iota_B={primal['iota']};iota_W={matching['iota']}"
                ] += root_multiplicity

            primal_local: Optional[int] = None
            matching_local: Optional[int] = None
            if local_radius is not None:
                primal_mark = landing_mark(
                    geometry, black_zero, local_radius, open_matching=False
                )
                matching_mark = landing_mark(
                    geometry, white_zero, local_radius, open_matching=True
                )
                if primal_mark != matching_mark:
                    failures["local_mark_preservation"] += root_multiplicity
                primal_local = primal_mark["h4"]
                matching_local = matching_mark["h4"]

            primal_physical = (
                geometry.periods.period_vector(primal["ell"])
                if primal["ell"] is not None
                else None
            )
            matching_physical = (
                geometry.periods.period_vector(matching["ell"])
                if matching["ell"] is not None
                else None
            )
            occupied_black = mask.bit_count()
            occupied_white = white_zero_mask.bit_count()
            _accumulate_channels(
                arrays["primal"],
                occupied_black,
                primal,
                root_multiplicity,
                physical_line=primal_physical,
                local_h4=primal_local,
            )
            _accumulate_channels(
                arrays["matching"],
                occupied_white,
                matching,
                root_multiplicity,
                physical_line=matching_physical,
                local_h4=matching_local,
            )

    polynomials = {
        side: {name: _polynomial_from_counts(counts) for name, counts in channels.items()}
        for side, channels in arrays.items()
    }
    parity_failures: dict[str, bool] = {}
    for channel_name in arrays["primal"]:
        expected_sign = (
            -1 if channel_name.endswith("_odd") or channel_name == "odd" else 1
        )
        transformed = _substitute_one_minus_p(
            polynomials["matching"][channel_name]
        )
        expected = [expected_sign * value for value in transformed]
        parity_failures[channel_name] = (
            polynomials["primal"][channel_name] != expected
        )
    if any(parity_failures.values()):
        raise AssertionError(f"coefficientwise complement parity failed: {parity_failures}")

    matching_derivative = _differentiate(_rank_polynomial(geometry, primal_states))
    if polynomials["primal"]["even"] != matching_derivative:
        raise AssertionError("even insertion channel is not M prime")
    rank_one_derivative = _differentiate(
        _rank_one_probability_polynomial(geometry, primal_states)
    )
    if polynomials["primal"]["odd"] != [-value for value in rank_one_derivative]:
        raise AssertionError("odd insertion channel is not -P1 prime")

    half = Fraction(1, 2)
    half_values = {
        name: _fraction_text(_evaluate(polynomial, half))
        for name, polynomial in polynomials["primal"].items()
    }
    return {
        "id": name,
        "N": geometry.n,
        "period_matrix": [list(row) for row in geometry.periods.matrix],
        "root_sampling": {
            "roots": list(roots),
            "multiplicity": root_multiplicity,
        },
        "paired_transition_counts": dict(sorted(paired_transition_counts.items())),
        "complement_failure_counts": dict(failures),
        "coefficient_parity_failures": parity_failures,
        "line_index_pair_counts": dict(sorted(line_index_pairs.items())),
        "bernstein_counts": {
            side: {
                channel: [_fraction_text(value) for value in values]
                for channel, values in channels.items()
            }
            for side, channels in arrays.items()
        },
        "power_basis_coefficients_ascending": {
            side: {
                channel: _polynomial_text(values)
                for channel, values in channels.items()
            }
            for side, channels in polynomials.items()
        },
        "p_equals_half_primal": half_values,
        "exact_identifications": {
            "even": "f_01+f_12=M_prime",
            "odd": "f_12-f_01=-partial_p Prob(rank=1)=partial_p(P0+P2)",
        },
    }


def rotation_character_certificate() -> dict[str, Any]:
    probes = ((1, 0), (1, 1), (2, 1), (3, -2))
    records = []
    for vector in probes:
        x, y = vector
        character = spin4_character(vector)
        sign_reversed = spin4_character((-x, -y))
        quarter_turn = spin4_character((-y, x))
        gaussian_pi_over_4 = spin4_character((x - y, x + y))
        if sign_reversed != character or quarter_turn != character:
            raise AssertionError("projective/C4 character invariance failed")
        if gaussian_pi_over_4 != (-character[0], -character[1]):
            raise AssertionError("pi/4 spin-four sign failed")
        records.append(
            {
                "vector": list(vector),
                "chi4": _character_text(character),
                "chi4_minus_vector": _character_text(sign_reversed),
                "chi4_R_pi_over_2": _character_text(quarter_turn),
                "chi4_times_1_plus_i": _character_text(gaussian_pi_over_4),
            }
        )
    return {
        "definition": "chi4(ell)=(x+i y)^4/(x^2+y^2)^2",
        "projective_sign": "chi4(-ell)=chi4(ell)",
        "C4": "chi4(R_pi/2 ell)=chi4(ell)",
        "gaussian_pi_over_4": "chi4((1+i)ell)=-chi4(ell)",
        "records": records,
    }


def build_artifact() -> dict[str, Any]:
    geometries = [_geometry_summary(**spec) for spec in geometry_specs()]
    return {
        "schema": "matching-one/rank-birth-parity-channels/v1",
        "issue": 215,
        "base_commit": "3881e88",
        "status": "tiny_exact_complement_parity_oracle",
        "exact_algebra": {
            "complement_reversal": "rW(white v=0)=2-rB(black v=1); rW(white v=1)=2-rB(black v=0)",
            "gate_exchange": "I_01^B=I_12^W and I_12^B=I_01^W",
            "even_channel": "S=I_01+I_12=Delta r; S_G(p)=S_Ghat(1-p)",
            "odd_channel": "D=I_12-I_01; D_G(p)=-D_Ghat(1-p)",
            "integrated": "sum E[S]=M_prime; sum E[D]=-partial_p P(rank=1)",
            "simultaneous_birth": "0_to_2 contributes 2 to S and 0 to D",
        },
        "line_refinement": {
            "complement": "ell is preserved while the birth gate is exchanged",
            "channels": "chi4(ell)S remains complement-even; chi4(ell)D remains complement-odd",
            "simultaneous_boundary": "0_to_2 has no ell and therefore forms an unpolarized even sector; it cancels from D",
        },
        "rotation_character": rotation_character_certificate(),
        "geometries": geometries,
        "mechanism_inference": {
            "exact_mass": "S and its H4 marks are the canonical matching-even local derivative channel carrying M_prime",
            "new_odd_field": (
                "D_H4 is a canonical matching-odd local field measuring the derivative of rank-one "
                "persistence rather than total matching susceptibility"
            ),
            "conditional_hypothesis": (
                "If the global anomalous matching-odd H4 response is carried by the essential-H1 lifetime mode, "
                "the line- or landing-marked D channel should couple more directly than the old untyped pivotal mark"
            ),
        },
        "claim_boundary": {
            "exact": "gate exchange and channel parity generally follow from digital rank complementarity; declared coefficient oracles are exhaustive",
            "tiny_only": "integral iota equality and local landing-mark preservation are only exhaustively certified on the declared controls",
            "not_claimed": "nonzero continuum overlap, an asymptotic exponent, or a canonical ell for simultaneous 0_to_2",
        },
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    lines = [
        "# Rank-birth complement-parity channels",
        "",
        "Complement plus insertion reversal exchanges `I01 <-> I12`. Therefore",
        "",
        "```text",
        "S = I01+I12 = Delta r       (even; carries M')",
        "D = I12-I01 = -d P1/dp      (odd; rank-one-lifetime derivative)",
        "```",
        "",
        "A direct `0->2` contributes two to S and zero to D. For non-simultaneous births,",
        "the primitive plateau line ell is preserved by complement. Its exact character",
        "`chi4=(x+iy)^4/|x+iy|^4` refines S and D without changing their complement parity.",
        "Under multiplication by `1+i`, `chi4` changes sign.",
        "",
        "| geometry | S(1/2) | D(1/2) | line-cos4 S | line-cos4 D | local-H4 S | local-H4 D |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in artifact["geometries"]:
        half = row["p_equals_half_primal"]
        lines.append(
            f"| {row['id']} | {half['even']} | {half['odd']} | "
            f"{half['line_cos4_even']} | {half['line_cos4_odd']} | "
            f"{half['local_h4_even']} | {half['local_h4_odd']} |"
        )
    lines.extend(
        [
            "",
            "Every complement, gate, line and coefficient-parity residual is zero on the",
            "declared exhaustive controls. The odd H4 channel is a conditional mechanism target,",
            "not yet an operator identification.",
            "",
        ]
    )
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    artifact = build_artifact()
    rendered = (
        json.dumps(artifact, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(artifact)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
