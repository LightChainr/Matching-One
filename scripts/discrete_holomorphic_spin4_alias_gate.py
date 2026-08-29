#!/usr/bin/env python3
"""Exact character-alias gate for direction-only C4 defect readouts."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable


Direction = tuple[int, int]
GaussianRational = tuple[Fraction, Fraction]


AXIS: Direction = (1, 0)
DIAGONAL: Direction = (1, 1)


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def spin4_phase(direction: Direction) -> GaussianRational:
    """Return exp(-4 i theta) exactly for the integer direction x+i y."""
    x, y = direction
    norm = x * x + y * y
    if norm == 0:
        raise ValueError("direction must be nonzero")
    # (x-i y)^2 = a+i b, then square once more and divide by |x+i y|^4.
    a = x * x - y * y
    b = -2 * x * y
    denominator = norm * norm
    return Fraction(a * a - b * b, denominator), Fraction(2 * a * b, denominator)


def c4_orbit(direction: Direction) -> tuple[Direction, ...]:
    if direction == (0, 0):
        raise ValueError("direction must be nonzero")
    orbit = []
    current = direction
    for _ in range(4):
        orbit.append(current)
        current = (-current[1], current[0])
    return tuple(orbit)


def character_matrix(direction: Direction) -> tuple[tuple[GaussianRational, GaussianRational], ...]:
    """Rows are the scalar and spin-4 characters on one C4 orbit."""
    scalar = (Fraction(1), Fraction(0))
    return tuple((scalar, spin4_phase(item)) for item in c4_orbit(direction))


def real_character_rank(direction: Direction) -> int:
    phases = [phase for _, phase in character_matrix(direction)]
    if any(imaginary != 0 for _, imaginary in phases):
        raise ValueError("this exact rank helper expects a real spin-4 phase")
    return 1 if len({real for real, _ in phases}) == 1 else 2


def orbit_averages(direction: Direction, defects: Iterable[Fraction]) -> dict[str, GaussianRational | Fraction]:
    values = tuple(Fraction(value) for value in defects)
    if len(values) != 4:
        raise ValueError("one C4 orbit requires exactly four defect values")
    phases = [spin4_phase(item) for item in c4_orbit(direction)]
    scalar = sum(values, Fraction(0)) / 4
    real = sum((value * phase[0] for value, phase in zip(values, phases)), Fraction(0)) / 4
    imaginary = sum((value * phase[1] for value, phase in zip(values, phases)), Fraction(0)) / 4
    return {"scalar": scalar, "spin4": (real, imaginary)}


def serialize(value: Any) -> Any:
    if isinstance(value, Fraction):
        return fraction_text(value)
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    return value


def build_artifact() -> dict[str, Any]:
    axis_orbit = c4_orbit(AXIS)
    diagonal_orbit = c4_orbit(DIAGONAL)
    axis_phases = tuple(spin4_phase(item) for item in axis_orbit)
    diagonal_phases = tuple(spin4_phase(item) for item in diagonal_orbit)

    # Orbit-averaged response to scalar + spin-4 amplitudes:
    # [A_axis]     [1  1] [a0]
    # [A_diagonal] [1 -1] [a4].
    two_orbit_matrix = ((Fraction(1), Fraction(1)), (Fraction(1), Fraction(-1)))
    determinant = two_orbit_matrix[0][0] * two_orbit_matrix[1][1] - two_orbit_matrix[0][1] * two_orbit_matrix[1][0]

    axis_constant = orbit_averages(AXIS, (Fraction(7),) * 4)
    diagonal_constant = orbit_averages(DIAGONAL, (Fraction(7),) * 4)

    assert axis_phases == ((Fraction(1), Fraction(0)),) * 4
    assert diagonal_phases == ((Fraction(-1), Fraction(0)),) * 4
    assert axis_constant["spin4"] == (axis_constant["scalar"], Fraction(0))
    assert diagonal_constant["spin4"] == (-diagonal_constant["scalar"], Fraction(0))
    assert determinant == -2

    return serialize({
        "schema": "matching-one/discrete-holomorphic-spin4-alias-gate/v1",
        "issue": 109,
        "data_class": "exact direction-character algebra only",
        "single_orbit": {
            "axis_directions": axis_orbit,
            "axis_spin4_phases": axis_phases,
            "axis_character_rank": real_character_rank(AXIS),
            "diagonal_directions": diagonal_orbit,
            "diagonal_spin4_phases": diagonal_phases,
            "diagonal_character_rank": real_character_rank(DIAGONAL),
            "constant_defect": "7",
            "axis_scalar_average": axis_constant["scalar"],
            "axis_naive_spin4_average": axis_constant["spin4"],
            "diagonal_scalar_average": diagonal_constant["scalar"],
            "diagonal_naive_spin4_average": diagonal_constant["spin4"],
        },
        "two_orbit_separation": {
            "response_matrix": two_orbit_matrix,
            "determinant": determinant,
            "rank": 2,
            "scalar_projection": "(axis_average + diagonal_average)/2",
            "spin4_projection": "(axis_average - diagonal_average)/2",
        },
        "decision": (
            "A direction-only defect on one C4 orbit cannot distinguish scalar from spin 4. "
            "At least two C4 orbits with unequal exp(-4 i theta0) phases, or additional typed/internal "
            "edge information, are required."
        ),
        "boundary": (
            "This gate does not implement Zhou edge observables, prove discrete harmonicity, measure "
            "percolation, establish L^-2 scaling, or identify a matching/KdV operator."
        ),
    })


def render_markdown(artifact: dict[str, Any]) -> str:
    single = artifact["single_orbit"]
    pair = artifact["two_orbit_separation"]
    return "\n".join([
        "# Spin-4 direction-character alias gate", "",
        "This is exact direction-character algebra, not percolation data.", "",
        "| direction orbit | `exp(-4 i theta)` | scalar/spin-4 character rank |",
        "|---|---:|---:|",
        f"| axial | `{single['axis_spin4_phases'][0]}` | `{single['axis_character_rank']}` |",
        f"| diagonal | `{single['diagonal_spin4_phases'][0]}` | `{single['diagonal_character_rank']}` |", "",
        "A constant defect of `7` therefore has naive axial spin-4 average `[7, 0]`, exactly equal to",
        "its scalar average. On the diagonal orbit the same readout is only the negative scalar average.", "",
        "Using both orbit averages gives the exact response matrix", "",
        f"`{pair['response_matrix']}`", "",
        f"with determinant `{pair['determinant']}` and rank `{pair['rank']}`. Hence", "",
        f"- scalar: `{pair['scalar_projection']}`;",
        f"- spin 4: `{pair['spin4_projection']}`.", "",
        "## Decision", "", artifact["decision"], "", "## Boundary", "", artifact["boundary"], "",
    ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = build_artifact()
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n" if args.format == "json" else render_markdown(artifact)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
