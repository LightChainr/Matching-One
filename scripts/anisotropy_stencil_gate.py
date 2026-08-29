#!/usr/bin/env python3
"""Exact spin-4 stencil gate for an improved percolation action.

For a planar edge vector ``z = dx + i dy``, the fourth angular moment

    |z|^4 exp(4 i arg z)

is exactly ``z**4``.  Integer stencil vectors therefore admit a rational,
roundoff-free oracle.  The module tests the cheapest proposal in Issue 106:
can the exactly-critical inhomogeneous square-bond family, which contains
horizontal and vertical edges only, tune this microscopic proxy through zero?

This is a geometry/proxy calculation.  It does not identify the weights that
multiply the continuum spin-4 operator in the renormalized action.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Tuple


@dataclass(frozen=True)
class EdgeOrbit:
    """One directed-vector representative with an explicit multiplicity."""

    dx: int
    dy: int
    weight: Fraction
    multiplicity: int = 1


def spin4(dx: int, dy: int) -> Tuple[int, int]:
    """Return the exact real and imaginary parts of ``(dx+i*dy)**4``."""

    real = dx**4 - 6 * dx * dx * dy * dy + dy**4
    imag = 4 * dx * dy * (dx * dx - dy * dy)
    return real, imag


def weighted_moment(orbits: Iterable[EdgeOrbit]) -> Tuple[Fraction, Fraction]:
    """Sum the exact weighted fourth moment of a finite stencil."""

    real = Fraction(0)
    imag = Fraction(0)
    for orbit in orbits:
        if orbit.weight < 0:
            raise ValueError("stencil weights must be nonnegative")
        if orbit.multiplicity < 1:
            raise ValueError("orbit multiplicity must be positive")
        orbit_real, orbit_imag = spin4(orbit.dx, orbit.dy)
        factor = orbit.weight * orbit.multiplicity
        real += factor * orbit_real
        imag += factor * orbit_imag
    return real, imag


def axis_probability_moment(t: Fraction) -> Tuple[Fraction, Fraction]:
    """Use critical probabilities ``p_h=t, p_v=1-t`` as proxy weights."""

    if not 0 <= t <= 1:
        raise ValueError("t must lie in [0,1]")
    return weighted_moment(
        [
            EdgeOrbit(1, 0, t, multiplicity=2),
            EdgeOrbit(0, 1, 1 - t, multiplicity=2),
        ]
    )


def axis_variance_moment(t: Fraction) -> Tuple[Fraction, Fraction]:
    """Use Bernoulli variances as a second nonnegative proxy weighting."""

    if not 0 <= t <= 1:
        raise ValueError("t must lie in [0,1]")
    weight = t * (1 - t)
    return weighted_moment(
        [
            EdgeOrbit(1, 0, weight, multiplicity=2),
            EdgeOrbit(0, 1, weight, multiplicity=2),
        ]
    )


def axis_diagonal_moment(
    axis_weight: Fraction, diagonal_weight: Fraction
) -> Tuple[Fraction, Fraction]:
    """Fourth moment of C4-complete unit-axis and integer-diagonal shells."""

    return weighted_moment(
        [
            EdgeOrbit(1, 0, axis_weight, multiplicity=2),
            EdgeOrbit(0, 1, axis_weight, multiplicity=2),
            EdgeOrbit(1, 1, diagonal_weight, multiplicity=2),
            EdgeOrbit(1, -1, diagonal_weight, multiplicity=2),
        ]
    )


def reflected_orbit_moment(dx: int, dy: int, weight: Fraction) -> Tuple[Fraction, Fraction]:
    """Moment of the four reflection images of a vector off the axes."""

    if dx == 0 or dy == 0:
        raise ValueError("use an off-axis vector for a four-element orbit")
    return weighted_moment(
        [
            EdgeOrbit(dx, dy, weight),
            EdgeOrbit(dx, -dy, weight),
            EdgeOrbit(-dx, dy, weight),
            EdgeOrbit(-dx, -dy, weight),
        ]
    )


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return "%d/%d" % (value.numerator, value.denominator)


def complex_record(moment: Tuple[Fraction, Fraction]) -> dict[str, str]:
    return {"real": fraction_text(moment[0]), "imag": fraction_text(moment[1])}


def build_artifact() -> dict[str, Any]:
    samples = [Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1)]
    probability_samples = {
        fraction_text(t): complex_record(axis_probability_moment(t)) for t in samples
    }
    variance_samples = {
        fraction_text(t): complex_record(axis_variance_moment(t)) for t in samples
    }

    axis_unit = spin4(1, 0)
    vertical_unit = spin4(0, 1)
    diagonal_unit = spin4(1, 1)
    antidiagonal_unit = spin4(1, -1)
    cancellation = axis_diagonal_moment(Fraction(4), Fraction(1))

    assert axis_unit == (1, 0)
    assert vertical_unit == (1, 0)
    assert diagonal_unit == (-4, 0)
    assert antidiagonal_unit == (-4, 0)
    assert all(axis_probability_moment(t) == (2, 0) for t in samples)
    assert cancellation == (0, 0)

    return {
        "schema": "matching-one/anisotropy-stencil-gate/v1",
        "issue": 106,
        "status": "exact_phase0_stencil_no_go",
        "definition": {
            "moment": "A4=sum_e w_e*(dx_e+i*dy_e)^4",
            "identity": "(dx+i*dy)^4=|v|^4*exp(4i*theta)",
            "arithmetic": "exact integer/rational",
            "weight_domain": "nonnegative local proxy weights",
        },
        "critical_axis_family": {
            "surface": "p_h+p_v=1",
            "parameterization": "p_h=t, p_v=1-t, 0<=t<=1",
            "probability_weighted_formula": "A4=2*t+2*(1-t)=2",
            "probability_weighted_samples": probability_samples,
            "variance_weighted_formula": "A4=4*t*(1-t)",
            "variance_weighted_samples": variance_samples,
            "general_nonnegative_formula": "A4=2*(w_h+w_v)",
            "no_go": (
                "Every nonzero horizontal/vertical stencil with nonnegative weights "
                "has strictly positive real A4; no interior critical parameter tunes it to zero."
            ),
        },
        "orbit_certificate": {
            "unit_horizontal": complex_record(
                (Fraction(axis_unit[0]), Fraction(axis_unit[1]))
            ),
            "unit_vertical": complex_record(
                (Fraction(vertical_unit[0]), Fraction(vertical_unit[1]))
            ),
            "integer_diagonal": complex_record(
                (Fraction(diagonal_unit[0]), Fraction(diagonal_unit[1]))
            ),
            "integer_antidiagonal": complex_record(
                (Fraction(antidiagonal_unit[0]), Fraction(antidiagonal_unit[1]))
            ),
        },
        "minimal_escape": {
            "necessary_sign_condition": (
                "With reflection-cancelled imaginary part and nonnegative weights, "
                "cancelling a retained positive axis contribution needs at least one "
                "orbit with cos(4*theta)<0."
            ),
            "axis_shell_formula": "A4_axis=4*w_axis",
            "diagonal_shell_formula": "A4_diagonal=-16*w_diagonal",
            "combined_formula": "A4=4*w_axis-16*w_diagonal",
            "exact_zero_ratio": "w_axis/w_diagonal=4",
            "certificate_at_axis4_diagonal1": complex_record(cancellation),
            "interpretation": (
                "An oblique/diagonal orbit is geometrically necessary for a positive-weight "
                "zero; this certificate does not establish an exactly-critical mixed-shell model."
            ),
        },
        "decision": {
            "reject": "blind zero search inside the axes-only inhomogeneous square-bond family",
            "retain": (
                "search for an exactly-critical star-triangle/isoradial family with a negative-phase "
                "edge orbit, then derive its physical coupling weights before simulation"
            ),
        },
        "claim_boundary": {
            "proved": (
                "the exact fourth-moment signs, the axes-only nonnegative-weight no-go, "
                "and the axis/diagonal cancellation ratio for the declared integer stencil"
            ),
            "not_proved": (
                "that any chosen proxy weight equals the renormalized spin-4 coupling, "
                "that the mixed axis/diagonal stencil is exactly critical, or that its measured H4 amplitude vanishes"
            ),
        },
        "sources": [
            {
                "id": "arXiv:1105.5535",
                "use": "critical surface p_h+p_v=1 for inhomogeneous square-bond percolation",
            },
            {
                "id": "arXiv:1204.0505",
                "use": "canonical critical isoradial families and star-triangle transport",
            },
        ],
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    critical = artifact["critical_axis_family"]
    escape = artifact["minimal_escape"]
    lines = [
        "# Exact spin-4 stencil feasibility gate",
        "",
        "## Decision",
        "",
        "Reject a blind zero search inside the horizontal/vertical-only exactly-critical square-bond family.",
        "For every nonzero nonnegative axis weighting, `A4` is strictly positive.",
        "",
        "## Exact calculations",
        "",
        "| stencil/proxy | exact `A4` | zero? |",
        "|---|---:|---|",
        "| axis shell, arbitrary nonnegative weights | `2(w_h+w_v)` | only the empty stencil |",
        "| critical axis family, probability weights | `2` | no |",
        "| critical axis family, Bernoulli-variance weights | `4t(1-t)` | only degenerate endpoints |",
        "| C4 axis + integer-diagonal shells | `4w_axis-16w_diagonal` | `w_axis/w_diagonal=4` |",
        "",
        "The unit axis orbit has spin-4 phase `+1`; the integer diagonal has moment `-4` because",
        "its squared length is two and its spin-4 phase is `-1`.",
        "",
        "## Geometric gate",
        "",
        escape["necessary_sign_condition"],
        "Thus an improved-action candidate with nonnegative microscopic weights must contain an oblique",
        "orbit in a negative spin-4 phase sector. The exact axis/diagonal cancellation certificate is",
        "`w_axis:w_diagonal = 4:1` for the declared integer-step normalization.",
        "",
        "## Evidence boundary",
        "",
        "This is an exact statement about the declared microscopic fourth-moment proxy, not a derivation",
        "of the renormalized spin-4 coupling. The calculation neither constructs an exactly-critical",
        "mixed axis/diagonal model nor predicts a vanishing measured H4 amplitude. The next admissible",
        "step is to find an exactly-critical star-triangle/isoradial family with a negative-phase orbit",
        "and derive the physical edge weights before simulation.",
        "",
        "## Reproduction",
        "",
        "```bash",
        "python scripts/anisotropy_stencil_gate.py --format json",
        "python scripts/anisotropy_stencil_gate.py --format markdown",
        "python -m unittest tests.test_anisotropy_stencil_gate",
        "```",
        "",
        "Sources: Grimmett--Manolescu `arXiv:1105.5535` and `arXiv:1204.0505`.",
        "",
    ]
    assert critical["probability_weighted_formula"] == "A4=2*t+2*(1-t)=2"
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    artifact = build_artifact()
    if args.format == "json":
        rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    else:
        rendered = render_markdown(artifact)

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
