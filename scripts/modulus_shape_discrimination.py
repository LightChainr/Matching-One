#!/usr/bin/env python3
"""What a rectangular-torus measurement of the spin-4 shape can and cannot decide.

The Q4-Jordan bridge predicts that the root-normalized spin-4 amplitude carries
the area-normalized weight-4 shape `g2(tau) = Im(tau)^2 E4(tau)`.  Before buying
production at a new torus shape it is worth knowing what the alternatives
predict for the same measurement, because a prediction that every alternative
also makes is not worth measuring.

Two things come out of the arithmetic and both change the design.

1.  Pairing each torus with the same torus at 45 degrees to the lattice — the
    axis/diagonal pair the orientation programme already uses — gives a
    difference channel that annihilates every modular weight congruent to 0
    mod 8 exactly.  All scalar corrections live there.  Weight 4 survives with
    a sign flip; so does weight 12, separated from weight 4 by a factor of ~12
    at aspect ratio 2.

2.  The one serious survivor is a spin-4 amplitude with plain area scaling and
    no modular correction.  Against that competitor the whole discriminating
    content is the single constant 1/E4(i) = 0.6869..., visible as 11/4 instead
    of 4 at aspect ratio 2.  E4(r i) -> 1 so fast that aspect ratios 3 and 4
    reproduce the naive r^2 law to 1e-5 and better.  Longer rectangles buy
    nothing against the competitor that matters.

The design conclusion is therefore that one rectangle is the whole experiment,
and its value is set by how well the ratio can be measured, not by how many
shapes are run.

Nothing here is a claim about the lattice.  It is arithmetic about candidate
shapes, and about which of them a measurement could tell apart.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp

from derive_rectangular_thermal_q4_hecke import exact_ratios

SCHEMA = "matching-one.modulus-shape-discrimination.v2"

# Eisenstein normalisation E_k = 1 - (2k/B_k) sum sigma_{k-1}(n) q^n.
EISENSTEIN_COEFFICIENT = {
    2: Fraction(-24),
    4: Fraction(240),
    6: Fraction(-504),
    8: Fraction(480),
    10: Fraction(-264),
    12: Fraction(65520, 691),
}

ASPECT_RATIOS = (2, 3, 4)


def divisor_power_sum(n: int, power: int) -> int:
    total = 0
    divisor = 1
    while divisor * divisor <= n:
        if n % divisor == 0:
            partner = n // divisor
            total += divisor**power
            if partner != divisor:
                total += partner**power
        divisor += 1
    return total


def eisenstein(tau: mp.mpc, weight: int, *, terms: int = 250) -> mp.mpc:
    coefficient = EISENSTEIN_COEFFICIENT[weight]
    scale = mp.mpf(coefficient.numerator) / mp.mpf(coefficient.denominator)
    q = mp.exp(2 * mp.pi * mp.j * tau)
    value = mp.mpc(1)
    for n in range(1, terms):
        value += scale * divisor_power_sum(n, weight - 1) * q**n
    return value


def discriminant(tau: mp.mpc, *, terms: int = 250) -> mp.mpc:
    q = mp.exp(2 * mp.pi * mp.j * tau)
    product = mp.mpc(1)
    for n in range(1, terms):
        product *= 1 - q**n
    return q * product**24


def lattice_amplitude(omega1: mp.mpc, omega2: mp.mpc, weight: int, shape) -> mp.mpc:
    """Area-normalized weight-k lattice amplitude for the lattice <omega1, omega2>.

    For a lattice L = omega1 Z + omega2 Z the weight-k Eisenstein sum scales as
    omega1^-k f(tau), and multiplying by Area(L)^(k/2) = |omega1|^k Im(tau)^(k/2)
    leaves (|omega1|/omega1)^k Im(tau)^(k/2) f(tau).  The prefactor is exp(-i k
    theta) with theta the lattice orientation: that is where the spin lives.
    """
    tau = omega2 / omega1
    if mp.im(tau) < 0:
        tau = -tau
    return (abs(omega1) / omega1) ** weight * mp.im(tau) ** (
        mp.mpf(weight) / 2
    ) * shape(tau)


def axis_and_diagonal(aspect: int, weight: int, shape) -> tuple[mp.mpc, mp.mpc]:
    """The same 1:aspect torus, realized axis-aligned and at 45 degrees."""
    axis = lattice_amplitude(mp.mpc(1, 0), mp.mpc(0, aspect), weight, shape)
    turn = mp.mpc(1, 1)
    diagonal = lattice_amplitude(turn, turn * mp.mpc(0, aspect), weight, shape)
    return axis, diagonal


SHAPES = (
    ("E4", 4, lambda t: eisenstein(t, 4), "the Q4-Jordan prediction"),
    ("E8", 8, lambda t: eisenstein(t, 8), "the next even weight; equals E4^2"),
    ("E12", 12, lambda t: eisenstein(t, 12), "weight-12 Eisenstein"),
    ("E4^3", 12, lambda t: eisenstein(t, 4) ** 3, "a weight-12 power of the same shape"),
    ("Delta", 12, discriminant, "the weight-12 cusp form"),
    ("E6", 6, lambda t: eisenstein(t, 6), "an odd-spin weight, for contrast"),
)

VANISHING_TOLERANCE = mp.mpf(10) ** -25


def production_design() -> dict:
    """The concrete N=290 realization of the axis/diagonal rectangular test.

    A square torus needs a Gaussian integer of norm N; a 1:2 rectangular torus of
    the same site count needs one of norm N/2.  Both families need at least two
    representations to separate the spin-4 amplitude from the scalar part.  N=290
    is the smallest size in the project's existing production range where both
    hold and the angular leverage is near its maximum -- and it is a size this
    repository has already run.
    """

    def cos4(a: int, b: int) -> Fraction:
        norm = a * a + b * b
        return Fraction(a**4 - 6 * a * a * b * b + b**4, norm * norm)

    square = [(17, 1), (13, 11)]
    rectangular = [(12, 1), (9, 8)]

    def family(reps, multiplier):
        rows = []
        for a, b in reps:
            value = cos4(a, b)
            rows.append(
                {
                    "gaussian_integer": f"{a}+{b}i",
                    "period_vectors": [[a, b], [-multiplier * b, multiplier * a]],
                    "sites": (a * a + b * b) * multiplier,
                    "cos4theta": f"{value.numerator}/{value.denominator}",
                }
            )
        leverage = cos4(*reps[0]) - cos4(*reps[1])
        return rows, leverage

    def cos8(a: int, b: int) -> Fraction:
        value = cos4(a, b)
        return 2 * value * value - 1

    def spin8_leakage(reps) -> Fraction:
        """What a spin-8 component contributes to the spin-4 projector.

        An exact 45-degree turn would cancel weight 8, but two Gaussian integers
        of the *same* norm can never differ by exactly 45 degrees -- multiplying
        by 1+i doubles the norm.  So the realized pair removes spin 0 exactly
        (any orientation pair does) and leaves spin 8 at this coefficient.
        """
        return (cos8(*reps[0]) - cos8(*reps[1])) / (cos4(*reps[0]) - cos4(*reps[1]))

    square_rows, square_leverage = family(square, 1)
    rectangular_rows, rectangular_leverage = family(rectangular, 2)
    square_leak = spin8_leakage(square)
    rectangular_leak = spin8_leakage(rectangular)

    return {
        "site_count": 290,
        "square_family": {
            "modulus": "i",
            "lattice": "<w, i w> with |w|^2 = 290",
            "members": square_rows,
            "angular_leverage": f"{square_leverage.numerator}/{square_leverage.denominator}",
        },
        "rectangular_family": {
            "modulus": "2i",
            "lattice": "<w, 2i w> with |w|^2 = 145",
            "members": rectangular_rows,
            "angular_leverage": f"{rectangular_leverage.numerator}/{rectangular_leverage.denominator}",
        },
        "leverages_are_equal": square_leverage == rectangular_leverage,
        "spin8_leakage": {
            "square": f"{square_leak.numerator}/{square_leak.denominator}",
            "rectangular": f"{rectangular_leak.numerator}/{rectangular_leak.denominator}",
            "equal_and_opposite": square_leak == -rectangular_leak,
            "meaning": (
                "a spin-8 component of amplitude A8 enters the square estimator "
                "as +0.0546 A8 and the rectangular one as -0.0546 A8, so it biases "
                "the ratio by roughly -0.055 (A8/A4 + A8'/A4'). The committed "
                "H4-beats-H8 results bound A8/A4 well below 1, but this is a "
                "systematic on the score and not a statistical error"
            ),
            "how_to_remove_it": (
                "three orientations per family determine C, A4 and A8 together. "
                "N=650 is the smallest size where both the square family "
                "(|w|^2=650) and the rectangular family (|w|^2=325) have three "
                "representations"
            ),
        },
        "why_they_are_equal": (
            "multiplication by 1+i maps the norm-145 representations to the "
            "norm-290 ones, which is the 45-degree turn; the two families "
            "therefore sample the same pair of cos4theta values with the roles "
            "exchanged, so the ratio estimator carries no relative variance "
            "penalty between its numerator and denominator"
        ),
        "estimator": (
            "within each family solve O(theta) = C + A cos4theta on its two "
            "members, so A = (O1-O2)/(cos4theta1-cos4theta2); the score is "
            "A_rectangular / A_square, predicted 11/4 against 4 for area scaling"
        ),
        "caveats": [
            "two members per family determine C and A exactly, leaving nothing "
            "over to check the cos4theta form itself at this size; that form "
            "rests on the existing orientation programme at other sizes",
            "the engine couples two period matrices per run, so the two families "
            "are two runs; treating them as independent for the ratio is "
            "conservative but loses the shared-field variance reduction",
            "the exact weight-8 cancellation belongs to an idealized 45-degree "
            "pair, which no two Gaussian integers of equal norm can realize; the "
            "realized pair cancels spin 0 exactly and leaves the spin-8 leakage "
            "recorded above",
        ],
    }


def render(dps: int = 40) -> dict:
    exact_two = exact_ratios()["E4hat_2i_over_E4hat_i"]

    with mp.workdps(dps):
        channels = []
        for label, weight, shape, reason in SHAPES:
            axis_one, diagonal_one = axis_and_diagonal(1, weight, shape)
            reference = axis_one - diagonal_one
            entry = {
                "shape": label,
                "modular_weight": weight,
                "on_the_list_because": reason,
                "diagonal_over_axis_at_square": mp.nstr(
                    mp.re(diagonal_one / axis_one), 12
                )
                if abs(axis_one) > VANISHING_TOLERANCE
                else None,
            }
            if abs(reference) <= VANISHING_TOLERANCE:
                entry["survives_difference_channel"] = False
                entry["why"] = (
                    "the 45-degree lattice gives the same amplitude, so the "
                    "axis-minus-diagonal difference is identically zero"
                    if abs(axis_one) > VANISHING_TOLERANCE
                    else "the shape vanishes at the square point"
                )
            else:
                entry["survives_difference_channel"] = True
                entry["difference_ratio"] = {
                    str(r): mp.nstr(
                        mp.re(
                            (lambda pair: pair[0] - pair[1])(
                                axis_and_diagonal(r, weight, shape)
                            )
                            / reference
                        ),
                        12,
                    )
                    for r in ASPECT_RATIOS
                }
            channels.append(entry)

        # The competitor that the difference channel does not remove: a spin-4
        # amplitude that simply scales with torus area.
        e4_at_square = mp.re(eisenstein(mp.mpc(0, 1), 4))
        naive = {str(r): str(r * r) for r in ASPECT_RATIOS}
        survivors = [row for row in channels if row["survives_difference_channel"]]
        q4 = next(row for row in survivors if row["shape"] == "E4")

        separations = []
        for r in ASPECT_RATIOS:
            predicted = mp.mpf(q4["difference_ratio"][str(r)])
            area = mp.mpf(r * r)
            separations.append(
                {
                    "aspect_ratio": r,
                    "q4_jordan": mp.nstr(predicted, 12),
                    "area_scaling": str(r * r),
                    "relative_distance": mp.nstr(abs(predicted - area) / predicted, 8),
                    "q4_over_area": mp.nstr(predicted / area, 12),
                }
            )

        # How much does a second aspect ratio add against that competitor?
        anchor = mp.mpf(separations[0]["q4_over_area"])
        marginal = [
            {
                "aspect_ratio": row["aspect_ratio"],
                "q4_over_area": row["q4_over_area"],
                "distance_from_the_r2_value": mp.nstr(
                    abs(mp.mpf(row["q4_over_area"]) - anchor) / anchor, 6
                ),
            }
            for row in separations[1:]
        ]

        nearest = min(
            (
                row
                for row in survivors
                if row["shape"] != "E4"
            ),
            key=lambda row: abs(
                mp.mpf(row["difference_ratio"]["2"]) - mp.mpf(q4["difference_ratio"]["2"])
            ),
        )

        return {
            "schema": SCHEMA,
            "claim_level": "C0_exact_design_input",
            "question": (
                "which candidate spin-4 shapes can a rectangular-torus measurement "
                "separate, and how many rectangles are worth running"
            ),
            "design": {
                "observable_pair": (
                    "the same 1:r torus realized twice, once axis-aligned and once "
                    "with the lattice at 45 degrees; period vectors (1,0),(0,r) and "
                    "(1,1),(-r,r)"
                ),
                "score": "D(r) = A_axis(r) - A_diagonal(r), reported as D(r)/D(1)",
                "why_the_difference": (
                    "the orientation prefactor is exp(-i k theta), so a 45-degree "
                    "turn multiplies a weight-k amplitude by exp(-i k pi/4); weights "
                    "congruent to 0 mod 8 are unchanged and cancel exactly in D"
                ),
                "null_channel": (
                    "A_axis + A_diagonal must be consistent with zero for a pure "
                    "weight-4 shape, and carries the cancelled scalar content"
                ),
            },
            "q4_jordan_prediction_at_2i": f"{exact_two.numerator}/{exact_two.denominator}",
            "production_design": production_design(),
            "E4_at_square_point": mp.nstr(e4_at_square, 12),
            "candidate_shapes": channels,
            "surviving_competitor_without_modular_structure": {
                "description": "spin-4 amplitude scaling with torus area and nothing else",
                "difference_ratio": naive,
            },
            "separation_from_area_scaling": separations,
            "marginal_value_of_extra_aspect_ratios": marginal,
            "nearest_modular_competitor_at_2i": {
                "shape": nearest["shape"],
                "difference_ratio_at_2i": nearest["difference_ratio"]["2"],
            },
            "design_conclusions": [
                "the axis/diagonal difference removes every weight congruent to 0 "
                "mod 8 exactly, which is where all scalar finite-size corrections "
                "live; this costs nothing and is not a statistical argument",
                "the nearest surviving modular competitor sits roughly a factor of "
                "12 away at aspect ratio 2, so modular weight is not the hard part",
                "the one competitor that survives is plain area scaling, and the "
                "entire discriminating content against it is the constant 1/E4(i); "
                "the measurement is 11/4 against 4, a gap of 45 percent of the "
                "predicted value",
                "aspect ratios 3 and 4 reproduce the area law to within 1e-5 and "
                "better, so they add essentially nothing against that competitor: "
                "one rectangle is the whole experiment",
                "3 sigma between 11/4 and 4 needs the ratio to about 15 percent "
                "relative, and that number, not the number of geometries, is what "
                "the production design has to hit",
            ],
            "not_established": [
                "that any lattice observable equals one of these shapes",
                "that the measured amplitude is the root-normalized log slope "
                "rather than a leading amplitude with the same symmetry",
                "that the difference channel is free of lattice-level spin-4 "
                "contamination that is not of modular origin",
                "that 11/4 is a prediction rather than a conditional statement: "
                "the additive shape A~(tau) is not fixed by the Jordan relation "
                "and the ratio is clean only if the normalization removes the "
                "same block, which is the open question in docs/astra/Q2",
                "any statement about achievable precision, which depends on the "
                "production design and not on this arithmetic",
            ],
            "effect_on_the_open_normalization_question": (
                "the axis/diagonal difference removes every weight congruent to "
                "0 mod 8, so any scalar part of the unfixed additive shape drops "
                "out of the score; a spin-4 part does not. The channel narrows "
                "docs/astra/Q2 without closing it, and the narrowing is exact "
                "rather than statistical."
            ),
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=40)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    text = json.dumps(render(arguments.dps), indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(text, end="")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
