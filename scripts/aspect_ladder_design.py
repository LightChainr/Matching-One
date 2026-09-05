#!/usr/bin/env python3
"""Design the aspect-ratio ladder r = 1, 2, 4 at a single site count.

Why this exists.  The N=290 fingerprint measured one number, the spin-4
amplitude ratio between a square torus and a 2:1 rectangular torus, and got
1.880 +/- 0.177 against a predicted 11/4.  That excluded the weight-4 shape, and
plain area scaling, and no dependence at all -- but it is a single point, and a
single point cannot distinguish a law from a coincidence.  The one reading it
happens to sit on, the bare aspect ratio r itself, has no standing: it was
noticed after the number came back.

This module freezes the design that gives it standing or takes it away.  Three
things make the ladder work at one site count:

* **r = 4 is where the competitors separate.**  At r = 2 the weight-4 shape
  predicts 2.75 and the linear law predicts 2.00, which the existing 9 percent
  measurement cannot cleanly split.  At r = 4 they predict 10.99 and 4.00.
* **N = 1300 carries all three moduli.**  A torus of modulus r*i with N sites
  needs a Gaussian integer of norm N/r, so the ladder needs 1300, 650 and 325 to
  each be a sum of two squares in at least three ways.  325 is the smallest
  integer with three essentially distinct representations, and 650 = (1+i)*325
  and 1300 = 2*325 inherit them.  Nothing smaller does this.
* **Three orientations fit the spin-8 leakage out.**  Two orientations determine
  a constant and a spin-4 amplitude, leaving spin 8 as an unremovable systematic
  -- which is exactly the caveat the N=290 run had to carry.  Three determine
  C, A4 and A8 together.

The r = 2 rung is therefore also a replication: it re-measures the surprising
N=290 number at a different site count, with the systematic fitted rather than
bounded.

Nothing here is a measurement.  This is the design and the competitor list,
frozen before any block is run.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Optional, Sequence

from mpmath import mp

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "aspect-ladder-design" / "latest.json"
SCHEMA = "matching-one/aspect-ladder-design/v1"

SITE_COUNT = 1300
ASPECTS = (1, 2, 4)
PRECISION = 40


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def representations(norm: int) -> list[tuple[int, int]]:
    """Essentially distinct a+bi with a^2+b^2 = norm and 0 <= a <= b.

    Sign and order changes are lattice symmetries: they rotate the torus by a
    multiple of 90 degrees, which leaves cos(4*theta) alone.
    """
    found = []
    a = 0
    while 2 * a * a <= norm:
        remainder = norm - a * a
        b = int(round(remainder ** 0.5))
        if b * b == remainder and b >= a:
            found.append((a, b))
        a += 1
    return found


def cos4(a: int, b: int) -> Fraction:
    norm = a * a + b * b
    return Fraction(a**4 - 6 * a * a * b * b + b**4, norm * norm)


def cos8(a: int, b: int) -> Fraction:
    value = cos4(a, b)
    return 2 * value * value - 1


def period_matrix(a: int, b: int, aspect: int) -> list[int]:
    """Row-major <w, aspect*i*w>, the engine's integer period matrix."""
    return [a, -aspect * b, b, aspect * a]


def design_matrix(members: Sequence[tuple[int, int]]) -> list[list[Fraction]]:
    return [[Fraction(1), cos4(a, b), cos8(a, b)] for a, b in members]


def invert3(matrix: Sequence[Sequence[Fraction]]) -> tuple[list[list[Fraction]], Fraction]:
    (a, b, c), (d, e, f), (g, h, i) = matrix
    determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
    _require(determinant != 0, "the three orientations do not determine C, A4 and A8")
    adjugate = [
        [(e * i - f * h), -(b * i - c * h), (b * f - c * e)],
        [-(d * i - f * g), (a * i - c * g), -(a * f - c * d)],
        [(d * h - e * g), -(a * h - b * g), (a * e - b * d)],
    ]
    return [[value / determinant for value in row] for row in adjugate], determinant


def family(aspect: int, site_count: int = SITE_COUNT) -> dict[str, Any]:
    _require(site_count % aspect == 0, f"aspect {aspect} does not divide {site_count}")
    norm = site_count // aspect
    members = representations(norm)
    _require(len(members) >= 3, f"norm {norm} has only {len(members)} orientations")
    members = members[:3]
    matrix = design_matrix(members)
    inverse, determinant = invert3(matrix)
    return {
        "aspect_ratio": aspect,
        "modulus": "i" if aspect == 1 else f"{aspect}i",
        "gaussian_norm": norm,
        "members": [
            {
                "gaussian_integer": f"{a}+{b}i",
                "period_matrix_row_major": period_matrix(a, b, aspect),
                "sites": (a * a + b * b) * aspect,
                "cos4theta": _text(cos4(a, b)),
                "cos8theta": _text(cos8(a, b)),
            }
            for a, b in members
        ],
        "design_determinant": _text(determinant),
        "variance_amplification": {
            "A4": _text(sum(value * value for value in inverse[1])),
            "A8": _text(sum(value * value for value in inverse[2])),
            "meaning": (
                "Var(A) = this factor times the per-orientation variance, for "
                "equal counts on the three orientations"
            ),
        },
    }


def shape_predictions(aspects: Sequence[int] = ASPECTS) -> dict[str, Any]:
    """Area-normalized weight-k amplitude ratios A(r*i)/A(i), from the shape module."""
    import sys

    sys.path.insert(0, str(ROOT / "scripts"))
    from modulus_shape_discrimination import (  # noqa: E402
        SHAPES,
        VANISHING_TOLERANCE,
        lattice_amplitude,
    )

    rows = []
    with mp.workdps(PRECISION):
        for label, weight, shape, reason in SHAPES:
            base = lattice_amplitude(mp.mpc(1, 0), mp.mpc(0, 1), weight, shape)
            if abs(base) <= VANISHING_TOLERANCE:
                rows.append({
                    "shape": label,
                    "modular_weight": weight,
                    "on_the_list_because": reason,
                    "usable": False,
                    "why_not": "the shape vanishes at the square point, so the ratio is undefined",
                })
                continue
            ratios = {}
            for aspect in aspects:
                value = lattice_amplitude(mp.mpc(1, 0), mp.mpc(0, aspect), weight, shape)
                ratios[str(aspect)] = mp.nstr(mp.re(value / base), 12)
            rows.append({
                "shape": label,
                "modular_weight": weight,
                "on_the_list_because": reason,
                "usable": True,
                "amplitude_ratio_by_aspect": ratios,
            })
    return {"rows": rows}


def non_modular_competitors(aspects: Sequence[int] = ASPECTS) -> list[dict[str, Any]]:
    """Laws that are not modular shapes at all, including the post-hoc one.

    ``bare_aspect_ratio`` is here because the N=290 point sits on it, and it was
    noticed *after* that point came back.  Naming it in a frozen design before
    the ladder is run is the whole reason this file exists: it converts a
    post-hoc reading into a prospective competitor, or kills it.
    """
    return [
        {
            "law": "bare_aspect_ratio",
            "formula": "A(r)/A(1) = r",
            "predicted_by_aspect": {str(r): str(float(r)) for r in aspects},
            "standing": (
                "post-hoc on the N=290 point (z = -0.68 there), prospective from "
                "this design forward"
            ),
        },
        {
            "law": "plain_area_scaling",
            "formula": "A(r)/A(1) = r^2",
            "predicted_by_aspect": {str(r): str(float(r * r)) for r in aspects},
            "standing": "prospective; excluded at 12 sigma by the N=290 run",
        },
        {
            "law": "no_modulus_dependence",
            "formula": "A(r)/A(1) = 1",
            "predicted_by_aspect": {str(r): "1.0" for r in aspects},
            "standing": "prospective; excluded at 5.0 sigma by the N=290 run",
        },
    ]


def build_result(site_count: int = SITE_COUNT) -> dict[str, Any]:
    families = [family(aspect, site_count) for aspect in ASPECTS]
    determinants = {row["design_determinant"] for row in families}
    amplifications = {row["variance_amplification"]["A4"] for row in families}
    return {
        "schema": SCHEMA,
        "status": "design_only_no_measurement",
        "site_count": site_count,
        "roadmap_item": 2,
        "families": families,
        "families_share_one_design_matrix": len(determinants) == 1,
        "families_share_one_variance_amplification": len(amplifications) == 1,
        "why_they_share_it": (
            "1300 = 2*650 = 4*325, and 1300 is divisible by 4, so both parts of every "
            "representation of it are even: the r=1 family is exactly 2 times the r=4 "
            "family and samples the identical three orientations. The r=2 family is "
            "(1+i) times the r=4 family, a 45 degree turn, under which cos8 is "
            "invariant and cos4 changes sign; the sign change is undone by the "
            "reversal of the sorted order, so all three design matrices have the same "
            "determinant, not merely the same magnitude. No family is the noisy one"
        ),
        "estimator": (
            "within each family solve O(theta) = C + A4 cos4theta + A8 cos8theta on "
            "its three members; the score is the vector (A4(2i)/A4(i), A4(4i)/A4(i))"
        ),
        "score_is_amplitude_free": (
            "both entries are ratios of A4 within one observable, so the unknown "
            "non-universal overlap constant and the additive constant C cancel"
        ),
        "modular_shape_predictions": shape_predictions(),
        "non_modular_competitors": non_modular_competitors(),
        "what_separates_at_r4": (
            "at r=2 the weight-4 shape predicts 2.75 and the bare aspect ratio 2.00, "
            "which the existing 9 percent measurement cannot split; at r=4 they "
            "predict 10.99 and 4.00, a factor of 2.7 apart, so about 16 percent "
            "relative precision on the r=4 amplitude ratio decides it either way"
        ),
        "what_the_r2_rung_adds": (
            "a replication of the N=290 number at a different site count with the "
            "spin-8 leakage fitted rather than bounded. The N=290 design could only "
            "bound it, because two orientations determine only C and A4"
        ),
        "cost_relative_to_the_n290_design": {
            "orientations_per_family": 3,
            "variance_amplification_ratio": "0.8933 / 0.5438 = 1.64",
            "runs_per_family_ratio": "3 / 2 = 1.5",
            "samples_for_the_same_sigma_on_A4": "about 2.5 times the N=290 design, per family, per site",
            "sites_per_sample_ratio": "1300 / 290 = 4.48",
        },
        "not_established_by_this_design": [
            "anything at all: no block has been run",
            "identification of the Q4 Jordan module, which the ratio cannot decide "
            "on its own (docs/astra/Q2-additive-shape-ambiguity.md)",
            "that the fitted A4 is the log slope rather than a leading amplitude",
        ],
        "before_running": [
            "size the pilot honestly. A 200000-sample pilot of the N=290 channel "
            "returned an amplitude five times the 20M value, and a sample count "
            "projected from it was wrong by about forty times",
            "freeze a prediction file naming every competitor above, including "
            "bare_aspect_ratio, and declare no optional stopping",
        ],
    }


def validate_result(result: Any, site_count: int = SITE_COUNT) -> dict[str, Any]:
    expected = build_result(site_count)
    if result != expected:
        raise ValueError("aspect ladder design does not exactly reproduce")
    return {"schema": SCHEMA, "status": "valid", "site_count": site_count}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--site-count", type=int, default=SITE_COUNT)
    arguments = parser.parse_args(argv)
    text = json.dumps(build_result(arguments.site_count), indent=2, sort_keys=True) + "\n"
    destination = arguments.output or DEFAULT_OUTPUT
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text, encoding="utf-8")
    print(f"wrote {destination.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
