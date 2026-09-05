#!/usr/bin/env python3
"""Design the aspect-ratio ladder r = 1, 2, 4 at a single site count.

Why this exists.  The N=290 fingerprint measured one number, the spin-4
amplitude ratio between a square torus and a 2:1 rectangular torus, and got
1.880 +/- 0.177 against a predicted 11/4.  That excluded the weight-4 shape, and
plain area scaling, and no dependence at all -- but it is a single point, and a
single point cannot distinguish a law from a coincidence.  The one reading it
happens to sit on, the bare aspect ratio r itself, has no standing: it was
noticed after the number came back.

This module chooses and freezes the design that gives that reading its one
chance to lose.

Three things fix the shape of the design:

* **r = 4 is where the competitors separate.**  At r = 2 the weight-4 shape
  predicts 2.75 and the linear law predicts 2.00, which the existing 9 percent
  measurement cannot cleanly split.  At r = 4 they predict 10.99 and 4.00.
* **One site count carries all three moduli.**  A torus of modulus r*i with N
  sites needs a Gaussian integer of norm N/r, so the ladder needs N, N/2 and N/4
  each to be a sum of two squares in at least two ways.  Since N = 4m makes all
  three inherit m's representations, the candidates are 4m for m a sum of two
  squares in two ways.
* **The choice among those candidates is not free.**  They differ by a factor of
  two in angular leverage and by a factor of thirty-six in spin-8 leakage, and
  the leakage signs across the three rungs decide whether the systematic cancels
  in the score or compounds in it.  ``rank_candidates`` searches them and states
  the objective rather than asserting a winner.

What made this concrete: a 1M pilot.  An earlier version of this design used
N=1300 and three orientations per family, to fit the spin-8 amplitude out
instead of bounding it.  The pilot killed it twice over.  The analysis path
returned zeros at N=1300 -- the binomial tail's recurrence underflows near 790
sites, now fixed in ``analyze_p48_retrospective`` -- and the measured
per-difference noise there was 2.0 times the N=290 noise against a signal some
6.5 times smaller, which puts a decisive run about three orders of magnitude
outside what we can spend.  N=580 needs no such fit: its r=1 and r=4 rungs carry
the *same* leakage, so the leading spin-8 bias cancels in the ratio that
discriminates.

Nothing here is a measurement.  This is the design and the competitor list,
frozen before any scoring block is run.
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
SCHEMA = "matching-one/aspect-ladder-design/v2"

SITE_COUNT = 580
ASPECTS = (1, 2, 4)
SEARCH_CEILING = 1000
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


def rung(aspect: int, site_count: int) -> Optional[dict[str, Any]]:
    """One family: a modulus, two orientations, a leverage and a leakage.

    Two orientations determine the constant and the spin-4 amplitude:

        (O_1 - O_2) / delta_cos4  =  A4  +  A8 * (delta_cos8 / delta_cos4).

    The bracketed factor is the leakage.  It does not shrink with samples, so it
    is a property of the geometry and belongs in the design rather than in the
    error bar.
    """
    if site_count % aspect:
        return None
    members = representations(site_count // aspect)
    if len(members) < 2:
        return None
    first, second = members[0], members[1]
    leverage = cos4(*first) - cos4(*second)
    if leverage == 0:
        return None
    leakage = (cos8(*first) - cos8(*second)) / leverage
    return {
        "aspect_ratio": aspect,
        "modulus": "i" if aspect == 1 else f"{aspect}i",
        "gaussian_norm": site_count // aspect,
        "first_rep": list(first),
        "second_rep": list(second),
        "first_matrix_row_major": period_matrix(*first, aspect),
        "second_matrix_row_major": period_matrix(*second, aspect),
        "delta_cos4": _text(leverage),
        "delta_cos8": _text(cos8(*first) - cos8(*second)),
        "spin8_leakage": _text(leakage),
        "orientations_available": len(members),
    }


def ladder(site_count: int) -> Optional[list[dict[str, Any]]]:
    rungs = [rung(aspect, site_count) for aspect in ASPECTS]
    if any(row is None for row in rungs):
        return None
    return rungs


def rank_candidates(ceiling: int = SEARCH_CEILING) -> list[dict[str, Any]]:
    """Every site count up to the ceiling that carries the ladder, best first.

    The objective, in order: **maximum angular leverage**, because it divides
    into every amplitude and so sets the sample budget; then **minimum spin-8
    leakage**, because that one does not shrink with samples at all.
    """
    scored = []
    for site_count in range(4, ceiling + 1, 4):
        rungs = ladder(site_count)
        if rungs is None:
            continue
        leverages = {row["delta_cos4"] for row in rungs}
        leakages = {abs(Fraction(row["spin8_leakage"])) for row in rungs}
        scored.append({
            "site_count": site_count,
            "leverage": _text(Fraction(next(iter(leverages)))) if len(leverages) == 1 else None,
            "leverage_is_shared": len(leverages) == 1,
            "max_abs_leakage": _text(max(leakages)),
            "leakage_cancels_in_the_r4_over_r1_ratio":
                Fraction(rungs[0]["spin8_leakage"]) == Fraction(rungs[2]["spin8_leakage"]),
        })
    scored.sort(key=lambda row: (
        -Fraction(row["leverage"]) if row["leverage_is_shared"] else Fraction(0),
        Fraction(row["max_abs_leakage"]),
        row["site_count"],
    ))
    return scored


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
    rungs = ladder(site_count)
    _require(rungs is not None, f"{site_count} does not carry the r = 1, 2, 4 ladder")
    leverages = {row["delta_cos4"] for row in rungs}
    ranked = rank_candidates()
    _require(ranked[0]["site_count"] == site_count,
             f"the search prefers {ranked[0]['site_count']} over the frozen {site_count}")
    return {
        "schema": SCHEMA,
        "status": "design_only_no_measurement",
        "site_count": site_count,
        "roadmap_item": 2,
        "rungs": rungs,
        "runs_total": len(rungs),
        "leverage_is_shared_across_rungs": len(leverages) == 1,
        "shared_leverage": next(iter(leverages)) if len(leverages) == 1 else None,
        "spin8_cancels_in_the_discriminating_ratio":
            Fraction(rungs[0]["spin8_leakage"]) == Fraction(rungs[2]["spin8_leakage"]),
        "why_that_matters": (
            "the estimator is A4 + A8 * leakage. When the r=1 and r=4 rungs carry "
            "the same leakage, the ratio A4(4i)/A4(i) -- the one that discriminates "
            "-- has that bias cancel to leading order instead of compounding, which "
            "is what the N=290 pair could not do: its two families had equal and "
            "*opposite* leakage, so the systematic entered the score twice"
        ),
        "candidate_search": {
            "ceiling": SEARCH_CEILING,
            "objective": (
                "maximum shared angular leverage first, because it divides into "
                "every amplitude and sets the sample budget; then minimum spin-8 "
                "leakage, because that one does not shrink with samples"
            ),
            "best_five": ranked[:5],
            "candidates_found": len(ranked),
        },
        "estimator": (
            "per rung, one paired run gives (O_first - O_second) / delta_cos4 = "
            "A4 + A8 * leakage; the score is the vector "
            "(A4(2i)/A4(i), A4(4i)/A4(i))"
        ),
        "score_is_amplitude_free": (
            "both entries are ratios of A4 within one observable, so the unknown "
            "non-universal overlap constant and the additive constant cancel"
        ),
        "modular_shape_predictions": shape_predictions(),
        "non_modular_competitors": non_modular_competitors(),
        "what_separates_at_r4": (
            "at r=2 the weight-4 shape predicts 2.75 and the bare aspect ratio 2.00, "
            "which the existing 9 percent measurement cannot split; at r=4 they "
            "predict 10.99 and 4.00, a factor of 2.7 apart, so about 30 percent "
            "relative precision on the r=1 amplitude decides it either way"
        ),
        "what_the_r2_rung_adds": (
            "a replication of the N=290 number at a different site count, in a "
            "geometry where the leakage enters the r=2/r=1 ratio once rather than "
            "twice"
        ),
        "rejected_alternative": {
            "site_count": 1300,
            "shape": "three orientations per rung, spin-8 fitted rather than bounded",
            "why_rejected": (
                "measured, not estimated. A 1M pilot at N=1300 gave a per-difference "
                "noise of 0.0131 against 0.0065 at N=290, a factor of 2.0 for a 4.5 "
                "times larger torus, while the amplitude falls roughly as N^-5/4. "
                "Reaching a decisive ratio there is about three orders of magnitude "
                "beyond what we can spend"
            ),
            "what_the_pilot_also_found": (
                "the analysis path returned zeros at N=1300: the binomial tail's "
                "recurrence starts at (1-p)^N, which underflows to exactly zero near "
                "790 sites at the percolation threshold, and then stays zero. Fixed "
                "in analyze_p48_retrospective by anchoring the recurrence at the mode"
            ),
        },
        "not_established_by_this_design": [
            "anything at all: no scoring block has been run",
            "identification of the Q4 Jordan module, which the ratio cannot decide "
            "on its own (docs/astra/Q2-additive-shape-ambiguity.md)",
            "that the fitted A4 is the log slope rather than a leading amplitude",
            "freedom from spin-8 in the r=2/r=1 entry, where the leakage does not "
            "cancel; only the r=4/r=1 entry is clean",
        ],
        "before_running": [
            "size the budget from a pilot at this site count, not from an "
            "extrapolation. A 200000-sample pilot of the N=290 channel returned an "
            "amplitude five times its 20M value; standard errors do extrapolate as "
            "one over root n, central values at small counts do not",
            "declare no optional stopping, and place the replica offset past every "
            "committed stream so the pilot can never be pooled into the run",
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
