#!/usr/bin/env python3
"""Rescore the committed N=580 aspect ladder without forming a ratio.

The frozen design scored one number, `A4(4i)/A4(i)`, and deliberately excluded
the r=2 rung from the decision because that rung carries a spin-8 leakage of the
opposite sign.  Two things follow from testing the whole response vector
instead, and neither is visible from the ratio.

**One.** The response is *concave* in the aspect ratio.  The second divided
difference on r = 1, 2, 4 is a linear functional of the response,

    f[1,2,4] = (m(4) - 3 m(2) + 2 m(1)) / 6,

so it needs no denominator and no covariance inversion to state.  Every frozen
competitor predicts it to be zero (the two linear-in-r families) or positive
(every modular family).  The measurement is negative.

**Two.** Reconciling any competitor with the r=2 rung requires a spin-8
amplitude comparable to or larger than the spin-4 amplitude it is supposed to
perturb.  The frozen design assumed the opposite -- `|A8/A4|` "well below 1" --
and that assumption traces to a model-selection result on a *homology*
character, `H4 0.4163/2` against `H8 16.0120/2`, not to any measurement of the
*angular* amplitude ratio.  The design's own `not_established` list already said
the `C + A4 cos4 + A8 cos8` form could not be checked at two orientations per
rung.  This is what it looks like when that unchecked assumption fails.

What this module does NOT do is claim which of those it is.  A large spin-8
amplitude and a wrong angular form are the same arithmetic here; separating them
needs three orientations per rung, which the design already costed at N=650.

The missing covariance.  The first scoring run stored ratios and pairwise
variance pieces, so `cov(r2, r4)` is not in the committed artifact.  Every number
below is therefore reported across the whole range of that entry which keeps the
covariance positive definite, and a conclusion that does not survive the range is
labelled as depending on it.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from projective_inference import ray_residual, subspace_residual


ROOT = Path(__file__).resolve().parents[1]
LADDER = ROOT / "results" / "aspect-ladder-n580" / "latest.json"
DEFAULT_OUTPUT = ROOT / "results" / "aspect-ladder-n580-projective" / "latest.json"
SCHEMA = "matching-one.aspect-ladder-n580-projective.v1"
CHANNEL = "P4_S_prime"
RUNGS = (1, 2, 4)

# |cos8 difference / cos4 difference| at this site count, identical in magnitude
# on all three rungs and opposite in sign on r=2.  From the frozen design.
LEAKAGE = 1148 / 21025

# The second divided difference on r = 1, 2, 4, as a linear functional.
CURVATURE_WEIGHTS = (2.0 / 6.0, -3.0 / 6.0, 1.0 / 6.0)


def load_response(path: Path = LADDER) -> dict[str, Any]:
    """The three-rung response vector and what is known of its covariance.

    Prefers the full matrix when the artifact carries one.  Falls back to the
    pairwise pieces, in which case ``cov(r2, r4)`` is reported as unknown rather
    than assumed to be zero -- assuming zero is a choice, and an undeclared one.
    """
    payload = json.loads(path.read_text(encoding="utf-8"))
    stored = payload.get("response", {}).get("channels", {}).get(CHANNEL)
    if stored is not None:
        return {
            "vector": list(stored["vector"]),
            "covariance": [list(row) for row in stored["covariance"]],
            "covariance_is_complete": True,
        }
    ratios = payload["ratios"]
    first = ratios["r4_over_r1"][CHANNEL]
    second = ratios["r2_over_r1"][CHANNEL]
    vector = [first["denominator_value"], second["numerator_value"], first["numerator_value"]]
    variances = [first["denominator_variance"], second["numerator_variance"],
                 first["numerator_variance"]]
    covariance = [[0.0] * 3 for _ in range(3)]
    for i in range(3):
        covariance[i][i] = variances[i]
    covariance[0][1] = covariance[1][0] = second["covariance"]
    covariance[0][2] = covariance[2][0] = first["covariance"]
    return {"vector": vector, "covariance": covariance, "covariance_is_complete": False}


def admissible_correlation_range(covariance: Sequence[Sequence[float]],
                                 steps: int = 400) -> tuple[float, float]:
    """The values of corr(r2, r4) that keep the 3x3 covariance positive definite."""
    scale = math.sqrt(covariance[1][1] * covariance[2][2])
    admissible = []
    for index in range(steps + 1):
        candidate = -1.0 + 2.0 * index / steps
        if _determinant(_with_corr(covariance, candidate, scale)) > 0.0:
            admissible.append(candidate)
    if not admissible:
        raise ValueError("no correlation keeps the covariance positive definite")
    return min(admissible), max(admissible)


def _with_corr(covariance, correlation: float, scale: float):
    filled = [list(row) for row in covariance]
    filled[1][2] = filled[2][1] = correlation * scale
    return filled


def _determinant(m) -> float:
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def curvature(vector: Sequence[float], covariance: Sequence[Sequence[float]]) -> dict[str, float]:
    """Second divided difference on r = 1, 2, 4, with its standard error.

    A linear functional, so this is exact -- no denominator, no matrix inverse,
    and nothing that degrades when one entry is close to zero.
    """
    value = sum(w * x for w, x in zip(CURVATURE_WEIGHTS, vector))
    variance = sum(
        CURVATURE_WEIGHTS[i] * covariance[i][j] * CURVATURE_WEIGHTS[j]
        for i in range(3) for j in range(3)
    )
    if variance <= 0.0:
        raise ValueError("non-positive variance for the curvature functional")
    error = math.sqrt(variance)
    return {"value": value, "standard_error": error, "z": value / error}


def required_spin8_ratio(vector: Sequence[float], covariance: Sequence[Sequence[float]],
                         ray: Sequence[float]) -> float:
    """How large |A8/A4| must be for this model to reach the r=2 rung.

    The frozen design records the leakage per rung with its sign: r=1 and r=4
    both carry -1148/21025 and r=2 carries +1148/21025.  So if the ratio
    rho = A8/A4 is the same at all three aspect ratios -- which is exactly what
    the design's own "bounded well below 1" presupposes when it is applied to
    every rung at once -- the measured quantities are

        m(1) = a v1 (1 - lambda rho)
        m(2) = a v2 (1 + lambda rho)
        m(4) = a v4 (1 - lambda rho)

    The amplitude is fitted from the two rungs that share a sign, which returns
    ``a (1 - lambda rho)`` rather than ``a``, and the middle rung then gives

        u = m(2) / (v2 * fitted)  =  (1 + lambda rho) / (1 - lambda rho),

    so ``rho = (u - 1) / (lambda (u + 1))`` exactly.  Dividing the raw gap by
    ``lambda`` instead -- the leading-order form this routine used until
    2026-09-06 -- omits the ``(u + 1)`` and overstates the requirement by that
    factor, which reaches 45 for the weight-12 rays.
    """
    pair_covariance = [[covariance[0][0], covariance[0][2]],
                       [covariance[2][0], covariance[2][2]]]
    fit = ray_residual([vector[0], vector[2]], pair_covariance, [ray[0], ray[2]])
    amplitude = fit["amplitudes"][0]
    predicted_middle = amplitude * ray[1]
    if predicted_middle == 0.0:
        return math.inf
    u = vector[1] / predicted_middle
    if u == -1.0:
        return math.inf
    return abs(u - 1.0) / (LEAKAGE * abs(u + 1.0))


def rescore(competitors: Mapping[str, Sequence[float]],
            response: Mapping[str, Any], steps: int = 40) -> dict[str, Any]:
    vector = response["vector"]
    covariance = response["covariance"]
    low, high = admissible_correlation_range(covariance)
    scale = math.sqrt(covariance[1][1] * covariance[2][2])
    grid = [low + (high - low) * i / steps for i in range(steps + 1)]

    rows: dict[str, Any] = {}
    for name, ray in competitors.items():
        centre = subspace_residual(vector, covariance, [list(ray)])
        sigmas = [
            subspace_residual(vector, _with_corr(covariance, c, scale), [list(ray)])
            ["equivalent_sigma"]
            for c in grid
        ]
        decisions = {value >= 3.0 for value in sigmas}
        rows[name] = {
            "ray": list(ray),
            "statistic": centre["statistic"],
            "degrees_of_freedom": centre["degrees_of_freedom"],
            "equivalent_sigma": centre["equivalent_sigma"],
            "sigma_over_admissible_correlations": [min(sigmas), max(sigmas)],
            "excluded_at_3_sigma": centre["equivalent_sigma"] >= 3.0,
            "verdict_survives_the_missing_covariance": len(decisions) == 1,
            "curvature_predicted": sum(
                w * x for w, x in zip(CURVATURE_WEIGHTS, ray)
            ),
            "required_abs_A8_over_A4_to_reach_r2": required_spin8_ratio(
                vector, covariance, ray
            ),
        }

    measured = curvature(vector, covariance)
    curvature_range = [curvature(vector, _with_corr(covariance, c, scale))["z"] for c in grid]
    return {
        "schema": SCHEMA,
        "source_artifact": str(LADDER.relative_to(ROOT)),
        "channel": CHANNEL,
        "rungs": list(RUNGS),
        "response_vector": list(vector),
        "covariance_is_complete": response["covariance_is_complete"],
        "admissible_corr_r2_r4": [low, high],
        "curvature": {
            "definition": "second divided difference on r = 1, 2, 4: (m4 - 3 m2 + 2 m1) / 6",
            "why_it_is_safe": (
                "a linear functional of the response, so it has no denominator and "
                "no matrix inverse; nothing about it degrades when one rung is "
                "close to zero, which is what broke the ratio test"
            ),
            **measured,
            "z_over_admissible_correlations": [min(curvature_range), max(curvature_range)],
            "sign_predicted_by_every_competitor": (
                "zero for the two families linear in r, strictly positive for every "
                "modular family; the measurement is negative"
            ),
        },
        "competitors": rows,
        "spin8_bound_provenance": {
            "assumed_by_the_frozen_design": "|A8/A4| well below 1",
            "how_the_requirement_is_solved": (
                "from the per-rung leakage signs the frozen design "
                "records -- r=1 and r=4 carry -1148/21025 and r=2 carries "
                "+1148/21025 -- so with rho = A8/A4 the same at all three rungs, "
                "u = m(2) / (v2 * amplitude fitted from r=1 and r=4) equals "
                "(1 + lambda rho) / (1 - lambda rho) and rho = (u-1)/(lambda(u+1))"
            ),
            "what_that_solution_assumes": (
                "that A8/A4 is the same at all three aspect ratios. That is what "
                "a single bound on |A8/A4| presupposes when it is applied to every "
                "rung at once, so it is the design's own assumption rather than a "
                "new one -- but it is an assumption, and a rung-dependent ratio "
                "would change these numbers"
            ),
            "traces_to": (
                "predictions/modulus_fingerprint_n290_v2_20260905.yaml, quoting the "
                "committed H4-beats-H8 results as H4 0.4163/2 against H8 16.0120/2"
            ),
            "what_that_actually_is": (
                "a model-selection result on a homology character, not a measurement "
                "of the angular spin-8 to spin-4 amplitude ratio. The inference from "
                "one to the other is a plausibility argument and is not quantified "
                "anywhere in the repository"
            ),
        },
        "what_this_does_not_separate": (
            "a spin-8 amplitude comparable to spin-4, and a C + A4 cos4 + A8 cos8 "
            "form that is simply wrong, produce the same arithmetic here. Two "
            "orientations per rung determine C and A4 exactly with nothing left "
            "over, so this run cannot tell them apart. Three orientations per rung "
            "can, and the design already costed that at N=650"
        ),
    }


def load_competitors(path: Path = ROOT / "predictions" / "aspect_ladder_n580_20260905.yaml"):
    """The frozen competitor list, as rays over the three rungs.

    Read from the frozen prediction rather than restated here, so that a
    competitor cannot be quietly added or dropped between the freeze and the
    rescore.
    """
    import yaml

    frozen = yaml.safe_load(path.read_text(encoding="utf-8"))["competing_predictions"]
    return {name: (1.0, float(pair[0]), float(pair[1])) for name, pair in frozen.items()}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ladder", type=Path, default=LADDER)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args(argv)
    payload = rescore(load_competitors(), load_response(arguments.ladder))
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(arguments.output)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
