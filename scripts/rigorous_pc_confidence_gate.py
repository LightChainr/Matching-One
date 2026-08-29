#!/usr/bin/env python3
"""Phase-0 certification and cost gate for rigorous pc confidence intervals.

The statistical core is deliberately small: for N independent block-event
trials and M successes, test the null ``pi <= p0`` with the exact upper tail
of ``Bin(N,p0)``.  All certification comparisons use ``Fraction`` arithmetic.
Floating-point calculations appear only in the explicitly heuristic cost table.
"""

from __future__ import annotations

import argparse
import json
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from typing import Any, Tuple


LEGACY_P0 = Fraction(8639, 10000)
MODERN_P0 = Fraction(8457, 10000)
FAMILYWISE_ALPHA = Fraction(1, 1_000_000)
SIDES = 2
ATTEMPTS_PER_SIDE = 3
PER_RUN_ALPHA = FAMILYWISE_ALPHA / (SIDES * ATTEMPTS_PER_SIDE)


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return "%d/%d" % (value.numerator, value.denominator)


def decimal_text(value: Fraction, significant_digits: int = 18) -> str:
    """Render a rational in normalized scientific notation."""

    with localcontext() as context:
        context.prec = significant_digits + 8
        decimal = Decimal(value.numerator) / Decimal(value.denominator)
        return format(decimal, ".%dE" % (significant_digits - 1))


def binomial_tail(n: int, minimum_successes: int, probability: Fraction) -> Fraction:
    """Return P[Bin(n, probability) >= minimum_successes] exactly."""

    if n < 0:
        raise ValueError("n must be nonnegative")
    if not 0 <= minimum_successes <= n + 1:
        raise ValueError("minimum_successes must lie in [0,n+1]")
    if not 0 <= probability <= 1:
        raise ValueError("probability must lie in [0,1]")
    if minimum_successes == 0:
        return Fraction(1)
    if minimum_successes == n + 1:
        return Fraction(0)
    if probability == 0:
        return Fraction(0)
    if probability == 1:
        return Fraction(1)

    failure = 1 - probability
    success_count = n
    term = probability**n
    tail = term
    while success_count > minimum_successes:
        previous = (
            term
            * success_count
            * failure
            / ((n - success_count + 1) * probability)
        )
        success_count -= 1
        term = previous
        tail += term
    return tail


def minimal_successes(n: int, null_probability: Fraction, alpha: Fraction) -> Tuple[int, Fraction]:
    """Smallest k with P[Bin(n,p0)>=k] <= alpha, using exact recurrence."""

    if not 0 < alpha < 1:
        raise ValueError("alpha must lie in (0,1)")
    if not 0 < null_probability < 1:
        raise ValueError("null_probability must lie in (0,1)")

    success_count = n
    term = null_probability**n
    tail = term
    if tail > alpha:
        return n + 1, Fraction(0)
    while success_count > 0:
        previous = (
            term
            * success_count
            * (1 - null_probability)
            / ((n - success_count + 1) * null_probability)
        )
        candidate = tail + previous
        if candidate > alpha:
            return success_count, tail
        success_count -= 1
        term = previous
        tail = candidate
    return 0, Fraction(1)


def certification_row(n: int, p0: Fraction, alpha: Fraction) -> dict[str, Any]:
    cutoff, tail = minimal_successes(n, p0, alpha)
    previous_tail = binomial_tail(n, cutoff - 1, p0) if cutoff > 0 else Fraction(1)
    return {
        "trials": n,
        "null_event_probability": fraction_text(p0),
        "per_run_alpha": fraction_text(alpha),
        "minimum_successes": cutoff,
        "tail_at_cutoff": decimal_text(tail),
        "tail_one_success_lower": decimal_text(previous_tail),
        "exact_decision": "tail_at_cutoff<=alpha<tail_one_success_lower",
    }


def power_row(n: int, cutoff: int, true_probability: Fraction) -> dict[str, Any]:
    power = binomial_tail(n, cutoff, true_probability)
    return {
        "true_event_probability": fraction_text(true_probability),
        "power": decimal_text(power),
    }


def heuristic_cost_row(base_width: float, target_width: float) -> dict[str, Any]:
    """Heuristic four-arm scaling only; never used in a certification claim."""

    if not 0 < target_width <= base_width:
        raise ValueError("target width must be positive and no larger than the baseline")
    ratio = base_width / target_width
    linear_scale = ratio ** (4.0 / 3.0)
    area_work = ratio ** (8.0 / 3.0)
    return {
        "target_width": format(target_width, ".1e"),
        "width_improvement_factor": ratio,
        "linear_scale_multiplier": linear_scale,
        "area_work_multiplier": area_work,
    }


def build_artifact() -> dict[str, Any]:
    legacy = certification_row(400, LEGACY_P0, PER_RUN_ALPHA)
    modern = certification_row(400, MODERN_P0, PER_RUN_ALPHA)
    legacy_paper_tail = binomial_tail(400, 378, LEGACY_P0)

    assert legacy["minimum_successes"] == 378
    assert modern["minimum_successes"] == 373
    assert legacy_paper_tail <= PER_RUN_ALPHA
    assert binomial_tail(400, 377, LEGACY_P0) > PER_RUN_ALPHA
    assert binomial_tail(400, 372, MODERN_P0) > PER_RUN_ALPHA

    power_probabilities = [Fraction(90, 100), Fraction(92, 100), Fraction(94, 100), Fraction(95, 100)]
    sample_sizes = [100, 200, 400, 800]
    sample_size_table = []
    for n in sample_sizes:
        cutoff, tail = minimal_successes(n, MODERN_P0, PER_RUN_ALPHA)
        sample_size_table.append(
            {
                "trials": n,
                "minimum_successes": cutoff,
                "null_tail": decimal_text(tail),
                "power": [power_row(n, cutoff, probability) for probability in power_probabilities],
            }
        )

    base_width = 5e-4
    cost_widths = [5e-4, 2.5e-4, 1e-4, 5e-5, 1e-5]
    cost_table = [heuristic_cost_row(base_width, width) for width in cost_widths]

    return {
        "schema": "matching-one/rigorous-pc-confidence-gate/v1",
        "issue": 112,
        "status": "phase0_methodologically_interesting_but_weak",
        "block_event": {
            "geometry": "two adjacent s-by-s cells Su,Sv with union Re of size 2s-by-s",
            "event": (
                "each half has a unique largest open cluster and the two largest clusters "
                "belong to one open cluster in Re"
            ),
            "locality": "the event depends only on sites in Re",
            "renormalized_dependence": (
                "events for vertex-disjoint block bonds use disjoint rectangles, hence define "
                "a 1-independent bond model on Z2"
            ),
            "lifting": "an infinite open block path forces an infinite original open cluster",
        },
        "deterministic_constants": {
            "legacy_2007_p0": fraction_text(LEGACY_P0),
            "modern_2022_upper_bound_on_pmax_Z2": fraction_text(MODERN_P0),
            "modern_use": "reject H0: block-event probability<=0.8457; strict rejection certifies >0.8457",
        },
        "error_budget": {
            "familywise_alpha": fraction_text(FAMILYWISE_ALPHA),
            "sides": SIDES,
            "attempts_per_side": ATTEMPTS_PER_SIDE,
            "per_run_alpha": fraction_text(PER_RUN_ALPHA),
            "allocation": "Bonferroni over three upper-bound and three matching-lattice lower-bound runs",
        },
        "legacy_reproduction": {
            **legacy,
            "paper_cutoff": 378,
            "paper_tail_recomputed": decimal_text(legacy_paper_tail),
        },
        "modern_protocol": {
            **modern,
            "sample_size_and_power_table": sample_size_table,
        },
        "two_sided_square_site_interval": {
            "upper": (
                "at parameter p, reject block-event probability<=p0 on square-site Z2; "
                "then report pc(square)<=p"
            ),
            "lower": (
                "at parameter q, run the same upper-bound protocol on the square-site matching graph; "
                "then pc(square)>=1-q by site matching duality"
            ),
            "final_confidence": "at least 1-familywise_alpha under predeclared tests and IID genuine Bernoulli trials",
        },
        "heuristic_cost_model": {
            "certification_status": "not part of the rigorous guarantee",
            "baseline_width": "5e-4",
            "assumption": "thermal window scales as s^(-3/4), per-trial work as area s^2",
            "formulas": {
                "linear_scale_multiplier": "(baseline_width/target_width)^(4/3)",
                "area_work_multiplier": "(baseline_width/target_width)^(8/3)",
            },
            "table": cost_table,
        },
        "feasibility_decision": {
            "class": "methodologically_interesting_but_weak",
            "reason": (
                "the method applies directly and the statistical certificate is cheap, but the block event "
                "needs a new open-boundary evaluator and geometric scale dominates cost; a 5e-4 interval "
                "was already achieved in 2007"
            ),
            "production_gate": (
                "first implement and independently validate the exact block event; use exploratory samples "
                "only to freeze s and p; then draw fresh final IID samples under the declared error budget"
            ),
            "deterministic_enumeration": (
                "naive enumeration is exponential in O(s^2), while frontier transfer is exponential in s; "
                "neither is a plausible route at production scales"
            ),
        },
        "randomness_boundary": {
            "statistical_theorem": "conditional on independent trials from genuine Bernoulli site fields",
            "pseudorandom_generator": (
                "software reproducibility and stream-domain separation do not by themselves prove genuine randomness"
            ),
            "deterministic_claim": "no deterministic new pc bound is produced by this artifact",
        },
        "sources": [
            {
                "id": "arXiv:math/0702232",
                "role": "block event, 1-independent reduction, exact binomial test, and 2007 interval",
            },
            {
                "id": "arXiv:2206.12335",
                "role": "rigorous improved bound p_max(Z2)<=0.8457",
            },
        ],
    }


def render_markdown(artifact: dict[str, Any]) -> str:
    legacy = artifact["legacy_reproduction"]
    modern = artifact["modern_protocol"]
    lines = [
        "# Certified finite-size bound feasibility gate",
        "",
        "## Decision",
        "",
        "**Methodologically interesting but weak.** The Riordan--Walters reduction applies directly",
        "to square-site percolation, and its binomial certificate is inexpensive. The dominant cost is",
        "making the block event sufficiently likely at a parameter close to `pc`; the repository does",
        "not yet implement this open-boundary event, and a width `5e-4` interval already exists.",
        "",
        "## Certifiable theorem shape",
        "",
        "Partition the plane into `s x s` cells. A renormalized bond is open when both adjacent cells",
        "have unique largest open clusters and those clusters connect in their `2s x s` union.",
        "Nonincident block bonds depend on disjoint site sets, so the block model is 1-independent.",
        "If its edge probability is strictly above `0.8457`, the modern deterministic bound on",
        "`p_max(Z^2)` implies percolation of the block model and hence of the original model.",
        "",
        "For `N` independent trials and `M` successes, reject `H0: pi<=0.8457` using the exact tail",
        "`P[Bin(N,0.8457)>=M]`. With family-wise error `1e-6`, three attempts on each of the original",
        "and matching lattices give per-run alpha `1/6000000`.",
        "",
        "## Exact statistical thresholds",
        "",
        "| protocol | N | required successes | null tail |",
        "|---|---:|---:|---:|",
        "| 2007 constant `0.8639` | %d | %d | `%s` |" % (
            legacy["trials"], legacy["minimum_successes"], legacy["tail_at_cutoff"]
        ),
        "| 2022 constant `0.8457` | %d | %d | `%s` |" % (
            modern["trials"], modern["minimum_successes"], modern["tail_at_cutoff"]
        ),
        "",
        "The first row exactly reproduces the published `400/378` test. The improved deterministic",
        "constant lowers the cutoff to `373`, but 400 trials have useful power only when the true block",
        "event probability is already high (about the mid-0.9 range).",
        "",
        "## Two-sided interval",
        "",
        "A successful square-site run at `p` gives `pc<=p`. A successful run on the square matching",
        "site graph at `q` gives `pc(square)>=1-q`. These are high-confidence statistical bounds, not",
        "deterministic inequalities, and require predeclared tests plus independent genuine Bernoulli trials.",
        "",
        "## Heuristic cost warning",
        "",
        "Using the paper's non-rigorous thermal-window scaling `delta p ~ s^(-3/4)`, reducing interval",
        "width by a factor `r` multiplies linear scale by `r^(4/3)` and area work by `r^(8/3)`.",
        "This estimate is for planning only and is not part of the certificate.",
        "",
        "| target width | linear-scale multiplier | area-work multiplier |",
        "|---:|---:|---:|",
    ]
    for row in artifact["heuristic_cost_model"]["table"]:
        lines.append(
            "| `%s` | `%.3f` | `%.3f` |"
            % (row["target_width"], row["linear_scale_multiplier"], row["area_work_multiplier"])
        )
    lines.extend(
        [
            "",
            "## Production gate",
            "",
            "Implement and independently validate the exact open-boundary block event first. Exploratory",
            "samples may choose `s,p`, but final certification samples must be fresh. Stream separation and",
            "reproducibility are necessary software controls; they do not turn a pseudorandom generator into",
            "a mathematical source of genuine randomness.",
            "",
            "Sources: Riordan--Walters `arXiv:math/0702232`; Balister--Johnston--Savery--Scott",
            "`arXiv:2206.12335`.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
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


if __name__ == "__main__":
    main()
