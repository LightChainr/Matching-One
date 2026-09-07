#!/usr/bin/env python3
"""#595: which readout resolves the angular amplitude best per sample?

#589's option-2 attempt failed not because the noncyclic offset is large but
because the *denominator* of ``rho = delta / A4`` is barely resolved: across the
#205 block ``|A4_z|`` runs from 0.06 to 2.18, and the sample count needed to
bound the Smith contamination scales as ``1 / A4_z^2``.  The spread between the
five channels ``obs()`` happens to return is a factor of ~9000 in cost, and no
production in this repository ever selected for it -- #205, #583 and #581 all
chose their observable before this figure of merit existed.

This script scans a predeclared family of readouts over the same archived
histograms and reports one number per readout:

    leverage = |A4| / se(A4)

with ``A4`` the two-cyclic-row angular amplitude and ``se`` from the same
aligned delete-one jackknife the published scorer uses.

WHAT THIS MAY AND MAY NOT REPORT (GOVERNANCE 2D, 2E)
----------------------------------------------------
This is an optimization over observables scored against a finite sample, so it
will find leverage that is partly noise.  Therefore:

  * it reports ONLY ``se(A4)``, ``|A4|``, and the amplification relative to the
    historical baseline ``M`` at the frozen ``p_ref``;
  * it NEVER reports ``delta``, ``rho``, a z-score against the noncyclic row, or
    any contamination bound for a channel it selected -- that would score a
    number the selection already saw;
  * a channel it selects has to be declared prospectively and bought with new
    samples before it may carry any verdict.

The deliverable is a design recommendation with a measured selection-optimism
penalty, not a measurement.

DEMOTED BY #596 / #595's own review, AND THAT IS CORRECT
--------------------------------------------------------
This script maximizes ``|A4| / se(A4)``.  That is **not** the objective the
downstream decision needs, and it should be read as a screen rather than as the
answer.  The decision #583 waits on is not "is A4 large" but "can the angular
direction be separated from the Smith nuisance well enough to bound the induced
contamination", and a channel can improve the first while damaging the second.

Worse, forming ``rho = delta / A4`` at all during design is invalid exactly where
this block sits: ``A4`` is 0.06-2.18 sigma from zero, and a delta-method error on
a ratio with an unresolved denominator understates a confidence set that is in
truth unbounded.  The correct object is the joint ``(A4, delta)`` covariance and
a Fieller/projective confidence set, scored by the *new sample* needed for that
set to clear both nuisance poles and land inside the declared band.

That calculation is ``scripts/p205_projective_channel_design.py`` (#595 v2).
This file is retained because it is computed, tested and honest about its own
optimism, and because raw ``A4`` SNR is still the fourth-ranked criterion there
-- but nothing here may be quoted as a design recommendation on its own.

The optimism is measured, not asserted: the 100 production batches are split
into interleaved halves, the readout is selected on one half, and its leverage
is reported on the other.  If the amplification does not transport, the answer
is that the channel spread is noise -- which would also mean ``Sp``'s 4.5x
advantage over ``M`` is not real, and is worth knowing.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import mpmath as mp

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from analyze_matching_parity_derivatives_fast import combine, obs, remove  # noqa: E402
from score_p205_derivative_smith_loading import (  # noqa: E402
    GEOMETRY_ORDER,
    angular_functional,
    quadratic_form,
)
from score_p205_norm5_conjugate_coalescence import (  # noqa: E402
    SIZES,
    aligned_rows,
    jackknife_covariance,
    load_contract,
    load_pair,
    pair_spec,
    sha256,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "p205-channel-leverage" / "latest.json"
SCHEMA = "matching-one.p205-channel-leverage.v1"
ISSUE = 595

#: ``D`` is excluded throughout: ``obs`` defines it as exactly ``M/2``, so it is
#: linearly dependent on ``M`` and would make the combination covariance
#: singular while adding nothing.
BASIS_CHANNELS: Tuple[str, ...] = ("M", "S", "Sp", "Dp")

#: The historical choice, and the baseline every amplification is measured
#: against: the fixed-p scalar at the frozen #205 probability.
BASELINE_CHANNEL = "M"

#: Predeclared probability grid, written before any leverage was computed.
#: Centred on the frozen p_ref and symmetric, so the grid itself cannot encode a
#: preference for one side of the threshold.
P_REF = "0.59274605079"
P_HALF_WIDTH = "0.035"
P_GRID_POINTS = 21

#: Interleaved half-split.  Interleaved rather than contiguous so both halves
#: span the same production stretch and the same replica-counter interval; a
#: contiguous split would confound selection optimism with drift.
SPLIT_PARITY = 2

#: Directions of the combination covariance below this fraction of the largest
#: eigenvalue are dropped rather than inverted.  The Hotelling direction is a
#: ratio of an estimated signal to an estimated covariance; inverting a
#: direction the sample cannot resolve manufactures leverage out of arithmetic.
COMBINATION_FLOOR = mp.mpf("1e-8")


def probability_grid(
    half_width: str = P_HALF_WIDTH, points: int = P_GRID_POINTS
) -> List[mp.mpf]:
    """A grid centred on ``p_ref`` with an odd point count, so ``p_ref`` is on it.

    Centring is not cosmetic: the baseline every amplification is measured
    against is the historical readout at exactly ``p_ref``, so that point has to
    be a grid point rather than an interpolation.
    """

    if points % 2 == 0:
        raise ValueError("the grid must have an odd point count so p_ref lies on it")
    centre = mp.mpf(P_REF)
    half = mp.mpf(half_width)
    step = 2 * half / (points - 1)
    return [centre - half + step * i for i in range(points)]


def channel_table(
    rows: Mapping[str, Sequence[Any]], batches: Sequence[int], p: mp.mpf
) -> Dict[str, Any]:
    """Point estimate and delete-one replicates for every basis channel at ``p``.

    One ``obs`` call yields all channels, so the scan costs the same as a
    single-channel scan.  The jackknife is rebuilt inside the given batch subset
    so that a half-sample is scored exactly the way the full sample is.
    """

    combined = {name: combine([rows[name][b] for b in batches]) for name in GEOMETRY_ORDER}
    point = {name: obs(combined[name], p) for name in GEOMETRY_ORDER}
    replicates = []
    for b in batches:
        replicates.append({
            name: obs(remove(combined[name], rows[name][b]), p)
            for name in GEOMETRY_ORDER
        })
    return {"point": point, "replicates": replicates}


def leverage(table: Mapping[str, Any], n: int, channel: str) -> Dict[str, mp.mpf]:
    """``|A4| / se(A4)`` for one channel, from the aligned delete-one jackknife."""

    weights = angular_functional(n)
    point = [table["point"][name][channel] for name in GEOMETRY_ORDER]
    amplitude = mp.fsum(weights[i] * point[i] for i in range(3))
    deleted = [
        [replicate[name][channel] for name in GEOMETRY_ORDER]
        for replicate in table["replicates"]
    ]
    covariance = jackknife_covariance(deleted)
    variance = quadratic_form(weights, covariance, weights)
    if variance <= 0:
        return {"A4": amplitude, "standard_error": mp.mpf("inf"), "leverage": mp.mpf(0)}
    error = mp.sqrt(variance)
    return {"A4": amplitude, "standard_error": error, "leverage": abs(amplitude) / error}


def amplitude_covariance(
    table: Mapping[str, Any], n: int
) -> Tuple[List[mp.mpf], List[List[mp.mpf]]]:
    """The vector of per-channel ``A4`` values and its jackknife covariance.

    The combination question is a Hotelling problem in this space: the best
    linear readout has ``leverage^2 = v^T Sigma^-1 v``.  It is also the most
    overfitting-prone thing here, which is why the held-out half exists.
    """

    weights = angular_functional(n)
    width = len(BASIS_CHANNELS)

    def project(source: Mapping[str, Any]) -> List[mp.mpf]:
        return [
            mp.fsum(weights[i] * source[name][channel] for i, name in enumerate(GEOMETRY_ORDER))
            for channel in BASIS_CHANNELS
        ]

    vector = project(table["point"])
    deleted = [project(replicate) for replicate in table["replicates"]]
    covariance = jackknife_covariance(deleted)
    _ = width
    return vector, covariance


def symmetric_eigen(matrix: Sequence[Sequence[mp.mpf]]) -> Tuple[List[mp.mpf], List[List[mp.mpf]]]:
    """Jacobi eigendecomposition of a small symmetric matrix.

    Written out rather than delegated because the environment has no numpy and
    because the eigenvalue floor below has to be applied to the *actual*
    spectrum, not to a pivot heuristic.
    """

    size = len(matrix)
    a = [[mp.mpf(matrix[i][j]) for j in range(size)] for i in range(size)]
    v = [[mp.mpf(1) if i == j else mp.mpf(0) for j in range(size)] for i in range(size)]
    for _ in range(100):
        off = mp.fsum(a[i][j] ** 2 for i in range(size) for j in range(size) if i != j)
        if off <= mp.mpf("1e-60") * mp.fsum(a[i][i] ** 2 for i in range(size)):
            break
        for p in range(size - 1):
            for q in range(p + 1, size):
                if a[p][q] == 0:
                    continue
                theta = (a[q][q] - a[p][p]) / (2 * a[p][q])
                sign = mp.mpf(1) if theta >= 0 else mp.mpf(-1)
                t = sign / (abs(theta) + mp.sqrt(theta**2 + 1))
                c = 1 / mp.sqrt(t**2 + 1)
                s = t * c
                for k in range(size):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p], a[k][q] = c * akp - s * akq, s * akp + c * akq
                for k in range(size):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k], a[q][k] = c * apk - s * aqk, s * apk + c * aqk
                for k in range(size):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p], v[k][q] = c * vkp - s * vkq, s * vkp + c * vkq
    values = [a[i][i] for i in range(size)]
    vectors = [[v[i][j] for i in range(size)] for j in range(size)]
    return values, vectors


def best_combination(
    vector: Sequence[mp.mpf], covariance: Sequence[Sequence[mp.mpf]]
) -> Dict[str, Any]:
    """Leverage of the optimal linear combination, on the resolvable subspace."""

    values, vectors = symmetric_eigen(covariance)
    largest = max(values)
    if largest <= 0:
        return {"leverage": mp.mpf(0), "direction": None, "retained": 0}
    total = mp.mpf(0)
    direction = [mp.mpf(0)] * len(vector)
    retained = 0
    for value, evector in zip(values, vectors):
        if value <= COMBINATION_FLOOR * largest:
            continue
        retained += 1
        overlap = mp.fsum(evector[i] * vector[i] for i in range(len(vector)))
        total += overlap**2 / value
        for i in range(len(vector)):
            direction[i] += evector[i] * overlap / value
    norm = mp.sqrt(mp.fsum(component**2 for component in direction))
    return {
        "leverage": mp.sqrt(total),
        "direction": [component / norm for component in direction] if norm > 0 else None,
        "retained": retained,
        "condition": mp.nstr(largest / min(v for v in values if v > 0), 6)
        if any(v > 0 for v in values) else "n/a",
    }


# --------------------------------------------------------------------------
# Scan
# --------------------------------------------------------------------------


def scan_size(ca: Any, cb: Any, grid: Sequence[mp.mpf]) -> Dict[str, Any]:
    n = ca.n
    rows = {
        "C": aligned_rows(ca, "first"),
        "A": aligned_rows(ca, "second"),
        "B": aligned_rows(cb, "second"),
    }
    count = len(rows["C"])
    full = list(range(count))
    halves = {
        "half_0": [b for b in full if b % SPLIT_PARITY == 0],
        "half_1": [b for b in full if b % SPLIT_PARITY == 1],
    }
    subsets = {"full": full, **halves}

    per_subset: Dict[str, Dict[str, Any]] = {}
    for label, batches in subsets.items():
        by_point: Dict[str, Any] = {}
        for p in grid:
            table = channel_table(rows, batches, p)
            channels = {name: leverage(table, n, name) for name in BASIS_CHANNELS}
            vector, covariance = amplitude_covariance(table, n)
            by_point[mp.nstr(p, 12)] = {
                "channels": channels,
                "combination": best_combination(vector, covariance),
            }
        per_subset[label] = by_point
    return {"N": n, "batches": count, "grid": [mp.nstr(p, 12) for p in grid], "subsets": per_subset}


def baseline_key(grid: Sequence[mp.mpf]) -> str:
    """The grid point closest to the frozen p_ref, which the baseline is read at."""

    target = mp.mpf(P_REF)
    best = min(grid, key=lambda p: abs(p - target))
    if abs(best - target) > mp.mpf("1e-12"):
        raise ValueError("the predeclared grid must contain the frozen p_ref exactly")
    return mp.nstr(best, 12)


def readouts(scan: Mapping[str, Any], subset: str) -> Dict[Tuple[str, str], mp.mpf]:
    """Every (readout, p) leverage in one subset, flattened for selection."""

    flat: Dict[Tuple[str, str], mp.mpf] = {}
    for key, entry in scan["subsets"][subset].items():
        for name, value in entry["channels"].items():
            flat[(name, key)] = value["leverage"]
        flat[("combination", key)] = entry["combination"]["leverage"]
    return flat


def selection_control(scan: Mapping[str, Any], grid: Sequence[mp.mpf]) -> Dict[str, Any]:
    """Select on one half, report on the other; both directions.

    The in-sample leverage of a selected readout is meaningless on its own -- it
    is the maximum of many noisy numbers.  What is reportable is whether the
    *choice* transports: does the readout picked on half 0 still beat the
    baseline on half 1?
    """

    anchor = baseline_key(grid)
    directions = []
    for select_on, report_on in (("half_0", "half_1"), ("half_1", "half_0")):
        chosen = max(readouts(scan, select_on).items(), key=lambda item: item[1])
        (name, point), in_sample = chosen
        held_out = readouts(scan, report_on)[(name, point)]
        baseline_here = readouts(scan, report_on)[(BASELINE_CHANNEL, anchor)]
        baseline_there = readouts(scan, select_on)[(BASELINE_CHANNEL, anchor)]
        directions.append({
            "selected_on": select_on,
            "reported_on": report_on,
            "readout": name,
            "p": point,
            "in_sample_leverage": mp.nstr(in_sample, 8),
            "held_out_leverage": mp.nstr(held_out, 8),
            "in_sample_amplification": mp.nstr(in_sample / baseline_there, 6)
            if baseline_there > 0 else "inf",
            "held_out_amplification": mp.nstr(held_out / baseline_here, 6)
            if baseline_here > 0 else "inf",
            "optimism_factor": mp.nstr(in_sample / held_out, 6) if held_out > 0 else "inf",
            "transports": bool(baseline_here > 0 and held_out > baseline_here),
        })
    return {
        "baseline": f"{BASELINE_CHANNEL} @ {anchor}",
        "directions": directions,
        "both_directions_transport": all(row["transports"] for row in directions),
        "why": (
            "A readout selected on a half-sample and validated on the other half "
            "carries no selection optimism into the reported amplification.  The "
            "in-sample column is shown only so the size of the optimism is visible."
        ),
    }


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------


def cost_table(scan: Mapping[str, Any], grid: Sequence[mp.mpf]) -> Dict[str, Any]:
    """Full-sample amplification per readout, and what it does to sample cost.

    Reported for the *predeclared* readouts only.  Cost scales as the inverse
    square of the leverage, so a 2x amplification is a 4x saving -- which is the
    whole reason this scan is worth running before buying anything.
    """

    anchor = baseline_key(grid)
    flat = readouts(scan, "full")
    baseline = flat[(BASELINE_CHANNEL, anchor)]
    best_per_channel = {}
    for name in list(BASIS_CHANNELS) + ["combination"]:
        candidates = {key[1]: value for key, value in flat.items() if key[0] == name}
        point, value = max(candidates.items(), key=lambda item: item[1])
        best_per_channel[name] = {
            "best_p": point,
            "leverage": mp.nstr(value, 8),
            "leverage_at_p_ref": mp.nstr(candidates[anchor], 8),
            "amplification_vs_baseline": mp.nstr(value / baseline, 6)
            if baseline > 0 else "inf",
            "sample_cost_ratio_vs_baseline": mp.nstr((baseline / value) ** 2, 6)
            if value > 0 else "inf",
        }
    return {
        "baseline": f"{BASELINE_CHANNEL} @ {anchor}",
        "baseline_leverage": mp.nstr(baseline, 8),
        "by_readout": best_per_channel,
        "caveat": (
            "These are in-sample maxima over the grid and carry selection "
            "optimism; the held-out control is the number to plan with."
        ),
    }


def decide(controls: Mapping[int, Any], predeclared: bool) -> Dict[str, Any]:
    """Read the control.  A post-hoc grid is refused a verdict, by code.

    Widening the grid after seeing that the optimum sits on its edge is exactly
    the selection this scan is disciplined against.  Such a run is still worth
    doing -- knowing *where* the optimum is changes the design -- but it may not
    carry a verdict, and the refusal is enforced here rather than left to
    whoever writes the note.
    """

    if not predeclared:
        return {
            "verdict": "POST_HOC_GRID__EXPLORATORY_ONLY_NO_VERDICT",
            "why": (
                "The probability grid was widened after the predeclared grid was "
                "scored.  Leverage numbers from this run describe where the "
                "optimum lies; they may not be used to justify a design choice "
                "on their own (GOVERNANCE 2D)."
            ),
        }
    transports = all(controls[n]["both_directions_transport"] for n in SIZES)
    held_out = [
        mp.mpf(row["held_out_amplification"])
        for n in SIZES for row in controls[n]["directions"]
        if row["held_out_amplification"] != "inf"
    ]
    worst = min(held_out) if held_out else mp.mpf(0)
    if transports and worst > 1:
        verdict = "CHANNEL_CHOICE_TRANSPORTS__SELECT_THE_READOUT_BEFORE_BUYING_SAMPLES"
    elif worst > 1:
        verdict = "CHANNEL_CHOICE_TRANSPORTS_AT_ONE_SIZE_ONLY__NOT_YET_A_DESIGN_RULE"
    else:
        verdict = "CHANNEL_SPREAD_IS_NOT_REPRODUCIBLE__TREAT_THE_9000x_AS_NOISE"
    return {
        "both_sizes_transport": transports,
        "worst_held_out_amplification": mp.nstr(worst, 6),
        "worst_held_out_sample_cost_ratio": mp.nstr(1 / worst**2, 6) if worst > 0 else "inf",
        "verdict": verdict,
        "what_this_cannot_do": (
            "This scan reports angular resolution only.  It deliberately does "
            "not compute delta, rho, or any contamination bound for a readout it "
            "selected: those would be numbers the selection had already seen "
            "(GOVERNANCE 2E).  A selected readout must be declared "
            "prospectively and bought with new samples before it carries a verdict."
        ),
    }


def grid_boundary_binds(scan: Mapping[str, Any], grid: Sequence[mp.mpf]) -> Dict[str, bool]:
    """Which readouts put their best leverage on an endpoint of the grid.

    A binding boundary means the true optimum is outside the range that was
    looked at, so the reported amplification is a *lower* bound -- and it must
    not be quietly fixed by widening the grid inside the same artifact.
    """

    edges = {mp.nstr(grid[0], 12), mp.nstr(grid[-1], 12)}
    flat = readouts(scan, "full")
    binds = {}
    for name in list(BASIS_CHANNELS) + ["combination"]:
        candidates = {key[1]: value for key, value in flat.items() if key[0] == name}
        best = max(candidates.items(), key=lambda item: item[1])[0]
        binds[name] = best in edges
    return binds


def assemble(
    pairs: Mapping[int, Any], contract: Mapping[str, Any],
    half_width: str = P_HALF_WIDTH, points: int = P_GRID_POINTS,
) -> Dict[str, Any]:
    predeclared = (half_width == P_HALF_WIDTH and points == P_GRID_POINTS)
    grid = probability_grid(half_width, points)
    scans = {n: scan_size(pairs[n][0], pairs[n][1], grid) for n in SIZES}
    controls = {n: selection_control(scans[n], grid) for n in SIZES}
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "status": "design scan over archived histograms; angular resolution only, no verdict-bearing statistic",
        "basis_channels": list(BASIS_CHANNELS),
        "baseline_channel": BASELINE_CHANNEL,
        "p_ref": P_REF,
        "grid": {
            "half_width": half_width,
            "points": points,
            "values": [mp.nstr(p, 12) for p in grid],
            "predeclared": predeclared,
            "verdict_bearing": predeclared,
        },
        "excluded_channel": {
            "D": "exactly M/2 in obs(), so linearly dependent and singular in the combination"
        },
        "by_size": {
            str(n): {
                "batches": scans[n]["batches"],
                "cost_table": cost_table(scans[n], grid),
                "selection_control": controls[n],
                "grid_boundary_binds": grid_boundary_binds(scans[n], grid),
                "full_sample_leverage": {
                    key: {
                        name: mp.nstr(value["leverage"], 8)
                        for name, value in entry["channels"].items()
                    } | {"combination": mp.nstr(entry["combination"]["leverage"], 8),
                         "combination_retained": entry["combination"]["retained"]}
                    for key, entry in scans[n]["subsets"]["full"].items()
                },
            }
            for n in SIZES
        },
        "decision": decide(controls, predeclared),
        "provenance": {
            "inputs": [
                {
                    "N": run.n,
                    "pair": f"C-{run.partner}",
                    "histogram_sha256": sha256(run.histogram_path),
                    "git_commit": run.metadata["git_commit"],
                }
                for n in SIZES for run in pairs[n]
            ],
        },
    }


def render(payload: Mapping[str, Any]) -> str:
    lines = [f"#{ISSUE} P205 channel leverage: |A4| / se(A4) over a predeclared grid", ""]
    for n in SIZES:
        block = payload["by_size"][str(n)]
        table = block["cost_table"]
        lines.append(f"N = {n}   baseline {table['baseline']}  leverage {table['baseline_leverage'][:8]}")
        lines.append(
            f"  {'readout':<12} {'best p':>13} {'leverage':>10} {'@p_ref':>10} {'x base':>8} {'x samples':>10}"
        )
        for name, row in table["by_readout"].items():
            lines.append(
                f"  {name:<12} {row['best_p']:>13} {row['leverage'][:10]:>10}"
                f" {row['leverage_at_p_ref'][:10]:>10}"
                f" {row['amplification_vs_baseline'][:8]:>8}"
                f" {row['sample_cost_ratio_vs_baseline'][:10]:>10}"
            )
        bound = [name for name, hit in block["grid_boundary_binds"].items() if hit]
        if bound:
            lines.append(
                "  grid boundary binds for: " + ", ".join(bound)
                + "  (optimum is outside the scanned range; amplification is a lower bound)"
            )
        lines.append("  held-out selection control (select on one half, report on the other):")
        for row in block["selection_control"]["directions"]:
            lines.append(
                f"    {row['selected_on']} -> {row['reported_on']}: {row['readout']} @ {row['p']}"
                f"  in {row['in_sample_amplification'][:6]}x  held-out {row['held_out_amplification'][:6]}x"
                f"  optimism {row['optimism_factor'][:6]}x  {'TRANSPORTS' if row['transports'] else 'does not beat baseline'}"
            )
        lines.append("")
    decision = payload["decision"]
    if not payload["grid"]["predeclared"]:
        lines += ["POST-HOC GRID: exploratory only, carries no verdict.", f"({decision['verdict']})"]
        return "\n".join(lines)
    lines += [
        f"worst held-out amplification: {decision['worst_held_out_amplification']}"
        f"  (sample cost x{decision['worst_held_out_sample_cost_ratio']})",
        f"VERDICT: {decision['verdict']}",
    ]
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair", action="append", type=pair_spec, required=True)
    parser.add_argument(
        "--prediction", type=Path,
        default=ROOT / "predictions/norm5_conjugate_coalescence_20260829.yaml",
    )
    parser.add_argument(
        "--experiment", type=Path,
        default=ROOT / "experiments/p205_norm5_conjugate_coalescence_20260829.yaml",
    )
    parser.add_argument("--dps", type=int, default=40)
    parser.add_argument(
        "--half-width", default=P_HALF_WIDTH,
        help="grid half-width; anything but the predeclared value marks the run post-hoc",
    )
    parser.add_argument("--points", type=int, default=P_GRID_POINTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    mp.mp.dps = args.dps

    contract = load_contract(args.prediction, args.experiment)
    loaded = {}
    for n, partner, histogram, moments, metadata in args.pair:
        loaded[(n, partner)] = load_pair(n, partner, histogram, moments, metadata, contract)
    expected = {(n, partner) for n in SIZES for partner in ("A", "B")}
    if set(loaded) != expected:
        raise SystemExit(f"pairs must be exactly {sorted(expected)}")
    pairs = {n: (loaded[(n, "A")], loaded[(n, "B")]) for n in SIZES}

    payload = assemble(pairs, contract, args.half_width, args.points)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(render(payload))
    print()
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
