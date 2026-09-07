#!/usr/bin/env python3
"""Type #582's structured remainder with pre-existing discrete labels -- issue #584 (Gate 3).

#582 established that the threshold law's finite-size shape flow is one dominant
transferable direction ``g_N(u)`` carrying ~99% of the covariance-weighted
residual, plus a small, strongly resolved, structured remainder ``r_N(u)``.  This
module asks #584's question with #584's own discipline: index that remainder by a
**pre-existing** discrete label before adding any rank, and score every candidate
against permutation / random-label nulls rather than against optimism.

Three things happen here, in order.

1. **Re-freeze.**  ``score_wasserstein_shape_flow`` does not persist the per-size
   9x9 delete-one covariance ``S_N`` (only its rank and condition), so this module
   re-reads the committed histograms and re-derives ``S_N``, the affine-orthogonal
   residuals, the consensus ``g_N``, and each transition's ``r_N``.  It reuses the
   committed pipeline, not a reimplementation of it.

2. **Taylor-curvature null** (the tightening added on #584).  A first finite-size
   displacement is trivially near rank one for a smooth one-parameter family; the
   second-difference curvature is the part that survives.  Where a lineage has
   three sizes this computes the curvature and asks how much of the remainder it
   predicts.  ``p50`` has two sizes and cannot supply a curvature.

3. **The label screen.**  Each candidate label is a partition of the five
   transitions.  For each partition the module reports within- vs between-class
   principal angles, an exact permutation null (every relabelling consistent with
   the class sizes), leave-one-transition prediction, and the residual norm left
   after a class-conditional one-direction correction.  The verdict is #584's
   decision table, not this module's opinion.

Nothing is sampled; every number is a reanalysis of committed productions.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

try:  # pragma: no cover
    from scripts.projective_inference import subspace_residual, spectral_pseudo_inverse
    from scripts import threshold_quantile_lineage as lineage
    import scripts.score_wasserstein_shape_flow as sf
except ModuleNotFoundError:  # pragma: no cover
    from projective_inference import subspace_residual, spectral_pseudo_inverse
    import threshold_quantile_lineage as lineage
    import score_wasserstein_shape_flow as sf

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "type582-residual" / "latest.json"
SCHEMA = "matching-one.type582-residual.v1"
ISSUE = 584

#: The parent Gaussian prime of each lineage, as a canonical ``(a, b)`` with
#: ``a >= b >= 0``.  gaussian_13 -> 3 + 2i (13 = 3^2 + 2^2), gaussian_17 ->
#: 4 + i (17 = 4^2 + 1^2), p50 -> 5 + 2i (29 = 5^2 + 2^2; the lineage's sizes
#: 145 = 5 * 29 and 290 = 2 * 5 * 29 both carry the norm-29 prime).
GAUSSIAN_PRIMES: dict[str, tuple[int, int]] = {
    "gaussian_13": (3, 2),
    "gaussian_17": (4, 1),
    "p50": (5, 2),
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _normalise(vector: Sequence[float]) -> list[float]:
    norm = math.sqrt(sum(value * value for value in vector))
    _require(norm > 0.0, "cannot normalise a zero direction")
    return [value / norm for value in vector]


def _dot(left: Sequence[float], right: Sequence[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def _acute_angle_degrees(left: Sequence[float], right: Sequence[float]) -> float:
    """Principal angle between two directions, folded into [0, 90]."""
    cosine = max(-1.0, min(1.0, _dot(_normalise(left), _normalise(right))))
    return math.degrees(math.acos(abs(cosine)))


def _average_transition_covariance(loaded: Mapping[int, Mapping[str, Any]],
                                   rows: Sequence[Mapping[str, Any]]) -> list[list[float]]:
    """Mean of the five transitions' ``(S_base + S_target)/log(m)^2``.

    Each transition has its own covariance; the screen needs one common
    weighting to define a single covariance-identifiable subspace, so the mean is
    used.  This is the same per-transition covariance ``shape_flow`` uses.
    """
    matrices = []
    for row in rows:
        base = loaded[row["base"]]["covariance"]
        target = loaded[row["target"]]["covariance"]
        scale = math.log(row["multiplier"]) ** 2
        matrices.append([[(b + t) / scale for b, t in zip(br, tr)]
                         for br, tr in zip(base, target)])
    size = len(matrices[0])
    return [[sum(m[i][j] for m in matrices) / len(matrices)
             for j in range(size)] for i in range(size)]


def _weighted_acute_angle(left: Sequence[float], right: Sequence[float], pinv) -> float:
    """Covariance-weighted principal angle, folded into [0, 90].

    ``pinv`` is ``spectral_pseudo_inverse`` of the common covariance, so the
    angle lives in the covariance-identifiable subspace (null directions of the
    covariance carry no weight) rather than in raw L2.
    """
    def quadratic(u, v):
        return float(sum(u[i] * sum(pinv[i, j] * v[j] for j in range(len(v)))
                         for i in range(len(u))))
    num = quadratic(left, right)
    den = math.sqrt(quadratic(left, left) * quadratic(right, right))
    if den <= 0.0:
        return 90.0
    cosine = max(-1.0, min(1.0, num / den))
    return math.degrees(math.acos(abs(cosine)))


# ---------------------------------------------------------------------------
# Part A -- re-freeze
# ---------------------------------------------------------------------------

def all_sizes() -> list[int]:
    return sorted({size for pair in sf.LINEAGES.values() for size in pair})


def refreeze() -> tuple[dict[int, dict[str, Any]], list[dict[str, Any]], list[float]]:
    """Load spin0 quantiles + covariance, score transitions, extract consensus g_N.

    Returns ``(loaded, rows, g_N)`` where ``rows`` carries each transition's
    affine-orthogonal ``shape_flow`` result under the key ``affine_only``.
    """
    loaded = sf.load_sizes(all_sizes())["spin0"]
    rows: list[dict[str, Any]] = []
    for meta in sf.transitions():
        base, target = loaded[meta["base"]], loaded[meta["target"]]
        rows.append({**meta, "affine_only": sf.shape_flow(base, target, meta["multiplier"])})
    residuals = [_normalise(row["affine_only"]["residual"]) for row in rows]
    consensus = sf._leading_direction(residuals)
    return loaded, rows, consensus


def residual_after(loaded: Mapping[int, Mapping[str, Any]],
                   row: Mapping[str, Any],
                   direction: Sequence[float]) -> dict[str, Any]:
    """Score a transition against ``span{1, Q_base, direction}``.

    The returned ``residual`` is the remainder after the affine tangent *and* the
    frozen direction are projected out; ``amplitudes[2]`` is that direction's
    coefficient.  This is #582's held-out prediction with a caller-chosen
    direction, so the consensus ``g_N`` gives ``r_N`` directly.
    """
    base = loaded[row["base"]]
    target = loaded[row["target"]]
    fit = sf.shape_flow(base, target, row["multiplier"], [list(direction)])
    fit["amplitude_of_direction"] = fit["amplitudes"][2]
    return fit


# ---------------------------------------------------------------------------
# Part B -- Taylor-curvature null
# ---------------------------------------------------------------------------

def _raw_velocity(base: Sequence[float], target: Sequence[float],
                  multiplier: float) -> list[float]:
    """The unprojected displacement ``(Q_M - Q_N) / log(m)``, per frozen level."""
    scale = math.log(multiplier)
    return [(t - b) / scale for b, t in zip(base, target)]


def taylor_curvature(loaded: Mapping[int, Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Second-difference curvature per three-size lineage, after affine projection.

    For sizes ``N0 < N1 < N2`` the two per-log-step velocities are
    ``v01 = (Q(N1)-Q(N0))/log(N1/N0)`` and ``v12 = (Q(N2)-Q(N1))/log(N2/N1)``.
    The curvature is ``(v12 - v01) / mean_step`` with
    ``mean_step = (log(N1/N0) + log(N2/N1)) / 2``, an unequal-step divided
    difference for ``d^2 Q / ds^2``.  It is projected against ``span{1, Q(N0)}``
    so only the affine-orthogonal curvature is reported.  Lineages with two sizes
    cannot supply a curvature and are skipped with a note.
    """
    out: list[dict[str, Any]] = []
    for name, sizes in sf.LINEAGES.items():
        if len(sizes) < 3:
            out.append({"lineage": name, "sizes": list(sizes),
                        "computable": False,
                        "reason": "two sizes give one step, not a second difference"})
            continue
        n0, n1, n2 = sizes
        q0 = loaded[n0]["quantiles"]
        q1 = loaded[n1]["quantiles"]
        q2 = loaded[n2]["quantiles"]
        h01 = math.log(n1 / n0)
        h12 = math.log(n2 / n1)
        mean_step = 0.5 * (h01 + h12)
        v01 = _raw_velocity(q0, q1, n1 / n0)
        v12 = _raw_velocity(q1, q2, n2 / n1)
        curvature = [(b - a) / mean_step for a, b in zip(v01, v12)]
        # Project out span{1, Q(N0)} so the reported curvature is shape-only.  The
        # projection uses the base size's covariance for the weighting; the full
        # three-size curvature covariance is not propagated here because this is a
        # direction estimate for a null comparison, not a hypothesis test.
        covariance = loaded[n0]["covariance"]
        basis = [[1.0] * len(curvature), list(q0)]
        fit = subspace_residual(curvature, covariance, basis)
        fitted = [sum(fit["amplitudes"][k] * basis[k][j] for k in range(len(basis)))
                  for j in range(len(curvature))]
        shape = [curvature[j] - fitted[j] for j in range(len(curvature))]
        out.append({
            "lineage": name,
            "sizes": [n0, n1, n2],
            "computable": True,
            "step_logs": [h01, h12],
            "curvature_shape": shape,
            "curvature_norm": math.sqrt(_dot(shape, shape)),
            "curvature_statistic": fit["statistic"],
            "curvature_degrees_of_freedom": fit["degrees_of_freedom"],
        })
    return out


# ---------------------------------------------------------------------------
# Part C -- labels
# ---------------------------------------------------------------------------

def _gaussian_norm(point: tuple[int, int]) -> int:
    return point[0] ** 2 + point[1] ** 2


def _label_table(loaded: Mapping[int, Mapping[str, Any]],
                 rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """One entry per transition, with every pre-existing discrete label.

    Labels are computed from exact arithmetic (the Gaussian primes and the
    committed per-size interpolation flags), never from the residuals, so the
    partition is fixed before any screen is run.
    """
    out = []
    for row in rows:
        lineage_name = row["lineage"]
        a, b = GAUSSIAN_PRIMES[lineage_name]
        cos4 = lineage.cos_four_theta(a, b)
        dominant = max(a, b)
        out.append({
            "transition": row["label"],
            "lineage": lineage_name,
            "multiplier": row["multiplier"],
            "gaussian_prime": [a, b],
            "prime_norm": _gaussian_norm((a, b)),
            "prime_dominant_parity": "odd" if dominant % 2 == 1 else "even",
            "prime_min_component": min(a, b),
            "prime_cos4": cos4,
            "prime_cos4_sign": "positive" if cos4 >= 0.0 else "negative",
            "target_interpolation": loaded[row["target"]]["spin4_correction_is_an_interpolation"],
        })
    return out


#: Candidate partitions.  Each is a callable mapping a label-table entry to a
#: class string.  The reference labels (lineage, multiplier) are kept so the
#: screen can confirm the #582 observation that they do not index the remainder;
#: the arithmetic labels are the genuinely new ones.
def curvature_alignment(curvatures: Sequence[Mapping[str, Any]],
                        rows: Sequence[Mapping[str, Any]],
                        r_vectors: Sequence[Sequence[float]],
                        g_N: Sequence[float]) -> list[dict[str, Any]]:
    """How the Taylor curvature aligns with ``g_N`` and with each remainder.

    The Taylor null asks whether the remainder is ordinary second-order curvature.
    If the remainder direction is close to its own lineage's curvature, the
    remainder is smooth curvature, not a discrete fiber; if the two lineages'
    curvatures point in unrelated directions while each remainder points at its
    own curvature, that is the smooth-curvature explanation.
    """
    out = []
    for curv in curvatures:
        if not curv["computable"]:
            out.append({"lineage": curv["lineage"], "computable": False})
            continue
        name = curv["lineage"]
        shape = curv["curvature_shape"]
        members = [(i, row) for i, row in enumerate(rows) if row["lineage"] == name]
        out.append({
            "lineage": name,
            "computable": True,
            "angle_to_g_deg": _acute_angle_degrees(shape, g_N),
            "remainder_angles_deg": [
                {"transition": row["label"],
                 "angle_to_curvature_deg": _acute_angle_degrees(shape, r_vectors[i])}
                for i, row in members],
        })
    return out


def _partition_functions() -> dict[str, Any]:
    return {
        "lineage": lambda e: e["lineage"],
        "multiplier": lambda e: f"m={e['multiplier']:g}",
        "prime_dominant_parity": lambda e: e["prime_dominant_parity"],
        "prime_min_component": lambda e: f"min={e['prime_min_component']}",
        "prime_cos4_sign": lambda e: e["prime_cos4_sign"],
        "target_interpolation": lambda e: ("interpolation" if e["target_interpolation"]
                                           else "extrapolation"),
    }


# ---------------------------------------------------------------------------
# Part D -- the screen
# ---------------------------------------------------------------------------

def _assign(entries: Sequence[Mapping[str, Any]], key) -> list[str]:
    return [key(entry) for entry in entries]


def _partition_summary(assignments: Sequence[str]) -> dict[str, Any]:
    classes = sorted(set(assignments))
    return {
        "classes": classes,
        "class_sizes": {c: assignments.count(c) for c in classes},
        "assignment": list(assignments),
    }


def _separation(angles: Sequence[Sequence[float]], assignments: Sequence[str]) -> dict[str, Any]:
    """Within- vs between-class principal angles for a fixed assignment."""
    n = len(assignments)
    within = [angles[i][j] for i in range(n) for j in range(i + 1, n)
              if assignments[i] == assignments[j]]
    between = [angles[i][j] for i in range(n) for j in range(i + 1, n)
               if assignments[i] != assignments[j]]
    mean_within = sum(within) / len(within) if within else float("nan")
    mean_between = sum(between) / len(between) if between else float("nan")
    return {
        "mean_within_deg": mean_within,
        "mean_between_deg": mean_between,
        "separation_deg": (mean_between - mean_within)
        if (within and between) else float("nan"),
        "within_pairs": len(within),
        "between_pairs": len(between),
    }


def permutation_null(angles: Sequence[Sequence[float]],
                     assignments: Sequence[str]) -> dict[str, Any]:
    """Exact permutation null over every relabelling with the same class sizes.

    With five transitions the number of distinct relabellings is at most 5! and
    is enumerated exactly, so the p-value needs no asymptotic shortcut.
    """
    n = len(assignments)
    labels = list(assignments)
    observed = _separation(angles, labels)["separation_deg"]
    counts: dict[str, int] = {}
    seen: set[tuple[str, ...]] = set()
    for perm in itertools.permutations(labels):
        seen.add(perm)
    better = 0
    total = len(seen)
    for perm in seen:
        value = _separation(angles, list(perm))["separation_deg"]
        if value >= observed:
            better += 1
    return {
        "observed_separation_deg": observed,
        "relabellings": total,
        "relabellings_at_least_as_separated": better,
        "p_value": better / total,
    }


def leave_one_out(angles: Sequence[Sequence[float]],
                  assignments: Sequence[str]) -> dict[str, Any]:
    """Leave-one-transition classification: nearest-class centroid by angle.

    The centroid of a class is its first principal direction (power iteration on
    its members' unit residuals is not needed here -- the class centroid is the
    direction minimising the sum of squared principal angles to its members, which
    for a pairwise-angle screen is approximated by the member with the smallest
    mean angle to the class).  The held-out transition is assigned to the class
    whose *other* members it aligns with best.
    """
    n = len(assignments)
    correct = 0
    per_transition = []
    for held in range(n):
        best_class, best_angle = None, float("inf")
        for cls in sorted(set(assignments)):
            others = [j for j in range(n) if assignments[j] == cls and j != held]
            if not others:
                continue
            mean_angle = sum(angles[held][j] for j in others) / len(others)
            if mean_angle < best_angle:
                best_angle, best_class = mean_angle, cls
        is_correct = best_class == assignments[held]
        correct += int(is_correct)
        per_transition.append({"held_out": held, "assigned": best_class,
                               "true": assignments[held],
                               "correct": is_correct, "mean_angle_deg": best_angle})
    return {"correct": correct, "total": n, "per_transition": per_transition}


def class_conditional_correction(loaded: Mapping[int, Mapping[str, Any]],
                                 rows: Sequence[Mapping[str, Any]],
                                 g_N: Sequence[float],
                                 assignments: Sequence[str]) -> list[dict[str, Any]]:
    """Residual norm left after a class-conditional one-direction correction.

    For each transition the correction direction is the leading direction of the
    *other* members of its class (leave-one-out), fitted on the affine-orthogonal
    residual; the report is how much of that transition's ``chi^2`` the
    class-conditional direction removes, compared to the frozen ``g_N``.
    """
    out = []
    n = len(rows)
    residuals = [row["affine_only"]["residual"] for row in rows]
    for i in range(n):
        cls = assignments[i]
        others = [j for j in range(n) if assignments[j] == cls and j != i]
        if not others:
            # A singleton class has no co-class direction to fit; there is no
            # class-conditional correction, which is itself the finding.
            affine = rows[i]["affine_only"]["statistic"]
            g_fit = sf.shape_flow(loaded[rows[i]["base"]], loaded[rows[i]["target"]],
                                  rows[i]["multiplier"], [list(g_N)])
            out.append({
                "transition": rows[i]["label"],
                "class": cls,
                "class_size": assignments.count(cls),
                "affine_statistic": affine,
                "class_conditional_statistic": None,
                "class_conditional_fraction_removed": None,
                "note": "singleton class, no co-class direction",
                "consensus_g_fraction_removed": 1.0 - g_fit["statistic"] / affine,
            })
            continue
        correction = sf._leading_direction([_normalise(residuals[j]) for j in others])
        fit = sf.shape_flow(loaded[rows[i]["base"]], loaded[rows[i]["target"]],
                            rows[i]["multiplier"], [list(correction)])
        g_fit = sf.shape_flow(loaded[rows[i]["base"]], loaded[rows[i]["target"]],
                              rows[i]["multiplier"], [list(g_N)])
        affine = rows[i]["affine_only"]["statistic"]
        out.append({
            "transition": rows[i]["label"],
            "class": cls,
            "class_size": assignments.count(cls),
            "affine_statistic": affine,
            "class_conditional_statistic": fit["statistic"],
            "class_conditional_fraction_removed": 1.0 - fit["statistic"] / affine,
            "consensus_g_fraction_removed": 1.0 - g_fit["statistic"] / affine,
        })
    return out


def screen(loaded: Mapping[int, Mapping[str, Any]],
           rows: Sequence[Mapping[str, Any]],
           g_N: Sequence[float],
           label_entries: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Run the full screen for every candidate partition."""
    n = len(rows)

    def angles_from(vectors: Sequence[Sequence[float]]) -> list[list[float]]:
        return [[_acute_angle_degrees(vectors[i], vectors[j]) for j in range(n)]
                for i in range(n)]

    # Angles of the *remainder* r_N (after the consensus g_N is removed) -- this
    # is the object #584 asks to type.  The primary screen uses the
    # covariance-identifiable (weighted) angle; the unweighted L2 angle is kept
    # as a transparent reference comparable to #582's note.
    r_N = [residual_after(loaded, row, g_N)["residual"] for row in rows]
    common_covariance = _average_transition_covariance(loaded, rows)
    pinv, _, _, _ = spectral_pseudo_inverse(common_covariance)
    unweighted_angles = angles_from(r_N)
    weighted_angles = [[_weighted_acute_angle(r_N[i], r_N[j], pinv) for j in range(n)]
                       for i in range(n)]

    results = []
    for name, key in _partition_functions().items():
        assignments = _assign(label_entries, key)
        summary = _partition_summary(assignments)
        sep = _separation(weighted_angles, assignments)
        sep_unweighted = _separation(unweighted_angles, assignments)
        null = permutation_null(weighted_angles, assignments)
        loo = leave_one_out(weighted_angles, assignments)
        correction = class_conditional_correction(loaded, rows, g_N, assignments)
        results.append({
            "label": name,
            "partition": summary,
            "within_between": sep,
            "within_between_unweighted": sep_unweighted,
            "permutation_null": null,
            "leave_one_out": loo,
            "class_conditional_correction": correction,
        })
    return results


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def assemble() -> dict[str, Any]:
    loaded, rows, g_N = refreeze()
    label_entries = _label_table(loaded, rows)
    r_fits = [residual_after(loaded, row, g_N) for row in rows]
    r_vectors = [fit["residual"] for fit in r_fits]

    transitions_out = []
    for row, entry, r_fit in zip(rows, label_entries, r_fits):
        transitions_out.append({
            "transition": row["label"],
            "lineage": row["lineage"],
            "base": row["base"],
            "target": row["target"],
            "multiplier": row["multiplier"],
            "affine_statistic": row["affine_only"]["statistic"],
            "affine_degrees_of_freedom": row["affine_only"]["degrees_of_freedom"],
            "residual_norm": math.sqrt(_dot(row["affine_only"]["residual"],
                                            row["affine_only"]["residual"])),
            "remainder_statistic_after_g": r_fit["statistic"],
            "remainder_degrees_of_freedom_after_g": r_fit["degrees_of_freedom"],
            "g_amplitude": r_fit["amplitude_of_direction"],
            "labels": entry,
        })

    curvature = taylor_curvature(loaded)
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "levels": list(lineage.FROZEN_LEVELS),
        "primary_weighting": "spin0",
        "consensus_g": list(g_N),
        "transitions": transitions_out,
        "pairwise_remainder_angles_deg": [
            [round(_acute_angle_degrees(r_vectors[i], r_vectors[j]), 4)
             for j in range(len(r_vectors))]
            for i in range(len(r_vectors))],
        "taylor_curvature": curvature,
        "taylor_curvature_alignment": curvature_alignment(
            curvature, rows, r_vectors, g_N),
        "screen": screen(loaded, rows, g_N, label_entries),
        "degenerate_or_missing_labels": [
            "Smith normal form / cyclic-noncyclic: degenerate for all eight sizes "
            "(every orientation is primitive with gcd=1, giving the cyclic Z x Z/N "
            "quotient; the noncyclic cases only appear at N=650/260/340)",
            "deck-group type / Gaussian cover word: no per-size label exists in the "
            "tree for these eight sizes",
            "primitive homology sector: available but one distinct value per "
            "orientation, i.e. too fine-grained to test with five transitions",
        ],
        "not_established": [
            "that a label beating the permutation null here is a state or a fiber; "
            "#588 has shown a structured remainder can be context or memory",
            "that the p50 lineage's absence of a second-difference curvature is "
            "innocuous; it limits the Taylor null to two of three lineages",
        ],
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    payload = assemble()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
