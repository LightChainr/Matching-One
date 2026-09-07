#!/usr/bin/env python3
"""P7 — declared-chart freeze, then the 4%.

Freeze ONE chart (the P7 default, not killed by P0-P6):

    declared chart  = span{1, Q_base, g_frozen}
    g_frozen        = the #582/#584 consensus (PR #614 commit f0981a98;
                      published in results/p582-amplitude-law/latest.json)
    transport       = the #612 identity, applied before any residual
    weightings      = spin0 AND equal, always together
    reconstruction  = production inverse-CDF (nine deciles, binomial-
                      profile CDF) -- P5 named this readout as load-bearing
    forbidden       = second exponent; chart change after seeing a residual

Then, and only then, recompute the leftover residual in that chart:

  - the scalar 4.25% analogue on the p50 triangle (145-290-725), both
    weightings;
  - the full residual vector (9 minus affine minus g = 6 free components)
    with the jackknife covariance;
  - a sign-flip null on batch labels (10^4 label Monte Carlo): the
    committed jackknife pseudo-values are sign-flipped batch-wise; the
    residual statistic's null distribution is what a "no shape" world
    looks like;
  - the projection of the residual onto span{1, Q, g} vs onto the
    orthogonal complement (the "perpendicular shape" question P3 raised).

If the 4% dies, say so and stop hunting.  If it survives as a
weighting-stable direction orthogonal to {1, Q, g}, hand the vector to the
owner UNNAMED: no fitting, no naming, no second exponent.
"""
from __future__ import annotations

import json
import math
import random
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import affine_gauge_common as ag  # noqa: E402
import threshold_quantile_lineage as lineage  # noqa: E402
import score_wasserstein_shape_flow as flow  # noqa: E402
import p582_amplitude_law as law  # noqa: E402
from p612_chart_identity import gls_fit, LEVELS  # noqa: E402

SCHEMA = "matching-one.probe-affine-gauge.p7.v1"
ISSUE = 617

P50_NODES = (145, 290, 725)
FROZEN_OMEGA = 0.9701786160832164
FROZEN_SCALE = 0.17834748516661017
N_PERMUTATIONS = 10_000
RANDOM_SEED = 6170425          # frozen before any run; recorded here


# ------------------------------------------------------------- machinery


def _raw_blocks(sizes: tuple[int, ...],
                weighting: str) -> dict[int, dict[str, Any]]:
    patched = dict(flow.SOURCES)
    patched[725] = ag.N725_SOURCE
    return {size: lineage.load_batch_histograms(ag.ROOT / patched[size])
            for size in sizes}


def _spin0_weights(block: dict[str, Any]) -> dict[str, float]:
    return lineage.spin_zero_weights(block["orientation_cos4theta"])


def _pseudo_values(block: dict[str, Any], weights: dict[str, float],
                   grid: tuple[float, ...]) -> dict[str, Any]:
    jack = lineage.jackknife_quantiles(block, weights,
                                       levels=list(grid))
    covariance = lineage.jackknife_covariance(jack["full"], jack["deleted"])
    return {"full": jack["full"], "deleted": jack["deleted"],
            "covariance": covariance, "batches": jack["batches"]}


def _second_difference_fit(readouts: list[dict[str, Any]],
                           nodes: tuple[float, ...],
                           attach: list[float],
                           direction: list[float]) -> dict[str, Any]:
    logs = [math.log(n) for n in nodes]
    weights = law.second_difference_weights(logs)
    rank = len(readouts[0]["full"])
    observation = [math.fsum(w * ro["full"][level]
                             for w, ro in zip(weights, readouts))
                   for level in range(rank)]
    covariance = [[math.fsum(w * w * ro["covariance"][i][j]
                             for w, ro in zip(weights, readouts))
                   for j in range(rank)] for i in range(rank)]
    fit = gls_fit(observation, covariance,
                  [[1.0] * rank, list(attach), list(direction)])
    return {
        "amplitude": fit["amplitudes"][2],
        "standard_error": fit["standard_errors"][2],
        "residual": fit["residual"],
        "statistic": fit["statistic"],
        "degrees_of_freedom": fit["covariance_rank"] - 3,
        "covectors": fit["covectors"],
        "second_difference_weights": list(weights),
    }


def _projection_split(residual: list[float],
                      attach: list[float],
                      direction: list[float],
                      covariance: list[list[float]]) -> dict[str, Any]:
    """Split the residual into the declared chart and its complement.

    The projection is in the COVARIANCE metric (the GLS metric): a vector
    is decomposed into its span{1, Q, g} part and the remainder, and each
    part's norm is reported in units of the jackknife noise of that part.
    """
    rank = len(residual)
    basis = [[1.0] * rank, list(attach), list(direction)]
    # Gram matrix in the covariance metric
    gram = [[math.fsum(basis[a][i] * covariance[i][j] * basis[b][j]
                       for i in range(rank) for j in range(rank))
             for b in range(3)] for a in range(3)]
    dual = []
    det = (gram[0][0] * (gram[1][1] * gram[2][2] - gram[1][2] * gram[2][1])
           - gram[0][1] * (gram[1][0] * gram[2][2] - gram[1][2] * gram[2][0])
           + gram[0][2] * (gram[1][0] * gram[2][1] - gram[1][1] * gram[2][0]))
    if abs(det) < 1e-300:
        return {"parallel_norm": None, "perpendicular_norm": None,
                "note": "chart Gram matrix singular"}
    adj = [[(gram[(a + 1) % 3][(b + 1) % 3]
             * gram[(a + 2) % 3][(b + 2) % 3]
             - gram[(a + 1) % 3][(b + 2) % 3]
             * gram[(a + 2) % 3][(b + 1) % 3]) / det
            for b in range(3)] for a in range(3)]
    for a in range(3):
        dual.append([math.fsum(adj[a][b] * basis[b][i]
                               for b in range(3))
                     for i in range(rank)])
    # coefficients c_a = <dual_a, r>
    coeffs = [math.fsum(dual[a][i] * covariance[i][j] * residual[j]
                        for i in range(rank) for j in range(rank))
              for a in range(3)]
    parallel = [math.fsum(coeffs[a] * basis[a][j] for a in range(3))
                for j in range(rank)]
    perpendicular = [residual[j] - parallel[j] for j in range(rank)]
    norm = lambda vec: math.sqrt(math.fsum(
        vec[i] * covariance[i][j] * vec[j]
        for i in range(rank) for j in range(rank)))
    parallel_norm = norm(parallel)
    perpendicular_norm = norm(perpendicular)
    total_norm = norm(residual)
    # noise of the perpendicular part: the residual covariance of a fit on
    # the complement -- approximate by the jackknife covariance's own
    # quadratic form with the projection's pseudoinverse; the honest scale
    # is the chi2 of the perpendicular part against the residual
    # covariance, df = rank - rank(chart).
    return {
        "parallel_norm": parallel_norm,
        "perpendicular_norm": perpendicular_norm,
        "total_norm": total_norm,
        "perpendicular_share": (perpendicular_norm / total_norm
                                if total_norm else None),
        "coefficients_on_chart": coeffs,
    }


# ------------------------------------------------------------------ main


def main() -> dict[str, Any]:
    raw_by_weighting = {}
    for weighting in ("spin0", "equal"):
        raw_by_weighting[weighting] = _raw_blocks(P50_NODES, weighting)

    committed = json.loads(
        (ag.OUT_DIR / "_reference" / "p612-n725-score-latest.json")
        .read_text())
    frozen_g = list(
        committed["weightings"]["spin0"]
        ["consensus_direction_frozen_from_five_committed_transitions"])

    out = {}
    for weighting in ("spin0", "equal"):
        raw = raw_by_weighting[weighting]
        weights = {size: _spin0_weights(block) if weighting == "spin0"
                   else {"first": 0.5, "second": 0.5}
                   for size, block in raw.items()}
        readouts = {size: _pseudo_values(raw[size], weights[size],
                                         (0.1 * i for i in range(1, 10)))
                    for size in P50_NODES}
        logs = [math.log(n) for n in P50_NODES]
        h0, h1 = logs[1] - logs[0], logs[2] - logs[1]
        factor = 2.0 / (h0 + h1)

        # ---- the declared-chart curvature (attachment Q_290) ----
        attach = readouts[290]["full"]
        curvature = _second_difference_fit(
            [readouts[145], readouts[290], readouts[725]], P50_NODES,
            attach, frozen_g)
        measured = curvature["amplitude"]

        # ---- the transported prediction (#612 identity, p50 version) ----
        def _transition(a: int, b: int, multiplier: float):
            step = math.log(multiplier)
            displacement = [(t - x) / step
                            for x, t in zip(readouts[a]["full"],
                                            readouts[b]["full"])]
            covariance = [[(p + q) / (step * step)
                           for p, q in zip(ra, rb)]
                          for ra, rb in zip(readouts[a]["covariance"],
                                            readouts[b]["covariance"])]
            fit = gls_fit(displacement, covariance,
                          [[1.0] * LEVELS, list(readouts[a]["full"]),
                           list(frozen_g)])
            return fit
        lower = _transition(145, 290, 2.0)
        upper = _transition(290, 725, 2.5)
        a0, a1 = lower["amplitudes"][2], upper["amplitudes"][2]
        beta0 = lower["amplitudes"][1]
        k = 1.0 + h0 * beta0
        model_lower = law.first_difference_image(
            FROZEN_SCALE, FROZEN_OMEGA, (logs[0] + logs[1]) / 2.0, h0)
        model_upper = law.first_difference_image(
            FROZEN_SCALE, FROZEN_OMEGA, (logs[1] + logs[2]) / 2.0, h1)
        chart_corrected = factor * (model_upper - model_lower / k)
        scalar_ratio = measured / chart_corrected
        scalar_excess_percent = 100.0 * (scalar_ratio - 1.0)

        # ---- the full residual vector, chart vs complement ----
        weights_sd = curvature["second_difference_weights"]
        cov_obs = [[math.fsum(w * w * readouts[size]["covariance"][i][j]
                              for w, size in zip(weights_sd, P50_NODES))
                    for j in range(LEVELS)] for i in range(LEVELS)]
        projection = _projection_split(curvature["residual"], attach,
                                       frozen_g, cov_obs)

        # ---- sign-flip null on batch labels (10^4 draws) ----
        # The null world: each batch's pseudo-value contribution flips
        # sign.  Operationally: rebuild the delete-one pseudo-values with a
        # random +/- on each batch, resample the curvature amplitude, and
        # collect the residual statistic's null distribution.  The
        # pseudo-value identity  pv_b = B*full - (B-1)*deleted_b  gives the
        # per-batch contribution; the null rebuilds full with random signs
        # on (pv_b - full) deltas.
        rng = random.Random(RANDOM_SEED)
        batch_count = readouts[145]["batches"]
        deltas = {
            size: [[(batch_count * readouts[size]["full"][level]
                     - (batch_count - 1) * deleted[level]
                     - readouts[size]["full"][level])
                    for level in range(LEVELS)]
                   for deleted in readouts[size]["deleted"].values()]
            for size in P50_NODES
        }
        null_statistics = []
        for _ in range(N_PERMUTATIONS):
            signs = [rng.choice((-1.0, 1.0)) for _ in range(batch_count)]
            null_readouts = []
            for size in P50_NODES:
                full = readouts[size]["full"]
                perturbed = []
                for level in range(LEVELS):
                    total = full[level] + math.fsum(
                        s * d[level] / batch_count
                        for s, d in zip(signs, deltas[size]))
                    perturbed.append(total)
                null_readouts.append({
                    "full": perturbed,
                    "covariance": readouts[size]["covariance"],
                })
            null_fit = _second_difference_fit(
                null_readouts, P50_NODES, attach, frozen_g)
            null_statistics.append(null_fit["statistic"])
        null_statistics.sort()
        observed_stat = curvature["statistic"]
        null_p = sum(1 for s in null_statistics if s >= observed_stat) \
            / N_PERMUTATIONS

        out[weighting] = {
            "declared_chart": {
                "attachment": "Q_290",
                "g_frozen": "consensus from five committed transitions "
                            "(PR #614, commit f0981a98)",
                "transport": "#612 identity",
                "reconstruction": "nine deciles, binomial-profile CDF "
                                  "(P5 named this readout)",
                "weightings": "spin0 AND equal (paired)",
            },
            "scalar": {
                "measured": measured,
                "measured_se": curvature["standard_error"],
                "chart_corrected": chart_corrected,
                "ratio_measured_over_chart_corrected": scalar_ratio,
                "excess_percent": scalar_excess_percent,
                "k_contraction": k,
            },
            "residual_vector": {
                "statistic": observed_stat,
                "degrees_of_freedom": curvature["degrees_of_freedom"],
                "projection": projection,
            },
            "sign_flip_null": {
                "permutations": N_PERMUTATIONS,
                "seed": RANDOM_SEED,
                "observed_statistic": observed_stat,
                "null_median": null_statistics[
                    N_PERMUTATIONS // 2] if null_statistics else None,
                "null_p_value": null_p,
                "null_max": null_statistics[-1] if null_statistics else None,
            },
        }

    # ---- verdicts ----
    spin0 = out["spin0"]
    equal = out["equal"]
    excess = [out[w]["scalar"]["excess_percent"] for w in ("spin0", "equal")]
    p_value = spin0["sign_flip_null"]["null_p_value"]
    perp_share = spin0["residual_vector"]["projection"][
        "perpendicular_share"]
    survives = p_value < 0.01
    weighting_stable = (abs(excess[0] - excess[1])
                        < 2.0 * spin0["scalar"]["measured_se"]
                        / spin0["scalar"]["chart_corrected"] * 100.0)
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "declared_chart": spin0["declared_chart"],
        "weightings": out,
        "excess_percent_both_weightings": excess,
        "residual_survives_null": survives,
        "residual_null_p_value": p_value,
        "perpendicular_share_spin0": perp_share,
        "weighting_stable": weighting_stable,
        "verdict": (
            "the 4%% SURVIVES the sign-flip null (p=%.4f) as a "
            "perpendicular shape component in the declared chart; deliver "
            "the vector UNNAMED -- no fit, no name, no second exponent."
            % p_value if survives else
            "the 4%% DIES under the sign-flip null (p=%.4f): the leftover "
            "is consistent with the null world of label noise.  Stop "
            "hunting; write the death and do not replace it." % p_value),
    }


if __name__ == "__main__":
    result = main()
    path = ag.dump("p7-residual", result)
    print(json.dumps({
        "excess_percent": result["excess_percent_both_weightings"],
        "null_p_value": result["residual_null_p_value"],
        "perpendicular_share": result["perpendicular_share_spin0"],
        "verdict": result["verdict"],
        "output": str(path)}, indent=2))
