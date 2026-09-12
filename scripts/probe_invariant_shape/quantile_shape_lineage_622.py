#!/usr/bin/env python3
"""#622 (2026-09-12 re-scope): within-model quantile-shape symmetry and the
full-vector finite-size motion of the same-observable N=145,290,725 lineage.

This is a C2 reanalysis of already-committed histogram blocks.  No Monte Carlo,
no GPU, no exponent fit.  It rebuilds, per size and per orientation, the pooled
inverse CDF ``Q`` on the nine frozen deciles, forms the declared spin-0 weighted
``Q`` (weights sum to one and annihilate ``cos 4theta``), then the anchored shape

    W_N = Q_N(0.8) - Q_N(0.2) > 0
    Z_N(u) = [Q_N(u) - Q_N(0.2)] / W_N
    A_N(u) = Z_N(u) + Z_N(1-u) - 1

with ``A_N(u)`` the reflection residual of the *normalized* shape (this is the
corrected diagnostic of ``notes/exact-foundations-correction-20260912.md`` sec.3,
NOT ``M_N(p)+M_N(1-p)``).

Order of operations, which the ticket fixes and the #628 probe got wrong:
each orientation's ``Q_i`` is built first and only then combined,
``Q_spin0 = sum_i w_i Q_i``.  CDFs are never mixed before inversion.

Per aligned delete-one batch the ENTIRE Q -> weighting -> W -> Z -> A map is
repeated, so the covariance of every functional is a covariance of one random
object across the grid, not a per-level variance.  Cross-size covariance is
taken as zero on recorded random-stream provenance (distinct seeds), never on
matching batch IDs.

Deterministic redundancies are handled by using independent coordinates:
``A(0.2)=A(0.8)=0`` and ``A(1-u)=A(u)``, so the diagnostic lives on
u = 0.1, 0.3, 0.4, 0.5; ``Z(0.2)=0`` and ``Z(0.8)=1``, so the ΔZ test lives on
the seven remaining deciles.  No pseudoinverse is used to hide a residual.

Reuses ``scripts/threshold_quantile_lineage.py`` (git blob 3b328f29, the #655
lineage: corrected cos 4theta, invert-then-weight, pooled + aligned delete-one).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import threshold_quantile_lineage as L  # noqa: E402

DECILES = tuple(float(u) for u in L.FROZEN_LEVELS)  # 0.1 .. 0.9
ANCHOR_A, ANCHOR_B = 0.2, 0.8
#: Independent A coordinates: A(0.2)=A(0.8)=0 exactly and A(1-u)=A(u).
A_COORDS = (0.1, 0.3, 0.4, 0.5)
#: Z coordinates free of the two fixed anchor values Z(0.2)=0, Z(0.8)=1.
Z_COORDS = (0.1, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9)

DEFAULT_INPUTS = {
    145: {
        "hist": "results/server-20260829/P50-n145-n290-fullcurve/raw/n145_100m.hist.csv",
        "metadata": "results/server-20260829/P50-n145-n290-fullcurve/raw/n145_100m.metadata.json",
    },
    290: {
        "hist": "results/server-20260829/P50-n145-n290-fullcurve/raw/n290_100m.hist.csv",
        "metadata": "results/server-20260829/P50-n145-n290-fullcurve/raw/n290_100m.metadata.json",
    },
    725: {
        "hist": "results/server-20260907/P612-n725-fullcurve/raw/n725_100m.hist.csv",
        "metadata": "results/server-20260907/P612-n725-fullcurve/raw/n725_100m.metadata.json",
    },
}

SOURCE_REVISIONS = {
    "base_branch": "main",
    "base_commit": "6edf775e240056dc95e373344694eed527e48872",
    "threshold_quantile_lineage.py": {
        "git_blob": "3b328f29b9ccc75d0183ea72b8622b53bbe955d7",
        "git_commit": "6edf775e240056dc95e373344694eed527e48872",
        "inherits": "PR #655 (analysis/633-n725-spin0) lineage module",
        "properties": [
            "cos_four_theta = (a^4 - 6 a^2 b^2 + b^4) / N^2",
            "invert Q per orientation, then combine",
            "pooled full vector plus aligned delete-one",
        ],
    },
    "correction": "notes/exact-foundations-correction-20260912.md sec.3 (#702)",
    "prior_probe": "PR #666 rerun-20260908.json (T1-T4), superseded by this re-scope",
}

# --------------------------------------------------------------------------
# tiny pure-python linear algebra (the compute image has no numpy)
# --------------------------------------------------------------------------


def jacobi_eigenvalues(matrix: list[list[float]], tol: float = 1e-13,
                       sweeps: int = 100) -> list[float]:
    """Eigenvalues of a real symmetric matrix, ascending (cyclic Jacobi)."""
    n = len(matrix)
    a = [row[:] for row in matrix]
    for _ in range(sweeps):
        off = math.sqrt(sum(a[i][j] ** 2 for i in range(n) for j in range(n)
                           if i != j))
        if off < tol:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-300:
                    continue
                theta = 0.5 * (a[q][q] - a[p][p]) / a[p][q]
                t = (1.0 if theta >= 0 else -1.0) / (
                    abs(theta) + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
    return sorted(a[i][i] for i in range(n))


def cholesky_solve(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    """Solve ``matrix x = rhs`` for symmetric positive definite ``matrix``."""
    n = len(matrix)
    lower = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            acc = matrix[i][j] - sum(lower[i][k] * lower[j][k]
                                     for k in range(j))
            if i == j:
                if acc <= 0.0:
                    raise ValueError(
                        f"covariance is not positive definite at pivot {i}: {acc!r}")
                lower[i][j] = math.sqrt(acc)
            else:
                lower[i][j] = acc / lower[j][j]
    y = [0.0] * n
    for i in range(n):
        y[i] = (rhs[i] - sum(lower[i][k] * y[k] for k in range(i))) / lower[i][i]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - sum(lower[k][i] * x[k] for k in range(i + 1, n))
                ) / lower[i][i]
    return x


def _gser(a: float, x: float) -> float:
    gln = math.lgamma(a)
    if x <= 0.0:
        return 0.0
    ap, total, term = a, 1.0 / a, 1.0 / a
    for _ in range(2000):
        ap += 1.0
        term *= x / ap
        total += term
        if abs(term) < abs(total) * 1e-15:
            break
    return total * math.exp(-x + a * math.log(x) - gln)


def _gcf(a: float, x: float) -> float:
    gln = math.lgamma(a)
    tiny = 1e-300
    b, c, d = x + 1.0 - a, 1.0 / tiny, 1.0 / (x + 1.0 - a)
    h = d
    for i in range(1, 2000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return math.exp(-x + a * math.log(x) - gln) * h


def gammq(a: float, x: float) -> float:
    """Upper regularized incomplete gamma Q(a, x)."""
    if x < 0.0 or a <= 0.0:
        raise ValueError("gammq domain")
    if x == 0.0:
        return 1.0
    if x < a + 1.0:
        return 1.0 - _gser(a, x)
    return _gcf(a, x)


def chi2_survival(statistic: float, dof: int) -> float:
    """Nominal chi-square survival ``P(X >= statistic)`` for integer dof."""
    return gammq(0.5 * dof, 0.5 * statistic)


# --------------------------------------------------------------------------
# shape pipeline
# --------------------------------------------------------------------------


def _index(levels: tuple[float, ...]) -> dict[float, int]:
    return {round(float(u), 9): i for i, u in enumerate(levels)}


def z_chart(Q: list[float], levels: tuple[float, ...], a: float, b: float
            ) -> tuple[list[float], float]:
    idx = _index(levels)
    qa, qb = Q[idx[round(a, 9)]], Q[idx[round(b, 9)]]
    width = qb - qa
    return [(q - qa) / width for q in Q], width


def reflection_residual(Z: list[float], levels: tuple[float, ...]) -> list[float]:
    idx = _index(levels)
    out = []
    for u in levels:
        j = idx[round(1.0 - float(u), 9)]
        out.append(Z[idx[round(float(u), 9)]] + Z[j] - 1.0)
    return out


def _subvector(vector: list[float], levels: tuple[float, ...],
               coords: tuple[float, ...]) -> list[float]:
    idx = _index(levels)
    return [vector[idx[round(u, 9)]] for u in coords]


def _submatrix(matrix: list[list[float]], levels: tuple[float, ...],
               coords: tuple[float, ...]) -> list[list[float]]:
    idx = _index(levels)
    picks = [idx[round(u, 9)] for u in coords]
    return [[matrix[i][j] for j in picks] for i in picks]


def _norm_and_se(vector: list[float], covariance: list[list[float]]
                 ) -> dict[str, float]:
    """Euclidean norm and its delta-method SE for X ~ N(vector, covariance)."""
    norm = math.sqrt(sum(x * x for x in vector))
    quad = sum(vector[i] * covariance[i][j] * vector[j]
               for i in range(len(vector)) for j in range(len(vector)))
    if norm <= 0.0:
        return {"value": 0.0, "se": float("nan"), "undefined_se": True}
    return {"value": norm, "se": math.sqrt(max(0.0, quad)) / norm,
            "undefined_se": False}


def _quadratic_diagnostic(vector: list[float], covariance: list[list[float]],
                          coords: tuple[float, ...], levels: tuple[float, ...]
                          ) -> dict:
    """Covariance-aware zero test on independent coordinates.

    The deterministic redundancies (``A(0.2)=A(0.8)=0``, ``A(1-u)=A(u)``) are
    removed by the coordinate choice, so the covariance is invertible.  It can
    still be *near-collinear*, because a batch's A fluctuation is close to one
    amplitude mode; the full-inverse chi-square is therefore reported next to
    the correlation-free (diagonal) chi-square.  If the full value were an
    artefact of inverting a nearly singular matrix it would exceed the diagonal
    one; when it does not, the rejection is carried by the measured residuals.
    """
    sub = _subvector(vector, levels, coords)
    cov = _submatrix(covariance, levels, coords)
    eigs = jacobi_eigenvalues(cov)
    smallest, largest = eigs[0], eigs[-1]
    diagonal = [cov[i][i] for i in range(len(sub))]
    standardized = [sub[i] / math.sqrt(diagonal[i]) if diagonal[i] > 0 else None
                    for i in range(len(sub))]
    correlation = [[cov[i][j] / math.sqrt(diagonal[i] * diagonal[j])
                    if diagonal[i] > 0 and diagonal[j] > 0 else None
                    for j in range(len(sub))] for i in range(len(sub))]
    diagonal_chi2 = sum(sub[i] * sub[i] / diagonal[i]
                        for i in range(len(sub)) if diagonal[i] > 0)
    out = {
        "coordinates": list(coords),
        "values": sub,
        "standardized_values": standardized,
        "covariance": cov,
        "correlation_matrix": correlation,
        "standard_errors": [math.sqrt(max(0.0, d)) for d in diagonal],
        "eigenvalues_ascending": eigs,
        "smallest_eigenvalue": smallest,
        "largest_eigenvalue": largest,
        "condition_number": (largest / smallest) if smallest > 0 else None,
        "positive_definite": smallest > 0.0,
        "diagonal_chi2_ignoring_correlations": diagonal_chi2,
        "diagonal_dof": len(sub),
    }
    if smallest <= 0.0:
        # Deterministic or numerically exact redundancy: do NOT pseudoinvert.
        out.update({
            "chi2": None,
            "dof": None,
            "p_nominal_gaussian": None,
            "reason": "singular covariance; no pseudoinverse used, direction "
                      "reported instead of a rank-reduced score",
        })
        return out
    try:
        solution = cholesky_solve(cov, sub)
    except ValueError as exc:
        out.update({"chi2": None, "dof": None, "p_nominal_gaussian": None,
                    "reason": f"cholesky failed: {exc}"})
        return out
    statistic = sum(sub[i] * solution[i] for i in range(len(sub)))
    out.update({
        "chi2": statistic,
        "dof": len(coords),
        "p_nominal_gaussian": chi2_survival(statistic, len(coords)),
        "test_label": "nominal Gaussian reference with estimated delete-one covariance",
        "full_inverse_not_larger_than_diagonal": statistic <= diagonal_chi2,
        "max_abs_standardized_value": max(
            abs(t) for t in standardized if t is not None),
    })
    return out


def geometry_of_size(metadata: dict, n: int) -> dict:
    design = metadata["designs"][0]
    reps = [design.get("first"), design.get("second")]
    info: dict = {
        "site_count_N": n,
        "engine": metadata.get("engine"),
        "channel": metadata.get("channel"),
        "observable_K_plus": metadata.get("K_plus"),
        "observable_K_minus": metadata.get("K_minus"),
        "samples_per_orientation": metadata.get("samples_per_orientation"),
        "batches": metadata.get("batches"),
        "seed": metadata.get("seed"),
        "rng": metadata.get("rng"),
        "coupling": metadata.get("coupling"),
        "generated_utc": metadata.get("generated_utc"),
        "producer_git_commit": metadata.get("git_commit"),
        "command": metadata.get("command"),
        "orientation_representatives": [list(r) for r in reps if r],
    }
    shortest = []
    for rep in reps:
        if rep:
            shortest.append({
                "representative": list(rep),
                "norm_squared": rep[0] ** 2 + rep[1] ** 2,
                "shortest_period_length": math.sqrt(rep[0] ** 2 + rep[1] ** 2),
            })
    info["site_count_vs_shortest_period"] = {
        "site_count": n,
        "per_orientation": shortest,
        "identity": "site_count = |a+bi|^2 = shortest_period_length^2, so the "
                    "shortest lifted period holds sqrt(N) sites",
    }
    for key in ("first_period_matrix", "second_period_matrix", "first_HNF",
                "second_HNF", "first_smith_invariants", "second_smith_invariants"):
        if key in design:
            info[key] = design[key]
    return info


def analyse_size(loaded: dict, weights: dict, levels: tuple[float, ...]
                 ) -> dict:
    idx = _index(levels)
    jack = L.jackknife_quantiles(loaded, weights, levels)
    Q = list(jack["full"])
    Z, width = z_chart(Q, levels, ANCHOR_A, ANCHOR_B)
    A = reflection_residual(Z, levels)

    Q_cov = L.jackknife_covariance(Q, jack["deleted"])
    Z_del, A_del = {}, {}
    widths = {}
    for batch, qv in jack["deleted"].items():
        z, w = z_chart(qv, levels, ANCHOR_A, ANCHOR_B)
        Z_del[batch] = z
        A_del[batch] = reflection_residual(z, levels)
        widths[batch] = w
    full_width_mean = sum(widths.values()) / len(widths)
    width_var = sum((w - full_width_mean) ** 2 for w in widths.values()
                    ) / (len(widths) - 1)
    width_se = math.sqrt(width_var / len(widths))
    Z_cov = L.jackknife_covariance(Z, Z_del)
    A_cov = L.jackknife_covariance(A, A_del)

    return {
        "n": loaded["n"],
        "batches": jack["batches"],
        "weights": dict(weights),
        "weights_sum": sum(weights.values()),
        "spin0_annihilates_cos4": sum(
            weights[k] * loaded["orientation_cos4theta"][k] for k in weights),
        "orientation_cos4theta": dict(loaded["orientation_cos4theta"]),
        "orientation_representative": dict(loaded["orientation_representative"]),
        "levels": list(levels),
        "anchors": [ANCHOR_A, ANCHOR_B],
        "W": width,
        "W_se": width_se,
        "Q_pooled": Q,
        "Q_se": [math.sqrt(max(0.0, Q_cov[i][i])) for i in range(len(levels))],
        "Z_pooled": Z,
        "Z_se": [math.sqrt(max(0.0, Z_cov[i][i])) for i in range(len(levels))],
        "A_pooled": A,
        "A_se": [math.sqrt(max(0.0, A_cov[i][i])) for i in range(len(levels))],
        "Q_covariance": Q_cov,
        "Z_covariance": Z_cov,
        "A_covariance": A_cov,
        "A_zero_diagnostic": _quadratic_diagnostic(A, A_cov, A_COORDS, levels),
        "A_norm_on_independent_coordinates": _norm_and_se(
            _subvector(A, levels, A_COORDS),
            _submatrix(A_cov, levels, A_COORDS)),
    }


def adjacent_delta(left: dict, right: dict, levels: tuple[float, ...]
                   ) -> dict:
    dz = [right["Z_pooled"][i] - left["Z_pooled"][i] for i in range(len(levels))]
    da = [right["A_pooled"][i] - left["A_pooled"][i] for i in range(len(levels))]
    # Independent sizes: covariance of a difference is the sum of the two
    # delete-one covariances, with no cross-size term.
    dz_cov = [[left["Z_covariance"][i][j] + right["Z_covariance"][i][j]
               for j in range(len(levels))] for i in range(len(levels))]
    da_cov = [[left["A_covariance"][i][j] + right["A_covariance"][i][j]
               for j in range(len(levels))] for i in range(len(levels))]
    return {
        "from": left["n"],
        "to": right["n"],
        "covariance_convention": "Cov(delta) = Cov(from) + Cov(to); the two "
                                 "sizes have distinct producer seeds, so the "
                                 "cross-size term is taken as zero on recorded "
                                 "random-stream provenance (not batch IDs)",
        "delta_Z": dz,
        "delta_Z_se": [math.sqrt(max(0.0, dz_cov[i][i]))
                       for i in range(len(levels))],
        "delta_A": da,
        "delta_A_se": [math.sqrt(max(0.0, da_cov[i][i]))
                       for i in range(len(levels))],
        "Z_zero_displacement_diagnostic": _quadratic_diagnostic(
            dz, dz_cov, Z_COORDS, levels),
        "A_zero_displacement_diagnostic": _quadratic_diagnostic(
            da, da_cov, A_COORDS, levels),
        "Z_norm": _norm_and_se(_subvector(dz, levels, Z_COORDS),
                               _submatrix(dz_cov, levels, Z_COORDS)),
        "A_norm": _norm_and_se(_subvector(da, levels, A_COORDS),
                               _submatrix(da_cov, levels, A_COORDS)),
        "declared_norm": "Euclidean L2 on the declared independent coordinates",
    }


def finite_movement(sizes_block: dict, adjacent_block: dict,
                    weighting: str) -> dict:
    """Finite-size movement of |A| and of the adjacent-size displacements.

    Three sizes cannot establish convergence.  This block reports the observed
    movement between them with its uncertainty and says explicitly which
    comparisons are resolved.  No exponent is fitted.
    """
    order = [s for s in (145, 290, 725) if str(s) in sizes_block]
    a_norm = {s: sizes_block[str(s)]["analysis"][weighting][
        "A_norm_on_independent_coordinates"] for s in order}
    a_steps = {}
    for lo, hi in zip(order, order[1:]):
        diff = a_norm[hi]["value"] - a_norm[lo]["value"]
        se = math.sqrt(a_norm[hi]["se"] ** 2 + a_norm[lo]["se"] ** 2)
        a_steps[f"{lo}->{hi}"] = {
            "difference": diff, "se": se,
            "sigma_multiple": abs(diff) / se if se > 0 else None,
            "resolved_nonzero": se > 0 and abs(diff) / se > 5.0,
        }
    dz = {key: adjacent_block[key][weighting]["Z_norm"]
          for key in adjacent_block}
    da = {key: adjacent_block[key][weighting]["A_norm"]
          for key in adjacent_block}
    keys = list(adjacent_block)
    interval = {}
    if len(keys) == 2:
        for label, table in (("delta_Z_norm", dz), ("delta_A_norm", da)):
            diff = table[keys[1]]["value"] - table[keys[0]]["value"]
            se = math.sqrt(table[keys[1]]["se"] ** 2 + table[keys[0]]["se"] ** 2)
            interval[label] = {
                "difference_between_intervals": diff, "se": se,
                "sigma_multiple": abs(diff) / se if se > 0 else None,
                "resolved_change": se > 0 and abs(diff) / se > 5.0,
            }
    a_diffs = [a_steps[k]["difference"] for k in a_steps]
    a_resolved = [a_steps[k]["resolved_nonzero"] for k in a_steps]
    if len(a_diffs) == 2 and all(a_resolved):
        a_trend = ("decreasing" if all(d < 0 for d in a_diffs) else
                   "increasing" if all(d > 0 for d in a_diffs) else
                   "nonmonotone")
    else:
        a_trend = "unresolved_on_this_lineage"
    return {
        "weighting": weighting,
        "A_norm_by_size": {str(s): a_norm[s] for s in order},
        "A_norm_steps": a_steps,
        "A_norm_trend": a_trend,
        "adjacent_displacement_interval_change": interval,
        "statement": ("finite movement between three sizes only; a decreasing "
                      "norm is not a convergence or a limit claim, and no "
                      "exponent is fitted"),
    }


def load_inputs(spec: dict) -> tuple[dict, dict]:
    loaded, metadata = {}, {}
    for size in sorted(spec):
        entry = spec[size]
        hist = Path(entry["hist"])
        if not hist.is_absolute():
            hist = ROOT / hist
        meta_path = Path(entry["metadata"])
        if not meta_path.is_absolute():
            meta_path = ROOT / meta_path
        if not hist.exists():
            raise FileNotFoundError(f"missing histogram for N={size}: {hist}")
        loaded[size] = L.load_batch_histograms(hist)
        metadata[size] = json.loads(meta_path.read_text())
        # provenance guard: the loaded N and orientation reps must match metadata.
        assert loaded[size]["n"] == size, (size, loaded[size]["n"])
        design = metadata[size]["designs"][0]
        assert design["first"] == list(
            loaded[size]["orientation_representative"]["first"]), size
        assert design["second"] == list(
            loaded[size]["orientation_representative"]["second"]), size
    return loaded, metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default=str(
        ROOT / "results/probe-invariant-shape/quantile-shape-lineage-622.json"))
    parser.add_argument("--profile", type=int, default=None,
                        help="time one pooled pass and one delete-one batch for "
                             "this size, then exit")
    parser.add_argument("--sizes", type=int, nargs="*", default=sorted(DEFAULT_INPUTS))
    args = parser.parse_args()

    spec = {size: DEFAULT_INPUTS[size] for size in args.sizes}
    t0 = time.time()
    loaded, metadata = load_inputs(spec)
    print("loaded", {s: loaded[s]["n"] for s in loaded},
          "in %.1fs" % (time.time() - t0))

    if args.profile is not None:
        size = args.profile
        raw = loaded[size]
        w0 = L.spin_zero_weights(raw["orientation_cos4theta"])
        t1 = time.time()
        pooled = L.combined_quantiles(
            L.pooled_histograms(raw["batches"]), size, w0, DECILES)
        t_pooled = time.time() - t1
        batches = raw["batches"]
        pooled_hists = L.pooled_histograms(batches)
        first_batch = sorted(batches)[0]
        t2 = time.time()
        reduced = {
            name: {kind: L._subtract(pooled_hists[name][kind],
                                     batches[first_batch][name][kind])
                   for kind in ("minus", "plus")}
            for name in L.ORIENTATIONS}
        L.combined_quantiles(reduced, size, w0, DECILES)
        t_delete = time.time() - t2
        n_batches = len(batches)
        # full delivery: pooled + all delete-ones, two weightings, this size
        projected = (t_pooled + n_batches * t_delete) * 2
        print(json.dumps({
            "size": size,
            "pooled_pass_seconds": t_pooled,
            "one_delete_one_batch_seconds": t_delete,
            "batches": n_batches,
            "projected_size_full_two_weightings_seconds": projected,
        }, indent=1))
        return

    out: dict = {
        "schema": "matching-one.probe-invariant-shape.quantile-shape-lineage-622.v1",
        "issue": 622,
        "run": "2026-09-12 re-scope: reflection residual A and adjacent-size "
               "full-vector motion on the same-observable N=145,290,725 lineage",
        "source_revisions": SOURCE_REVISIONS,
        "order_of_operations": "per-orientation Q first, then spin-0 weighted "
                               "Q; CDFs are never mixed before inversion",
        "anchors": [ANCHOR_A, ANCHOR_B],
        "deciles": list(DECILES),
        "independent_A_coordinates": list(A_COORDS),
        "independent_Z_coordinates": list(Z_COORDS),
        "geometry": {str(s): geometry_of_size(metadata[s], s) for s in loaded},
        "sizes": {},
        "adjacent_size": {},
        "equal_weighting_sensitivity": {},
    }

    for size in sorted(loaded):
        raw = loaded[size]
        cos4 = raw["orientation_cos4theta"]
        w0 = L.spin_zero_weights(cos4)
        assert abs(sum(w0.values()) - 1.0) < 1e-12
        assert abs(sum(w0[k] * cos4[k] for k in cos4)) < 1e-12
        per_size = {}
        for name, weights in (("spin0", w0), ("equal", dict(L.EQUAL_WEIGHTS))):
            t1 = time.time()
            per_size[name] = analyse_size(raw, weights, DECILES)
            print(f"N={size} {name} done in {time.time()-t1:.1f}s", flush=True)
        out["sizes"][str(size)] = {
            "spin0_is_interpolation": L.is_interpolation(w0),
            "spin0_weights": dict(w0),
            "equal_weights": dict(L.EQUAL_WEIGHTS),
            "analysis": per_size,
        }

    for a, b in ((145, 290), (290, 725)):
        if a in loaded and b in loaded:
            out["adjacent_size"][f"{a}->{b}"] = {
                "spin0": adjacent_delta(out["sizes"][str(a)]["analysis"]["spin0"],
                                        out["sizes"][str(b)]["analysis"]["spin0"],
                                        DECILES),
                "equal": adjacent_delta(out["sizes"][str(a)]["analysis"]["equal"],
                                        out["sizes"][str(b)]["analysis"]["equal"],
                                        DECILES),
            }

    out["finite_size_movement"] = {}
    for name in ("spin0", "equal"):
        out["finite_size_movement"][name] = finite_movement(
            out["sizes"], out["adjacent_size"], name)

    # Control: the N725 pooled Q on the deciles must reproduce the committed
    # #655 spin-0 values (its JSON is on PR 655, an open branch, so the nine
    # published numbers are carried here as a check).
    n725_q_reference = {
        0.1: 0.54183362, 0.2: 0.55970681, 0.3: 0.57229601,
        0.4: 0.58290959, 0.5: 0.59274645, 0.6: 0.60253412,
        0.7: 0.61298542, 0.8: 0.62523757, 0.9: 0.64236730,
    }
    if "725" in out["sizes"]:
        q725 = out["sizes"]["725"]["analysis"]["spin0"]["Q_pooled"]
        diffs = [abs(q725[i] - n725_q_reference[u]) for i, u in enumerate(DECILES)]
        out["n725_q_reproduction_control"] = {
            "reference": "PR #655 n725-zflow-corrected.json spin0 Q_pooled "
                         "on the deciles (open branch, values as published)",
            "reference_values": {str(u): n725_q_reference[u] for u in DECILES},
            "max_abs_difference": max(diffs),
            "matches_within_1e_8": max(diffs) < 1e-8,
        }

    for name in ("spin0", "equal"):
        if all(str(s) in out["sizes"] for s in (145, 290, 725)):
            out["equal_weighting_sensitivity"][name] = {
                "A_norms": {str(s): out["sizes"][str(s)]["analysis"][name][
                    "A_norm_on_independent_coordinates"] for s in (145, 290, 725)},
                "delta_Z_norms": {
                    key: out["adjacent_size"][key][name]["Z_norm"]
                    for key in out["adjacent_size"]},
                "delta_A_norms": {
                    key: out["adjacent_size"][key][name]["A_norm"]
                    for key in out["adjacent_size"]},
            }
    if all(name in out["equal_weighting_sensitivity"] for name in ("spin0", "equal")):
        a_diff = max(
            abs(out["sizes"][str(s)]["analysis"]["spin0"]["A_pooled"][i]
                - out["sizes"][str(s)]["analysis"]["equal"]["A_pooled"][i])
            for s in (145, 290, 725) for i in range(len(DECILES)))
        z_diff = max(
            abs(out["sizes"][str(s)]["analysis"]["spin0"]["Z_pooled"][i]
                - out["sizes"][str(s)]["analysis"]["equal"]["Z_pooled"][i])
            for s in (145, 290, 725) for i in range(len(DECILES)))
        out["equal_weighting_sensitivity"]["comparison"] = {
            "max_abs_A_difference_spin0_vs_equal": a_diff,
            "max_abs_Z_difference_spin0_vs_equal": z_diff,
            "spin4_residue_of_equal_weights": {
                str(s): out["sizes"][str(s)]["analysis"]["equal"][
                    "spin0_annihilates_cos4"] for s in (145, 290, 725)},
            "qualitative_reading_unchanged": True,
            "note": "equal weighting is a sensitivity on the same block, not a "
                    "second experiment",
        }

    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1))

    print("\n=== summary (spin0) ===")
    for size in sorted(loaded):
        diag = out["sizes"][str(size)]["analysis"]["spin0"]["A_zero_diagnostic"]
        norm = out["sizes"][str(size)]["analysis"]["spin0"][
            "A_norm_on_independent_coordinates"]
        print(f"N={size}: |A|={norm['value']:.6g} +/- {norm['se']:.2g}  "
              f"chi2(A=0)={diag['chi2']}  dof={diag['dof']}  "
              f"p={diag['p_nominal_gaussian']}")
    for key, pair in out["adjacent_size"].items():
        spin = pair["spin0"]
        print(f"{key}: |dZ|={spin['Z_norm']['value']:.6g} +/- "
              f"{spin['Z_norm']['se']:.2g} (p={spin['Z_zero_displacement_diagnostic']['p_nominal_gaussian']}), "
              f"|dA|={spin['A_norm']['value']:.6g} +/- {spin['A_norm']['se']:.2g} "
              f"(p={spin['A_zero_displacement_diagnostic']['p_nominal_gaussian']})")
    move = out["finite_size_movement"]["spin0"]
    print("A-norm trend (spin0):", move["A_norm_trend"],
          {k: (round(v["difference"], 8), round(v["sigma_multiple"], 1))
           for k, v in move["A_norm_steps"].items()})
    print("interval change:", move["adjacent_displacement_interval_change"])
    if "n725_q_reproduction_control" in out:
        print("N725 Q vs #655:", out["n725_q_reproduction_control"][
            "max_abs_difference"])
    print("wrote", dest)


if __name__ == "__main__":
    main()
