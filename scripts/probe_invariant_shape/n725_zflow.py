#!/usr/bin/env python3
"""Per-batch anchored-shape flow of the N=725 threshold law, for #622.

Reads (read-only) the committed N=725 full-curve histogram of PR #614
(``results/server-20260907/P612-n725-fullcurve/raw/n725_100m.hist.csv``,
pulled to a temp path via ``git show`` — nothing in the source tree is
touched) and reconstructs, per batch and per orientation combination:

    F_batch(p) = sum_j B(N, j, p) G_batch(j),   G_batch = P(K <= j)

exactly as ``scripts/threshold_quantile_lineage.py`` defines it (same binomial
window, same rank_cdf mixture weights), and extracts

    Q_batch(u)  on the frozen decile grid (u = 0.1, ..., 0.9)
    Z_batch     with anchors a = 1/4, b = 3/4 (added to the grid)

The scientific outputs:

  * mean Z and its per-batch standard error at N=725 — is the shape pinned?
  * the three per-batch spacings Z(0.5)-Z(0.3), Z(0.7)-Z(0.5), Z(0.9)-Z(0.7)
    and their SEs — three nearly-independent shape readings;
  * the batch-level correlations between the spacings — a common mode across
    spacings is the signature of a *shape* degree of freedom (as opposed to
    independent level noise), and it is what a limit object S requires to
    even be measurable;
  * per-level batch SEs of Q (already committed by #612 as ~3e-4 at N=725;
    re-derived here in passing for self-consistency, not as a deliverable).

Conventions kept: ``equal`` and ``spin0`` orientation weightings, exactly as
in the lineage module; the spin0 combination is the primary (declared in
#582's weighting systematic), equal is the sensitivity.

No Monte Carlo is started; no block is created; this is a read of one
committed histogram plus exact algebra.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

N725 = 725
A, B = 0.25, 0.75
LEVELS = (0.1, 0.2, 0.3, A, 0.5, B, 0.7, 0.8, 0.9)
LOG_WEIGHT_FLOOR = 700.0
QUANTILE_TOLERANCE = 1e-12
HIST = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/n725_100m.hist.csv")


def load_batch_histograms(path: Path) -> dict:
    counts = defaultdict(lambda: defaultdict(lambda: {"minus": defaultdict(int),
                                                      "plus": defaultdict(int)}))
    reps = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            orientation = row["orientation"]
            rep = (int(row["a"]), int(row["b"]))
            reps[orientation] = rep
            counts[int(row["batch"])][orientation][row["kind"]][int(row["k"])] += \
                int(row["count"])
    assert set(reps) == {"first", "second"}
    return {
        "n": N725,
        "orientation_representative": reps,
        "batches": {b: {o: {"minus": dict(p[o]["minus"]),
                            "plus": dict(p[o]["plus"])}
                        for o in reps}
                    for b, p in sorted(counts.items())},
    }


def cos_four_theta(a: int, b: int) -> float:
    n = a * a + b * b
    return float(a * a - b * b) / n


def rank_cdf(minus, plus, n: int) -> list[float]:
    minus_total = sum(minus.values())
    plus_total = sum(plus.values())
    assert minus_total > 0 and plus_total > 0
    assert minus_total == plus_total
    denom = 2.0 * minus_total
    cum, out = 0.0, [0.0]
    for rank in range(1, n + 1):
        cum += (minus.get(rank, 0) + plus.get(rank, 0)) / denom
        out.append(cum)
    assert abs(out[-1] - 1.0) < 1e-9
    return out


def binom_weights(n: int, p: float):
    lp, lq = math.log(p), math.log1p(-p)
    lfn = math.lgamma(n + 1)

    def log_pmf(j):
        return lfn - math.lgamma(j + 1) - math.lgamma(n - j + 1) \
            + j * lp + (n - j) * lq

    mode = min(n, max(0, int((n + 1) * p)))
    peak = log_pmf(mode)
    low = mode
    while low > 0 and peak - log_pmf(low - 1) < LOG_WEIGHT_FLOOR:
        low -= 1
    high = mode
    while high < n and peak - log_pmf(high + 1) < LOG_WEIGHT_FLOOR:
        high += 1
    return low, [math.exp(log_pmf(j) - peak) for j in range(low, high + 1)]


def profile_cdf(g: list[float], n: int, p: float) -> float:
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    low, w = binom_weights(n, p)
    total = 0.0
    mass = 0.0
    for off, wt in enumerate(w):
        mass += wt
        total += wt * g[low + off]
    return total / mass


def quantile(g: list[float], n: int, level: float) -> float:
    lo, hi = 0.0, 1.0
    while hi - lo > QUANTILE_TOLERANCE:
        mid = 0.5 * (lo + hi)
        if profile_cdf(g, n, mid) < level:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def combine(minus_first, plus_first, minus_second, plus_second, n, w1, w2):
    """Weighted mixture of the two orientations' rank CDFs (weights on the
    *histograms*, matching the lineage module's combination semantics:
    spin0 weights are on the orientation laws, applied at the histogram
    level before the common binomial pass)."""
    g1 = rank_cdf(minus_first, plus_first, n)
    g2 = rank_cdf(minus_second, plus_second, n)
    return [w1 * a + w2 * b for a, b in zip(g1, g2)]


def main() -> None:
    loaded = load_batch_histograms(HIST)
    reps = loaded["orientation_representative"]
    c1 = cos_four_theta(*reps["first"])
    c2 = cos_four_theta(*reps["second"])
    gap = c1 - c2
    w_first_spin0 = -c2 / gap
    w_second_spin0 = c1 / gap
    weightings = {
        "spin0": (w_first_spin0, w_second_spin0),
        "equal": (0.5, 0.5),
    }
    n = N725
    out = {
        "schema": "matching-one.probe-invariant-shape.n725-zflow.v1",
        "source": "PR #614 results/server-20260907/P612-n725-fullcurve/raw/"
                  "n725_100m.hist.csv (read via git show; read-only)",
        "n": n,
        "anchors": [A, B],
        "levels_including_anchors": list(LEVELS),
        "orientation_representative": {k: list(v) for k, v in reps.items()},
        "orientation_cos4theta": {"first": c1, "second": c2},
        "weightings": {},
    }
    for name, (w1, w2) in weightings.items():
        qrows = []
        for batch, payload in loaded["batches"].items():
            g = combine(payload["first"]["minus"], payload["first"]["plus"],
                        payload["second"]["minus"], payload["second"]["plus"],
                        n, w1, w2)
            qs = [quantile(g, n, u) for u in LEVELS]
            qa, qb = qs[3], qs[5]
            z = [(q - qa) / (qb - qa) for q in qs]
            qrows.append({"batch": batch, "Q": qs, "Z": z})
        # stats
        def mean(xs):
            return sum(xs) / len(xs)
        def se(xs):
            m = mean(xs)
            var = sum((x - m) ** 2 for x in xs) / (len(xs) - 1)
            return math.sqrt(var / len(xs))
        zcols = list(zip(*[r["Z"] for r in qrows]))
        qcols = list(zip(*[r["Q"] for r in qrows]))
        z_mean = [mean(c) for c in zcols]
        z_se = [se(c) for c in zcols]
        q_mean = [mean(c) for c in qcols]
        q_se = [se(c) for c in qcols]
        # spacings on the Z grid at levels (0.3, 0.5), (0.5, 0.7), (0.7, 0.9)
        # -> indices 2,4 4,6 6,8 in the LEVELS order above
        spac = {f"Z({LEVELS[i]},{LEVELS[j]})":
                [r["Z"][j] - r["Z"][i] for r in qrows]
                for (i, j) in ((2, 4), (4, 6), (6, 8))}
        spac_stats = {k: {"mean": mean(v), "se": se(v)} for k, v in spac.items()}
        # correlations between spacings (batch level)
        keys = list(spac)
        cors = {}
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                x, y = spac[keys[i]], spac[keys[j]]
                mx, my = mean(x), mean(y)
                sx = math.sqrt(sum((a - mx) ** 2 for a in x))
                sy = math.sqrt(sum((b - my) ** 2 for b in y))
                cors[f"{keys[i]} vs {keys[j]}"] = \
                    sum((a - mx) * (b - my) for a, b in zip(x, y)) / (sx * sy)
        # common-mode test: first principal direction of the 3 spacings
        # (power iteration on the 3x3 covariance)
        import itertools as it
        vec = [1.0, 1.0, 1.0]
        cov = [[0.0] * 3 for _ in range(3)]
        cols = [spac[k] for k in keys]
        mus = [mean(c) for c in cols]
        for a_ in range(3):
            for b_ in range(3):
                cov[a_][b_] = sum((x - mus[a_]) * (y - mus[b_])
                                  for x, y in zip(cols[a_], cols[b_])) / (len(cols[0]) - 1)
        for _ in range(200):
            nv = [sum(cov[i][j] * vec[j] for j in range(3)) for i in range(3)]
            norm = math.sqrt(sum(x * x for x in nv))
            vec = [x / norm for x in nv]
        out["weightings"][name] = {  # type: ignore[index]
            "weights": {"first": w1, "second": w2},
            "batches": len(qrows),
            "Q_mean": q_mean, "Q_se": q_se,
            "Z_mean": z_mean, "Z_se": z_se,
            "spacings": spac_stats,
            "spacing_correlations": cors,
            "spacing_pc1_direction": vec,
            "z_per_batch_first3": [r["Z"] for r in qrows[:3]],
        }
    dest = ROOT / "results" / "probe-invariant-shape" / "n725-zflow.json"
    dest.write_text(json.dumps(out, indent=1))
    for name in weightings:
        w = out["weightings"][name]  # type: ignore[index]
        print(f"== {name}")
        print("  Z_mean:", [round(v, 4) for v in w["Z_mean"]])
        print("  Z_se  :", ["%.2e" % v for v in w["Z_se"]])
        print("  spacings:", {k: (round(v['mean'], 4), '%.2e' % v['se'])
                              for k, v in w["spacings"].items()})
        print("  corr:", {k: round(v, 3) for k, v in w["spacing_correlations"].items()})
        print("  pc1  :", [round(v, 3) for v in w["spacing_pc1_direction"]])


if __name__ == "__main__":
    main()
