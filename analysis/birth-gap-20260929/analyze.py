#!/usr/bin/env python3
"""Joint-birth finite-size readout; stdlib, no fitting and no new random samples.

Conditional on D=J2-J1, the uniform-label gap has Beta(D,N+1-D)
law (a point mass at zero for D=0). Integer-parameter CDFs are evaluated
by the exact binomial-tail identity, using floating-point recurrences.
"""
import argparse
import collections
import gzip
import json
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
EPSILONS = (0.025, 0.05, 0.1, 0.2)


def read_batch(path):
    with gzip.open(path, "rt") as f:
        return json.load(f)


def combine(batches):
    hist = collections.Counter()
    for b in batches:
        for j1, j2, count in b["histogram"]:
            hist[j1, j2] += count
    return hist


def quantile(hist, u):
    target = u * sum(hist.values())
    running = 0
    for value, count in sorted(hist.items()):
        running += count
        if running >= target:
            return value
    raise ValueError("empty quantile histogram")


def beta_gap_cdfs(n, x):
    """Return P[Beta(d,n+1-d)<=x]; omitted right tail underflows to zero.

    Bin(n,x) survival, accumulated backwards to avoid cancellation near zero.
    This recurrence is not a general extreme-tail special-function library.
    Refuse numerical underflow at start. Terminating the already-underflowed
    right tail makes L=512 inexpensive without a scientific approximation.
    """
    if x <= 0:
        return [1.0] + [0.0] * n
    if x >= 1:
        return [1.0] * (n + 1)
    pmf = [math.exp(n * math.log1p(-x))]
    if pmf[0] == 0:
        raise ValueError("binomial recurrence start underflow; out of pilot range")
    ratio = x / (1 - x)
    for k in range(n):
        next_value = pmf[k] * ((n - k) / (k + 1)) * ratio
        if next_value == 0:
            break
        pmf.append(next_value)
    total = math.fsum(pmf)
    if abs(total - 1) > 1e-9:
        raise ArithmeticError(f"binomial recurrence normalization {total}")
    survival = [0.0] * len(pmf)
    running = 0.0
    for k in range(len(pmf) - 1, 0, -1):
        running += pmf[k]
        survival[k] = running / total
    survival[0] = 1.0
    return survival


def metrics(hist, L):
    n, count = L * L, sum(hist.values())
    pooled = collections.Counter()
    gaps = collections.Counter()
    s1 = s2 = s11 = s22 = s12 = 0
    for (j1, j2), c in hist.items():
        pooled[j1] += c
        pooled[j2] += c
        gaps[j2 - j1] += c
        s1 += c * j1; s2 += c * j2
        s11 += c * j1 * j1; s22 += c * j2 * j2; s12 += c * j1 * j2
    width = quantile(pooled, .75) - quantile(pooled, .25)
    if width <= 0:
        raise ValueError("birth-mixture IQR is zero")
    scale = L ** 1.25
    mean_d = sum(d * c for d, c in gaps.items()) / count
    atom = gaps[0] / count
    var1 = s11 / count - (s1 / count) ** 2
    var2 = s22 / count - (s2 / count) ** 2
    cov12 = s12 / count - s1 * s2 / count ** 2
    result = {
        "direct_atom": atom,
        "mean_D": mean_d,
        "mean_D_over_L54": mean_d / scale,
        "pooled_birth_IQR_J": float(width),
        "pooled_birth_IQR_J_over_L54": width / scale,
        "mean_D_over_IQR": mean_d / width,
        "mean_label_gap": mean_d / (n + 1),
        "birth_J_correlation": cov12 / math.sqrt(var1 * var2),
    }
    for eps in EPSILONS:
        tag = f"{eps:g}"
        result[f"P_D_L54_le_{tag}"] = sum(c for d, c in gaps.items() if d <= eps * scale) / count
        result[f"P_D_IQR_le_{tag}"] = sum(c for d, c in gaps.items() if d <= eps * width) / count
        result[f"Z_D_L54_{tag}"] = math.fsum(c * max(0.0, 1 - d / (eps * scale)) for d, c in gaps.items()) / count
        result[f"Z_D_IQR_{tag}"] = math.fsum(c * max(0.0, 1 - d / (eps * width)) for d, c in gaps.items()) / count
        for key, x in (("label_L34", eps / L ** .75),
                       ("label_IQR", eps * width / (n + 1))):
            probabilities = beta_gap_cdfs(n, x)
            result[f"P_{key}_le_{tag}"] = math.fsum(c * probabilities[d] for d, c in gaps.items() if d < len(probabilities)) / count
            tilted = beta_gap_cdfs(n + 1, x)
            # E[(1-B/x)+], B~Beta(d,n+1-d), exact truncated first
            # moment. Tilted beta has parameters d+1,n+1-d (sum n+2).
            result[f"Z_{key}_{tag}"] = math.fsum(
                c * max(0.0, (probabilities[d] if d < len(probabilities) else 0.0)
                        - d / ((n + 1) * x) * (tilted[d + 1] if d + 1 < len(tilted) else 0.0))
                for d, c in gaps.items()) / count
    return result


def wilson(count, n):
    z = 1.959963984540054
    p = count / n
    center = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [max(0.0, center - half), min(1.0, center + half)]


def summarize(batch_data, L):
    hist = combine(batch_data)
    full = metrics(hist, L)
    names = list(full)
    b = len(batch_data)
    leave = [metrics(combine(batch_data[:k] + batch_data[k+1:]), L) for k in range(b)]
    means = {key: sum(x[key] for x in leave) / b for key in names}
    covariance = [[(b - 1) / b * sum((x[a] - means[a]) * (x[c] - means[c]) for x in leave)
                   for c in names] for a in names]
    count = sum(hist.values())
    atom = sum(c for (j1, j2), c in hist.items() if j1 == j2)
    return {
        "L": L, "N": L*L, "samples": count, "batch_count": b,
        "point": full,
        "jackknife_se": {key: math.sqrt(max(0, covariance[i][i])) for i, key in enumerate(names)},
        "direct_atom_count": atom,
        "direct_atom_wilson_two_sided_95": wilson(atom, count),
        "direct_atom_zero_count_one_sided_95_upper": -math.expm1(math.log(.05)/count) if atom == 0 else None,
        "metric_order": names,
        "delete_one_batch_covariance": covariance,
        "delete_one_batch_metrics": leave,
        "uncertainty": f"paired complete-batch delete-one; re-estimates mixture IQR in each delete-one; {b} batches, exploratory uncertainty only",
        "positive_near_diagonal_rule": "subtract direct_atom from any P_* (zero is otherwise included)",
    }


def main():
    global EPSILONS
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=pathlib.Path, default=HERE / "data")
    ap.add_argument("--output-json", type=pathlib.Path, default=HERE / "summary.json")
    ap.add_argument("--output-md", type=pathlib.Path, default=HERE / "SUMMARY.md")
    ap.add_argument("--epsilons", type=float, nargs="+", default=list(EPSILONS),
                    help="fixed dimensionless CDF/Z grid, choose before reading the new block")
    args = ap.parse_args()
    if any(not math.isfinite(e) or e <= 0 for e in args.epsilons):
        ap.error("epsilons must be finite and positive")
    EPSILONS = tuple(sorted(set(args.epsilons)))
    manifest = json.loads((args.data / "run.json").read_text())
    if manifest["status"] != "completed":
        raise ValueError("run is incomplete")
    groups = collections.defaultdict(list)
    controls = []
    for ref in manifest["controls"]:
        b = read_batch(args.data / ref["file"])
        controls.append({"lattice": b["lattice"], "L": 3, "permutations": b["samples"],
                         "direct_atom_count": b["direct_atom_count"],
                         "mean_D": b["sum_D"] / b["samples"],
                         "mean_D2": b["sum_D2"] / b["samples"]})
    for ref in manifest["batches"]:
        b = read_batch(args.data / ref["file"])
        groups[b["lattice"], b["L"]].append(b)
    cells = []
    for (lattice, L), batches in sorted(groups.items()):
        cell = summarize(batches, L)
        cell["lattice"] = lattice
        cells.append(cell)
    report = {
        "schema": "matching-one.joint-birth-pilot.v1",
        "epsilons": EPSILONS,
        "definitions": {
            "J1_J2": "first occupied-site counts at global ambient homology rank >=1 and =2, SAME uniform permutation",
            "D": "J2-J1",
            "L54": "D/L^(5/4), prescribed nearcritical benchmark, no exponent fit; not proved square-site scaling",
            "label_L34": "L^(3/4)*(T2-T1), conditional Beta(D,N+1-D) mixing of uniform labels",
            "pooled_birth_IQR_J": "Q75-Q25 of equal-weight mixture of J1 and J2; empirical lower quantiles",
            "label_IQR": "(N+1)*(T2-T1)/(pooled_birth_IQR_J); exponent-free self-normalized readout",
            "Z": "E[(1-G/delta)+] with delta on the same epsilon grid; exact histogram or conditional-Beta integration, added on theory-agent request after acquisition",
            "triangular_geometry": "nearest neighbours (+/-1,0),(0,+/-1),(1,1),(-1,-1); LxL primitive-coordinate rhombus, not square conformal modulus",
        },
        "controls": controls, "exact_control_reference": manifest.get("exact_control_reference"), "cells": cells,
        "production_role": manifest.get("production_role", "pilot"),
        "total_production_samples": sum(c["samples"] for c in cells),
        "production_seconds_sum": manifest["aggregate_engine_seconds"],
        "acquisition_wall_seconds_including_build_control_benchmark": manifest["wall_seconds"],
        "limitations": [
            "Finite-size pilot, not a noncoalescence proof or arm-exponent estimate.",
            "A finite epsilon grid does not implement epsilon->0 after L->infinity.",
            "No cross-lattice equality asserted: primitive square and triangular tori have different conformal moduli.",
            "All readouts within one cell use the SAME sample block; Beta reconstruction adds no independent evidence.",
            "Finite batch counts support exploratory uncertainty only; no p-values, model selection or data-chosen exponent.",
            "A stable typical gap does not exclude an additional limiting diagonal atom.",
        ]}
    args.output_json.write_text(json.dumps(report, indent=2) + "\n")
    batch_counts = sorted(set(c["batch_count"] for c in cells))
    sample_counts = sorted(set(c["samples"] for c in cells))
    platform_label = f"{manifest['machine']} / Python {manifest['python']}"
    lines = ["# Executed joint-birth block — 2026-09-29", "",
             f"{report['total_production_samples']:,} independent-filtration samples; per-cell batch counts {batch_counts}, sample counts {sample_counts}.",
             f"Acquisition ({platform_label}) including compilation, requested controls and benchmarks: {manifest['wall_seconds']:.2f} s, at most {manifest['workers']} concurrent jobs.", "",
             f"Production seed namespace: `{report['production_role']}`.", "",
             "## Scale and direct jump", "",
             "Errors after ± are complete-batch jackknife SE, not confidence limits.", "",
             "| Lattice | L | Direct 0→2 / samples | E[D]/L^(5/4) | E[D]/mixture IQR |", "|---|---:|---:|---:|---:|"]
    for c in cells:
        p, e = c["point"], c["jackknife_se"]
        lines.append(f"| {c['lattice']} | {c['L']} | {c['direct_atom_count']} / {c['samples']} | {p['mean_D_over_L54']:.4f} ± {e['mean_D_over_L54']:.4f} | {p['mean_D_over_IQR']:.4f} ± {e['mean_D_over_IQR']:.4f} |")
    grid_header = "| Lattice | L | " + " | ".join(f"ε/δ={eps:g}" for eps in EPSILONS) + " |"
    grid_separator = "|---|---:|" + "---:|" * len(EPSILONS)
    for key, title in (("label_L34", "Continuous-label near diagonal: P[L^(3/4)(T2−T1) ≤ ε]"),
                       ("label_IQR", "Exponent-free near diagonal: P[(N+1)(T2−T1)/W ≤ ε]"),
                       ("D_L54", "Discrete near diagonal: P[D/L^(5/4) ≤ ε]"),
                       ("D_IQR", "Discrete exponent-free near diagonal: P[D/W ≤ ε]")):
        lines += ["", "## " + title, "", "W is the pooled J1/J2 interquartile width; zero atom included. Subtract the preceding direct-atom frequency for strictly positive mass.", "",
                  grid_header, grid_separator]
        for c in cells:
            values = [f"{c['point'][f'P_{key}_le_{eps:g}']:.5f} ± {c['jackknife_se'][f'P_{key}_le_{eps:g}']:.5f}" for eps in EPSILONS]
            lines.append(f"| {c['lattice']} | {c['L']} | " + " | ".join(values) + " |")
    for key, title in (("label_L34", "Soft near-diagonal Z(δ), G=L^(3/4)(T2−T1)"),
                       ("label_IQR", "Soft exponent-free Z(δ), G=(N+1)(T2−T1)/W")):
        lines += ["", "## " + title, "", "Z(δ)=E[(1−G/δ)+]; same samples and covariance, not additional independent evidence.", "",
                  grid_header, grid_separator]
        for c in cells:
            values = [f"{c['point'][f'Z_{key}_{eps:g}']:.5f} ± {c['jackknife_se'][f'Z_{key}_{eps:g}']:.5f}" for eps in EPSILONS]
            lines.append(f"| {c['lattice']} | {c['L']} | " + " | ".join(values) + " |")
    control_text = ("L=3 was enumerated once per lattice, all 9! permutations. Square: P(D=0)=3/35, E[D]=3/2, E[D²]=43/14. Triangular: 2/35, 3/2, 81/28. Every exact identity passed before this block."
                    if controls else "This block did not repeat exact enumeration; supplied prior-control reference: " + str(manifest.get("exact_control_reference")))
    lines += ["", "## Exact control and uncertainty", "", control_text, "",
              "JSON preserves full delete-one covariance, all leave-one estimates, atom Wilson intervals and an exact one-sided 95% upper bound if an atom count is zero. Zero observations give only the finite-sample upper bound 1−0.05^(1/n), never zero probability.", "",
              "The Beta CDF uses P[Beta(d,N+1−d)≤x]=P[Bin(N,x)≥d]. Backward binomial-tail sums avoid subtraction near zero; d=0 is retained as an atom. This integrates over labels without adding Monte Carlo noise or evidence.", "",
              "## What this does and does not decide", "",
              "Read the accompanying research note for the interpretation. This pilot distinguishes shrinking single-step double births from finite near-diagonal occupancy. It does not establish absence of a diagonal atom in any scaling limit.", ""]
    args.output_md.write_text("\n".join(lines))
    print(args.output_md.read_text())


if __name__ == "__main__":
    main()
