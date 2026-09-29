#!/usr/bin/env python3
"""One preselected three-time rank-memory contrast on existing paired births.

Retrospective archive analysis; no resampling, exponent fit, cutoff scan,
new percolation simulation, or independent-validation claim. Standard library.
The source cloud block has already been used for gap statistics.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import gzip
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import platform
import sys
import time


HERE = Path(__file__).resolve().parent
QUANTILES = (0.25, 0.50, 0.75)


def lower_quantile(hist, q):
    target = q * sum(hist.values())
    total = 0
    for k, count in sorted(hist.items()):
        total += count
        if total >= target:
            return k
    raise ValueError("empty birth distribution")


def combine(batches):
    hist = Counter()
    for b in batches:
        hist.update({(j1, j2): count for j1, j2, count in b["histogram"]})
    return hist


def table_metrics(table, samples):
    # Rows: early/late first birth. Columns: survives past c / exits by c.
    n11, n10 = table[0]
    n01, n00 = table[1]
    early, late = n11 + n10, n01 + n00
    risk = early + late
    if min(early, late) <= 0:
        raise ValueError("empty early/late rank-one risk group")
    pe, pl = n11 / early, n01 / late
    covariance = (n11 * n00 - n10 * n01) / risk**2
    contrast = pe - pl
    # Same quantity by its conditional-covariance identity.
    assert math.isclose(covariance, early * late / risk**2 * contrast,
                        rel_tol=1e-12, abs_tol=1e-15)
    return {
        "rank_one_risk_fraction": risk / samples,
        "early_fraction_within_risk": early / risk,
        "survival_early": pe,
        "survival_late": pl,
        "survival_difference": contrast,
        "conditional_covariance": covariance,
    }


def summarize_hist(hist):
    marginal = Counter()
    for (j1, j2), count in hist.items():
        marginal[j1] += count
        marginal[j2] += count
    a, b, c = [lower_quantile(marginal, q) for q in QUANTILES]
    if not a < b < c:
        raise ValueError("distinct three-time cutoffs not available")
    table = [[0, 0], [0, 0]]
    for (j1, j2), count in hist.items():
        if j1 <= b < j2:
            table[0 if j1 <= a else 1][0 if j2 > c else 1] += count
    samples = sum(hist.values())
    point = table_metrics(table, samples)
    n11, n10 = table[0]
    n01, n00 = table[1]
    # These are empirical kernel numerators, all over the same sample count.
    kernel_counts = {
        "H_a_b": n11 + n10,
        "H_a_c": n11,
        "H_b_c": n11 + n01,
        "H_b_b": n11 + n10 + n01 + n00,
    }
    return {
        "cutoffs_insertion_counts": [a, b, c],
        "risk_table_early_late_by_survive_exit": table,
        "kernel_counts": kernel_counts,
        "point": point,
    }


def summarize_cell(batches, refs, lattice, L):
    b = len(batches)
    if b < 3 or len({x["samples"] for x in batches}) != 1:
        raise ValueError("this delete-one contract requires equal complete batches")
    full = summarize_hist(combine(batches))
    leave = [summarize_hist(combine(batches[:i] + batches[i+1:]))
             for i in range(b)]
    names = list(full["point"])
    centers = {n: math.fsum(r["point"][n] for r in leave)/b for n in names}
    cov = [[(b-1)/b * math.fsum((r["point"][a]-centers[a]) *
                              (r["point"][c]-centers[c]) for r in leave)
            for c in names] for a in names]
    return {
        "lattice": lattice, "L": L, "N": L*L,
        "samples": sum(x["samples"] for x in batches), "batches": b,
        "role": "primary_largest_available_size" if L == 512 else "size_context",
        **full,
        "metric_order": names,
        "jackknife_se": {n: math.sqrt(max(0, cov[i][i])) for i, n in enumerate(names)},
        "delete_one_batch_covariance": cov,
        "delete_one": [{"omitted_batch": refs[i]["batch"], **r}
                       for i, r in enumerate(leave)],
        "source_batches": [{k: r[k] for k in
                            ("file", "batch", "samples", "sha256", "dependency_group")}
                           for r in refs],
    }


def exact_control():
    # One four-cell independence example and one dependence example;
    # not a physical simulation or a rerun of any earlier census.
    independent, dependent = [[1, 1], [1, 1]], [[3, 1], [1, 3]]
    assert table_metrics(independent, 4)["survival_difference"] == 0
    got = table_metrics(dependent, 8)
    assert Fraction(got["survival_difference"]) == Fraction(1, 2)
    assert Fraction(got["conditional_covariance"]) == Fraction(1, 8)
    # Birth-pair cells exactly realize the four risk-table categories at 2,3,4.
    pairs = {(1, 5): 3, (1, 4): 1, (3, 5): 1, (3, 4): 3}
    table = [[0, 0], [0, 0]]
    for (j1, j2), count in pairs.items():
        if j1 <= 3 < j2:
            table[int(j1 > 2)][int(j2 <= 4)] += count
    assert table == dependent
    return {"independent_table": independent, "dependent_table": dependent,
            "dependent_contrast": "1/2", "dependent_covariance": "1/8",
            "scope": "one exact four-cell arithmetic and inequality-boundary control"}


def exact_archive_analysis(data_path):
    """Analyze all triples of the existing 9! tables; never regenerate them."""
    source = json.loads((data_path / "run.json").read_text())
    results = []
    for ref in source["controls"]:
        raw = (data_path / ref["file"]).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == ref["sha256"]
        b = json.loads(gzip.decompress(raw))
        assert b["L"] == 3 and b["N"] == 9 and b["mode"] == "exact"
        assert sum(row[2] for row in b["histogram"]) == b["samples"] == math.factorial(9)
        def h(s, t):
            return sum(n for j1, j2, n in b["histogram"] if j1 <= s and j2 > t)
        count = 0
        violations = []
        for a, middle, c in combinations(range(b["N"]+1), 3):
            count += 1
            A, B, C, D = h(a, middle), h(a, c), h(middle, c), h(middle, middle)
            determinant = B*D-A*C
            if determinant:
                assert D > A > 0
                violations.append({
                    "cutoffs": [a, middle, c],
                    "risk_table": [[B, A-B], [C-B, D-A-C+B]],
                    "conditional_covariance": str(Fraction(determinant, D**2)),
                    "survival_early": str(Fraction(B, A)),
                    "survival_late": str(Fraction(C-B, D-A)),
                    "survival_difference": str(Fraction(B, A)-Fraction(C-B, D-A)),
                })
        # Exact uniform-label bridge: counts at p<=q are multinomial,
        # independent of the permutation. Fixed rational times, no search.
        def label_h(p, q):
            total = Fraction(0)
            n = b["N"]
            for k in range(n+1):
                for ell in range(k, n+1):
                    mass = (math.comb(n, k)*math.comb(n-k, ell-k) *
                            p**k * (q-p)**(ell-k) * (1-q)**(n-ell))
                    total += mass * Fraction(h(k, ell), b["samples"])
            return total
        a, middle, c = Fraction(1, 3), Fraction(1, 2), Fraction(2, 3)
        A, B, C, D = label_h(a, middle), label_h(a, c), label_h(middle, c), label_h(middle, middle)
        label_covariance = (B*D-A*C)/D**2
        label_difference = B/A-(C-B)/(D-A)
        results.append({
            "lattice": b["lattice"], "L": 3, "permutations": b["samples"],
            "source_file": str(data_path / ref["file"]), "source_sha256": ref["sha256"],
            "triples_checked": count, "violations": violations,
            "rank_insertion_process_is_markov": not violations,
            "exact_uniform_label_bridge": {
                "times": [str(a), str(middle), str(c)],
                "H_a_b": str(A), "H_a_c": str(B), "H_b_c": str(C), "H_b_b": str(D),
                "conditional_covariance": str(label_covariance),
                "survival_difference": str(label_difference),
                "survival_difference_decimal": float(label_difference),
                "markov_identity_violated_at_this_triple": label_covariance != 0,
                "method": "exact multinomial clock mixture of stored complete-permutation H_J; Fraction arithmetic, no random labels drawn",
            },
            "boundary": "exact finite archived histogram under its engine convention; not a label-clock, large-L or continuum claim",
        })
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path,
                    default=HERE.parent / "birth-gap-20260929/cloud-data-5k")
    ap.add_argument("--output-json", type=Path, default=HERE / "result.json")
    ap.add_argument("--output-md", type=Path, default=HERE / "RESULT.md")
    ap.add_argument("--exact-data", type=Path, default=HERE.parent / "birth-gap-20260929/data")
    args = ap.parse_args()
    started = time.perf_counter()
    control = exact_control()
    exact_archives = exact_archive_analysis(args.exact_data)
    manifest_path = args.data / "run.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest["status"] != "completed":
        raise ValueError("incomplete source block")
    groups = defaultdict(list)
    seeds, dependencies = set(), set()
    for ref in manifest["batches"]:  # controls and timing benchmarks never enter
        raw = (args.data / ref["file"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != ref["sha256"]:
            raise ValueError("source hash mismatch: " + ref["file"])
        batch = json.loads(gzip.decompress(raw))
        assert batch["histogram_columns"] == ["J1", "J2", "count"]
        assert batch["samples"] == ref["samples"]
        assert (batch["lattice"], batch["L"]) == (ref["lattice"], ref["L"])
        assert batch["seed"] == ref["seed"]
        assert ref["seed"] not in seeds and ref["dependency_group"] not in dependencies
        seeds.add(ref["seed"]); dependencies.add(ref["dependency_group"])
        rows = batch["histogram"]
        assert len({(j1, j2) for j1, j2, count in rows}) == len(rows)
        assert all(1 <= j1 <= j2 <= batch["N"] and count > 0 for j1, j2, count in rows)
        assert sum(count for j1, j2, count in rows) == batch["samples"]
        groups[batch["lattice"], batch["L"]].append((ref, batch))
    cells = []
    for (lattice, L), items in sorted(groups.items()):
        items.sort(key=lambda x: x[0]["batch"])
        refs, batches = zip(*items)
        cells.append(summarize_cell(list(batches), list(refs), lattice, L))
    result = {
        "schema": "matching-one.birth-markov-kernel.v1",
        "as_of": "2026-09-29",
        "question": "Does entry history still predict completion after conditioning on present rank one?",
        "contract": {
            "time": "occupied-site insertion count, not iid-label time at finite N",
            "cutoffs": "lower quantiles .25/.50/.75 of the equal-weight J1/J2 mixture, by cell",
            "quantiles": QUANTILES,
            "risk": "J1<=b<J2",
            "early": "J1<=a", "late": "a<J1<=b",
            "future": "J2>c",
            "contrast": "P(future|risk,early)-P(future|risk,late)",
            "null": "any time-inhomogeneous Markov description using current rank alone predicts zero",
            "selection": "one triple selected before this scoring, primary L512; lower sizes are context, no cutoff scan",
            "status": "retrospective exploratory reanalysis of an already-inspected random block",
            "uncertainty": "complete equal-batch delete-one; recompute all three cutoffs in every replicate; SE not confidence bounds",
        },
        "source_manifest": str(manifest_path),
        "source_manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "source_production_role": manifest["production_role"],
        "source_engine_commit": manifest["git_base_commit"],
        "source_samples": sum(c["samples"] for c in cells),
        "no_new_samples": True,
        "environment": {"python": sys.version, "platform": platform.platform()},
        "exact_control": control,
        "existing_exact_archives": exact_archives,
        "cells": cells,
        "limitations": [
            "Not prospective or independent model-elimination certification; this block already supplied gap readouts.",
            "Finite-size conditional dependence does not prove non-Markovianity of a continuum limit.",
            "Conditioning on rank one tests the Markov property; it is not a causal effect of earlier birth.",
            "A zero contrast at one triple would not prove Markovianity; the theorem requires all triples.",
            "Square and triangular primitive-coordinate tori have different moduli; no universality comparison.",
            "No fitted exponent, pooled cross-lattice significance, adaptive cutoff search or hidden-state identification.",
        ],
        "wall_seconds": time.perf_counter()-started,
    }
    args.output_json.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Does rank one retain birth-history memory?", "",
             "2026-09-29; retrospective reuse of the completed 560k cloud block. No new samples.", "",
             "One contrast, selected before this scoring: cutoffs a,b,c are the .25/.50/.75 lower quantiles of the pooled J1/J2 births.",
             "Among paths with J1<=b<J2, compare P(J2>c) for early (J1<=a) versus late (a<J1<=b) first birth.",
             "Every current-rank-only time-inhomogeneous Markov process predicts difference zero.", "",
             "## Result", "",
             "L512 is the primary available-size readout; smaller sizes supply context. Errors are 14-batch delete-one SE, not confidence intervals.",
             "All thresholds are re-estimated in every deletion; the JSON retains their values and the complete metric covariance.", "",
             "| Lattice | L | Rank-one risk count | Early survival | Late survival | Difference (percentage points) |",
             "|---|---:|---:|---:|---:|---:|"]
    for cell in cells:
        p, e = cell["point"], cell["jackknife_se"]
        risk = sum(map(sum, cell["risk_table_early_late_by_survive_exit"]))
        lines.append(f"| {cell['lattice']} | {cell['L']} | {risk}/{cell['samples']} | "
                     f"{p['survival_early']:.4f} | {p['survival_late']:.4f} | "
                     f"{100*p['survival_difference']:+.3f} ± {100*e['survival_difference']:.3f} |")
    lines += ["", "## Existing exact L3 archive, a different evidential object", "",
              "No 9! enumeration rerun: the stored complete tables are evaluated with integer/Fraction arithmetic.",
              "All 120 count-time triples 0<=a<b<c<=9 are checked here. This finite exhaustive identity check is separate from the ONE fixed production triple; it is not a production cutoff search.", ""]
    for r in exact_archives:
        lines.append(f"- {r['lattice']} L3: {len(r['violations'])} violating triples out of {r['triples_checked']}; "
                     f"rank-only insertion process Markov under the full-kernel criterion: {r['rank_insertion_process_is_markov']}.")
        for v in r["violations"]:
            lines.append(f"  - {v['cutoffs']}: early/late survival {v['survival_early']} / {v['survival_late']}; "
                         f"difference {v['survival_difference']}; conditional covariance {v['conditional_covariance']}.")
        bridge = r["exact_uniform_label_bridge"]
        lines.append(f"  - Exact iid-label bridge at p=(1/3,1/2,2/3): contrast {bridge['survival_difference']} "
                     f"({100*bridge['survival_difference_decimal']:+.4f} percentage points); "
                     f"Markov identity violated: {bridge['markov_identity_violated_at_this_triple']}.")
    lines += ["", "The exact clock bridge is a multinomial mixture, not a replacement of count k by its expectation Np. Insertion and label clocks have separately evaluated properties; neither finite result is extrapolated to larger sizes or a continuum limit.",
              "", "## Scope", ""] + ["- "+x for x in result["limitations"]]
    lines += ["", "Source: `../birth-gap-20260929/cloud-data-5k/run.json`; all 112 production batches only.",
              "No local-pilot pooling or benchmark rows. Grain: one uniform-permutation filtration, represented by weighted (J1,J2) histogram rows.",
              f"One run took {result['wall_seconds']:.2f} s on Python {platform.python_version()} ({platform.machine()}).",
              "One exact four-cell control is recorded in result.json; old physical enumeration was not rerun.", "",
              "Reproduce from the repository root:", "", "```sh",
              "python3 analysis/birth-markov-kernel-20260929/analyze.py", "```", ""]
    args.output_md.write_text("\n".join(lines))
    print(json.dumps({"seconds": result["wall_seconds"], "samples": result["source_samples"],
                      "primary": [{"lattice": c["lattice"], "point": c["point"], "se": c["jackknife_se"]}
                                  for c in cells if c["L"] == 512]}, indent=2))


if __name__ == "__main__":
    main()
