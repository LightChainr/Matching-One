#!/usr/bin/env python3
"""Original-H2 / collective birth readout of stored whole-pair exact clocks.

No geometry construction, reliability DP, sampling, or suffix reclassification.
All failed whole-pair replacements remain unclassified original observations.
"""
from collections import Counter
import gzip
import hashlib
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/p334-paired-clock-loading"
OUT = ROOT / "results/p334-direct-collective-population-loading"
SOURCE_COMMIT = "0d1e586dafbade5e7d1f9bfc598170d0c881e337"
ORIENTATIONS = ("first", "second")
READOUTS = ("canonical", "integrated")
COMPONENTS = ("original_H2_direct", "collective", "unclassified_original_Y")
LABELS = [f"{orientation}_{readout}_{component}"
          for orientation in ORIENTATIONS for readout in READOUTS
          for component in COMPONENTS]


def birth_channels(coefficients, h):
    """Integer numerators before float division; j is 1..d, never truncated."""
    d = len(coefficients) - 1
    direct, collective = [], []
    for j in range(1, d + 1):
        previous, current = coefficients[j - 1:j + 1]
        denominator = comb(d, j - 1) * (d - j + 1)
        direct_numerator = h * previous
        collective_numerator = (d - j + 1 - h) * previous - j * current
        if direct_numerator < 0 or collective_numerator < 0:
            raise ValueError("Negative exact source-channel count")
        direct.append(direct_numerator / denominator)
        collective.append(collective_numerator / denominator)
    return np.array([direct, collective])


def projection_summary(means, delta):
    """Descriptive means and allocation bounds, not confidence intervals."""
    result = {}
    for ri, readout in enumerate(READOUTS):
        first, second = means[:, ri, :]
        direct, collective, unknown = (first - second) / delta
        intervals = {}
        for ci, component in enumerate(COMPONENTS[:2]):
            difference = first[ci] - second[ci]
            endpoints = [(difference - second[2]) / delta,
                         (difference + first[2]) / delta]
            intervals[component] = sorted(endpoints)
        abs_channels = abs(direct) + abs(collective)
        result[readout] = {
            "H4_components": dict(zip(COMPONENTS, [direct, collective, unknown])),
            "H4_total": float(direct + collective + unknown),
            "classified_direct_collective_opposite_sign": bool(direct * collective < 0),
            "classified_absolute_loading_cancelled_fraction":
                float(1 - abs(direct + collective) / abs_channels) if abs_channels else None,
            "unknown_allocation_envelopes_H4": intervals,
            "envelope_interpretation": "Only allocation of this hybrid estimator's unclassified nonnegative Y mass to direct or collective; neither population uncertainty nor a confidence interval. Source allocations are coupled, not two independent free intervals."}
    return result


def main():
    source_score = json.loads((SOURCE / "score.json").read_text())
    contract = source_score["contract"]
    p_ref = contract["p_ref"]
    result = {
        "source_commit": SOURCE_COMMIT,
        "source_path": "results/p334-paired-clock-loading",
        "source_score_sha256": hashlib.sha256((SOURCE / "score.json").read_bytes()).hexdigest(),
        "p_ref": p_ref, "orientations": ORIENTATIONS, "readouts": READOUTS,
        "components": COMPONENTS, "labels": LABELS,
        "formula": "P(T=j,V_final in original D)=H2*f[j-1]/(C(d,j-1)*(d-j+1)); collective=((d-j+1-H2)*f[j-1]-j*f[j])/(C(d,j-1)*(d-j+1))",
        "source_semantics": "Original D is the checkpoint's H2 singleton-trigger set. Collective includes all other final sites, including new singleton triggers made by previous safe insertions.",
        "fallback_policy": "Every whole_pair_fallback keeps both original Y orientations entirely in unclassified, ignoring any partially saved exact clock. Outside-rank-one defined contributions are zero.",
        "new_MC": 0, "new_DP": 0, "new_network_constructions": 0,
        "sizes": {}}
    report = ["# Direct and collective sources of the paired R1 population loading", "",
              "Descriptive source decomposition of the fixed hybrid estimator; all 20000 counters remain in each size's denominator.", ""]
    for n in (325, 425):
        old = source_score["sizes"][str(n)]
        k0 = 193 if n == 325 else 252
        d = n - k0
        thresholds = k0 + np.arange(1, d + 1)
        kernels = np.array([binom.sf(thresholds - 1, n, p_ref),
                            (n - thresholds + 1) / (n + 1)])
        batch_vectors, batch_hybrid = [], []
        status_counts, exact_orientation_counts = Counter(), Counter()
        source_hashes = {}
        max_record_additivity_error = 0.
        max_total_birth_mass_error = 0.
        for batch in range(20):
            path = SOURCE / "batches" / f"N{n}.batch{batch:02d}.json.gz"
            blob = path.read_bytes()
            source_hashes[str(path.relative_to(ROOT))] = hashlib.sha256(blob).hexdigest()
            records = json.loads(gzip.decompress(blob))["records"]
            if len(records) != 1000:
                raise ValueError("Original 1000-counter batch is incomplete")
            values = np.zeros((1000, 2, 2, 3))
            hybrid = np.array([r["Y"] for r in records])
            for index, record in enumerate(records):
                status = record["status"]
                status_counts[status] += 1
                if status == "whole_pair_fallback":
                    for o in range(2):
                        values[index, o, :, 2] = hybrid[index, [o, 2 + o]]
                elif status == "exact_pair":
                    for o, row in enumerate(record["source_rows"]):
                        if row is None:
                            continue
                        exact_orientation_counts[ORIENTATIONS[o]] += 1
                        if row["k0"] != k0:
                            raise ValueError("Unexpected checkpoint outside the declared source")
                        channels = birth_channels(record["clocks"][o]["safe_coefficients"], row["H2"])
                        max_total_birth_mass_error = max(max_total_birth_mass_error,
                                                         abs(float(channels.sum()) - 1))
                        values[index, o, :, :2] = kernels @ channels.T
                elif status != "outside_rank_one":
                    raise ValueError("Unrecognized whole-pair policy status")
            totals = values.sum(axis=-1)
            reconstructed = np.column_stack([totals[:, 0, 0], totals[:, 1, 0],
                                              totals[:, 0, 1], totals[:, 1, 1]])
            max_record_additivity_error = max(max_record_additivity_error,
                                               float(np.max(np.abs(reconstructed - hybrid))))
            batch_vectors.append(values.mean(axis=0).reshape(-1))
            batch_hybrid.append(hybrid.mean(axis=0))
        vectors = np.array(batch_vectors)
        means = vectors.mean(axis=0).reshape(2, 2, 3)
        delta = old["delta_cos4"]
        projections = projection_summary(means, delta)
        hybrid_batches = np.array(batch_hybrid)
        # Root's existing columns: X-six then Y-six; select four unprojected Y.
        original_saved = np.array(old["joint_20_batch_means_X_then_Y"])[:, [6, 7, 9, 10]]
        h4_vectors = (vectors.reshape(20, 2, 2, 3)[:, 0] -
                      vectors.reshape(20, 2, 2, 3)[:, 1]) / delta
        result["sizes"][str(n)] = {
            "k0": k0, "d": d, "original_counter_count": 20000,
            "batch_ids": list(range(20)), "counter_count_per_batch": 1000,
            "pair_status_counts": dict(status_counts),
            "exact_R1_orientation_counts": dict(exact_orientation_counts),
            "source_batch_sha256": source_hashes,
            "delta_cos4": delta, "means": means.reshape(-1).tolist(),
            "joint_20_batch_means_orientation_readout_source": vectors.tolist(),
            "hybrid_fields": ["first_canonical", "second_canonical", "first_integrated", "second_integrated"],
            "hybrid_20_batch_means": hybrid_batches.tolist(),
            "H4_labels": [f"{readout}_{component}" for readout in READOUTS for component in COMPONENTS],
            "H4_20_batch_means": h4_vectors.reshape(20, 6).tolist(),
            "projection_summary": projections,
            "max_record_component_sum_minus_original_hybrid_abs": max_record_additivity_error,
            "max_exact_clock_total_birth_mass_minus_one_abs": max_total_birth_mass_error,
            "max_saved_batch_hybrid_reproduction_abs": float(np.max(np.abs(hybrid_batches - original_saved)))}
        report += [f"## N{n}", "", f"Whole-pair statuses: {dict(status_counts)}.", "",
                   "| Orientation/readout | Original H2 direct | Collective | Unclassified original Y | Total |",
                   "|---|---:|---:|---:|---:|"]
        for o, orientation in enumerate(ORIENTATIONS):
            for ri, readout in enumerate(READOUTS):
                v = means[o, ri]
                report.append(f"| {orientation} {readout} | {v[0]:.10g} | {v[1]:.10g} | {v[2]:.10g} | {v.sum():.10g} |")
        report += ["", "| H4-normalized contrast | Direct | Collective | Unclassified | Total | Classified cancellation |",
                   "|---|---:|---:|---:|---:|---:|"]
        for readout, projection in projections.items():
            v = projection["H4_components"]
            report.append(f"| {readout} | {v[COMPONENTS[0]]:.10g} | {v[COMPONENTS[1]]:.10g} | {v[COMPONENTS[2]]:.10g} | {projection['H4_total']:.10g} | {projection['classified_absolute_loading_cancelled_fraction']:.6g} |")
        report += ["", "Unknown-source allocation envelopes (not confidence intervals):",
                   json.dumps({r: p["unknown_allocation_envelopes_H4"] for r, p in projections.items()}), ""]
        print(f"N{n}", json.dumps(projections), flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "score.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    (OUT / "REPORT.md").write_text("\n".join(report) + "\n")


if __name__ == "__main__":
    main()
