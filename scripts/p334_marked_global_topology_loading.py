#!/usr/bin/env python3
"""Mark full A_top by original direct / collective completion on the safe gate."""
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.stats import binom

from p334_direct_collective_population_loading import birth_channels

ROOT = Path(__file__).resolve().parents[1]
CLOCKS = ROOT / "results/p334-paired-clock-loading"
OUT = ROOT / "results/p334-marked-global-topology-loading"
BIRTH_REV = "9c495ab13e65f2bc93dc0849ee3b73f88724c4b1"
CLOCK_REV = "0d1e586dafbade5e7d1f9bfc598170d0c881e337"
ORIENTATIONS = ("first", "second")
READOUTS = ("canonical", "integrated")
COMPONENTS = ("original_H2_direct_A", "collective_A", "remainder_A")
LABELS = [f"{o}_{r}_{s}" for o in ORIENTATIONS for r in READOUTS for s in COMPONENTS]
PART_LABELS = [f"{o}_{r}_{s}" for o in ORIENTATIONS for r in READOUTS for s in COMPONENTS[:2]]


def main():
    old = json.loads((CLOCKS / "score.json").read_text())
    p_ref = old["contract"]["p_ref"]
    output = {"birth_commit": BIRTH_REV, "clock_commit": CLOCK_REV, "p_ref": p_ref,
              "labels": LABELS, "positive_F2_and_first_birth_debt_labels": PART_LABELS,
              "source_clock_labels": [f"{o}_{q}" for o in ORIENTATIONS for q in ("pi_direct", "pi_collective", "tau_direct", "tau_collective")],
              "gate": "Both orientations rank>=1; original exact_pair clocks required for every R1. Both R2 are identity acceptance. Otherwise keep entire raw global vector.",
              "source_rule": "Only globally accepted R1 directions enter direct or collective; all other complete A_top enters signed remainder.",
              "formula": "A_s=positive_F2_s-pi_s*(1-g(K1)); integral A_direct=(H2*mean_wait-K1*pi_direct)/(N+1)",
              "new_MC": 0, "new_DP": 0, "new_networks": 0, "sizes": {}}
    report = ["# Marked full topology loading on the global safe paired gate", "",
              "All 20000 original counters per size remain. Remainder includes known-past R2 and entire global fallback, and may be signed.", ""]
    for n in (325, 425):
        path = f"results/p334-full-birth-archive/N{n}.csv"
        blob = subprocess.check_output(["git", "show", BIRTH_REV + ":" + path], cwd=ROOT)
        births = [{key: int(value) for key, value in row.items()}
                  for row in csv.DictReader(io.StringIO(blob.decode()))]
        index = {r["counter"]: i for i, r in enumerate(births)}
        k1 = np.array([[r[o + "_k1"] for o in ORIENTATIONS] for r in births])
        k2 = np.array([[r[o + "_k2"] for o in ORIENTATIONS] for r in births])
        ranks = np.array([[r[o + "_rank"] for o in ORIENTATIONS] for r in births])
        f1 = np.stack([binom.sf(k1 - 1, n, p_ref), (n - k1 + 1) / (n + 1)], axis=-1)
        raw_a = np.stack([binom.sf(k1 - 1, n, p_ref) + binom.sf(k2 - 1, n, p_ref) - 1,
                          1 - (k1 + k2) / (n + 1)], axis=-1)
        k0 = births[0]["k0"]
        d = n - k0
        steps = np.arange(1, d + 1)
        kernels = np.array([binom.sf(k0 + steps - 1, n, p_ref),
                            (n - k0 - steps + 1) / (n + 1)])
        batch_values, batch_raw, batch_hybrid = [], [], []
        batch_positive, batch_debt, batch_source_clock = [], [], []
        batch_gate_counts = []
        gate_counts, rank_counts, marked_counts = Counter(), Counter(), Counter()
        source_hashes = {}
        max_additivity_error = max_direct_integral_identity_error = 0.
        for batch in range(20):
            cp = CLOCKS / "batches" / f"N{n}.batch{batch:02d}.json.gz"
            cblob = cp.read_bytes()
            source_hashes[str(cp.relative_to(ROOT))] = hashlib.sha256(cblob).hexdigest()
            records = json.loads(gzip.decompress(cblob))["records"]
            ix = np.array([index[r["counter"]] for r in records])
            raw = raw_a[ix]
            hybrid = raw.copy()
            values = np.zeros((1000, 2, 2, 3))
            values[:, :, :, 2] = raw
            positive = np.zeros((1000, 2, 2, 2))
            debt = np.zeros_like(positive)
            source_clock = np.zeros((1000, 2, 4))
            local_counts = Counter()
            for j, record in enumerate(records):
                ii = ix[j]
                pair_ranks = ranks[ii]
                rank_counts["".join(str(r) for r in pair_ranks)] += 1
                if min(pair_ranks) == 0:
                    status = "fallback_contains_R0"
                elif min(pair_ranks) == 2:
                    status = "accepted_both_R2_identity"
                elif record["status"] == "exact_pair":
                    status = "accepted_exact_R1_pair"
                else:
                    status = "fallback_original_clock_policy"
                gate_counts[status] += 1
                local_counts[status] += 1
                if status != "accepted_exact_R1_pair":
                    continue
                for o, rank in enumerate(pair_ranks):
                    if rank != 1:
                        continue
                    row, clock = record["source_rows"][o], record["clocks"][o]
                    channels = birth_channels(clock["safe_coefficients"], row["H2"])
                    pi = channels.sum(axis=1)
                    tau = channels @ steps
                    source_clock[j, o] = np.r_[pi, tau]
                    positive[j, o] = kernels @ channels.T
                    debt[j, o] = (1 - f1[ii, o])[:, None] * pi[None, :]
                    values[j, o, :, :2] = positive[j, o] - debt[j, o]
                    values[j, o, :, 2] = 0.
                    hybrid[j, o] = f1[ii, o] + np.array(clock["conditional"]) - 1
                    shortcut = (row["H2"] * clock["mean_wait"] - k1[ii, o] * pi[0]) / (n + 1)
                    max_direct_integral_identity_error = max(max_direct_integral_identity_error,
                        abs(float(values[j, o, 1, 0]) - shortcut))
                    marked_counts[ORIENTATIONS[o]] += 1
            max_additivity_error = max(max_additivity_error,
                float(np.max(np.abs(values.sum(axis=-1) - hybrid))))
            batch_values.append(values.mean(axis=0).reshape(-1))
            batch_positive.append(positive.mean(axis=0).reshape(-1))
            batch_debt.append(debt.mean(axis=0).reshape(-1))
            batch_source_clock.append(source_clock.mean(axis=0).reshape(-1))
            batch_raw.append(raw.mean(axis=0).reshape(-1))
            batch_hybrid.append(hybrid.mean(axis=0).reshape(-1))
            batch_gate_counts.append(dict(local_counts))
        vectors = np.array(batch_values)
        means = vectors.mean(axis=0).reshape(2, 2, 3)
        positives = np.array(batch_positive)
        debts = np.array(batch_debt)
        delta = old["sizes"][str(n)]["delta_cos4"]
        h4 = (vectors.reshape(20, 2, 2, 3)[:, 0] - vectors.reshape(20, 2, 2, 3)[:, 1]) / delta
        h4_positive = (positives.reshape(20, 2, 2, 2)[:, 0] - positives.reshape(20, 2, 2, 2)[:, 1]) / delta
        h4_debt = (debts.reshape(20, 2, 2, 2)[:, 0] - debts.reshape(20, 2, 2, 2)[:, 1]) / delta
        output["sizes"][str(n)] = {
            "k0": k0, "d": d, "counter_count": 20000, "batch_ids": list(range(20)),
            "counter_count_per_batch": 1000, "delta_cos4": delta,
            "birth_path": path, "birth_sha256": hashlib.sha256(blob).hexdigest(),
            "clock_batch_sha256": source_hashes,
            "gate_counts": dict(gate_counts), "joint_rank_counts": dict(rank_counts),
            "marked_R1_orientation_counts": dict(marked_counts), "batch_gate_counts": batch_gate_counts,
            "means": means.reshape(-1).tolist(),
            "joint_20_batch_means_orientation_readout_source": vectors.tolist(),
            "positive_F2_20_batch_means": positives.tolist(),
            "first_birth_debt_20_batch_means": debts.tolist(),
            "source_clock_20_batch_means": np.array(batch_source_clock).tolist(),
            "global_A_fields": [f"{o}_{r}" for o in ORIENTATIONS for r in READOUTS],
            "raw_A_20_batch_means": np.array(batch_raw).tolist(),
            "safe_global_hybrid_A_20_batch_means": np.array(batch_hybrid).tolist(),
            "H4_labels": [f"{r}_{s}" for r in READOUTS for s in COMPONENTS],
            "H4_20_batch_means": h4.reshape(20, 6).tolist(),
            "H4_means": h4.mean(axis=0).tolist(),
            "H4_positive_F2_means": h4_positive.mean(axis=0).tolist(),
            "H4_first_birth_debt_means": h4_debt.mean(axis=0).tolist(),
            "max_record_component_sum_minus_safe_hybrid_abs": max_additivity_error,
            "max_direct_integral_shortcut_abs": max_direct_integral_identity_error}
        report += [f"## N{n}", "", f"Global gate: {dict(gate_counts)}; marked R1 counts: {dict(marked_counts)}.", "",
                   "| H4 contrast | Full marked direct A | Full marked collective A | Remainder A | Full safe-global A |",
                   "|---|---:|---:|---:|---:|"]
        for ri, readout in enumerate(READOUTS):
            v = h4.mean(axis=0)[ri]
            report.append(f"| {readout} | {v[0]:.10g} | {v[1]:.10g} | {v[2]:.10g} | {v.sum():.10g} |")
        report += ["", "Same-gate marked-source parts (full marked A = positive F2 - first-birth debt):", "",
                   "| H4 readout/source | Positive completion F2 | First-birth debt | Full marked A |",
                   "|---|---:|---:|---:|"]
        for ri, readout in enumerate(READOUTS):
            for si, component in enumerate(COMPONENTS[:2]):
                pos, charge = h4_positive.mean(axis=0)[ri, si], h4_debt.mean(axis=0)[ri, si]
                report.append(f"| {readout} {component} | {pos:.10g} | {charge:.10g} | {pos-charge:.10g} |")
        report.append("")
        print(f"N{n} gate={dict(gate_counts)} H4={h4.mean(axis=0).tolist()} max_additivity={max_additivity_error}", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "score.json").write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    (OUT / "REPORT.md").write_text("\n".join(report))


if __name__ == "__main__":
    main()
