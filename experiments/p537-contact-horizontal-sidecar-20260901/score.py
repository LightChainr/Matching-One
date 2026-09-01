#!/usr/bin/env python3
"""Gauge-horizontal sidecar for the frozen P537 N65 contact carrier.

This is post-processing only.  It consumes the already frozen contact-stage
tables and the independent N65 baseline; it does not alter the canonical
selection or draw random configurations.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.stats import binom


N = 65
BATCHES = 100
DELTA = 1152 / 845
FIELDS = (
    "count", "sum_q0", "sum_E0", "sum_a16_0", "sum_q0_a16_0",
    "sum_E0_a16_0", "sum_q1", "sum_E1", "sum_a16_1",
    "sum_q1_a16_1", "sum_E1_a16_1",
)
OUTPUTS = (
    "p", "M_t", "T_full", "T_contact_canonical", "C_contact",
    "beta_C_contact", "T_contact_horizontal", "T_remainder_canonical",
    "T_remainder_horizontal", "contact_identity_residual",
)


def load_baseline(path: Path) -> tuple[np.ndarray, np.ndarray]:
    raw = np.zeros((BATCHES, 2, 2, N + 1), dtype=np.float64)
    samples = np.zeros(BATCHES, dtype=np.float64)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if int(row["n"]) != N:
                continue
            b = int(row["batch"])
            g = ("first", "second").index(row["orientation"])
            t = ("minus", "plus").index(row["kind"])
            raw[b, g, t, int(row["k"])] += int(row["count"])
            samples[b] = int(row["samples"])
    if np.any(samples <= 0):
        raise ValueError("incomplete N65 baseline")
    return raw, samples


def load_tables(paths: list[Path]) -> tuple[dict, float, dict]:
    rows: dict[tuple, np.ndarray] = defaultdict(
        lambda: np.zeros(len(FIELDS), dtype=np.int64)
    )
    metadata = []
    for path in paths:
        meta: dict[str, str] = {}
        body = []
        for line in path.read_text().splitlines():
            if line.startswith("# "):
                key, value = line[2:].split("=", 1)
                meta[key] = value
            elif not line.startswith("#"):
                body.append(line)
        metadata.append(meta)
        for row in csv.DictReader(body, delimiter="\t"):
            key = (
                int(row["batch"]), row["kind"],
                ("axis", "tilted").index(row["geometry"]),
                int(row["dx"]), int(row["dy"]), row["stage"],
                int(row["contact_mask"]), int(row["k"]),
            )
            rows[key] += np.array([int(row[field]) for field in FIELDS], dtype=np.int64)
    shards = int(metadata[0]["shard_count"])
    samples = int(metadata[0]["samples"])
    seed = metadata[0]["seed"]
    proposal = float(metadata[0]["proposal_p"])
    if {int(item["shard_index"]) for item in metadata} != set(range(shards)):
        raise ValueError("incomplete shards")
    if any(
        (int(item["samples"]), item["seed"], float(item["proposal_p"]))
        != (samples, seed, proposal)
        for item in metadata
    ):
        raise ValueError("mixed shard contracts")
    return rows, proposal, {
        "samples": samples, "shards": shards, "seed": seed,
        "proposal_p": proposal,
    }


def baseline_at(raw: np.ndarray, samples: np.ndarray, p: float, omit: int | None):
    hist = raw.sum(axis=0) - (raw[omit] if omit is not None else 0)
    total = samples.sum() - (samples[omit] if omit is not None else 0)
    cumulative = np.cumsum(hist, axis=-1)
    q = (-total + cumulative[:, 0] + cumulative[:, 1]) / total
    e = (total - cumulative[:, 0] + cumulative[:, 1]) / total
    k = np.arange(N + 1)
    weight = binom.pmf(k, N, p)
    score = k - N * p
    score2 = score * score - N * p * (1 - p)
    return {
        "q": q @ weight, "e": e @ weight,
        "qt": q @ (weight * score), "et": e @ (weight * score),
        "qtt": q @ (weight * score2), "ett": e @ (weight * score2),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--tables", required=True, nargs=4, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)

    baseline_raw, baseline_samples = load_baseline(args.baseline)
    rows, proposal, run = load_tables(args.tables)
    records = list(rows.items())
    displacements = sorted({key[3:5] for key, _ in records if key[1] == "global"})
    displacement_index = {value: i for i, value in enumerate(displacements)}
    dcount = len(displacements)

    def arrays(kind: str):
        selected = [(key, value) for key, value in records if key[1] == kind]
        return (
            np.array([key[0] for key, _ in selected], dtype=np.int16),
            np.array([key[2] for key, _ in selected], dtype=np.int8),
            np.array([displacement_index[key[3:5]] for key, _ in selected], dtype=np.int16),
            np.array([key[7] for key, _ in selected], dtype=np.int16),
            np.stack([value for _, value in selected]),
        )

    global_b, global_g, global_d, global_k, global_values = arrays("global")
    carrier_b, carrier_g, carrier_d, carrier_k, carrier_values = arrays("carrier")
    probe = 0
    probe_mask = (global_g == 0) & (global_d == probe)
    batch_n = np.bincount(
        global_b[probe_mask], weights=global_values[probe_mask, 0], minlength=BATCHES
    )
    if np.any(batch_n <= 0):
        raise ValueError("empty producer batch")

    def evaluate(base_omit: int | None = None, source_omit: int | None = None,
                 detail: bool = False):
        def at(p: float):
            return baseline_at(baseline_raw, baseline_samples, p, base_omit)

        p = brentq(lambda value: at(value)["q"].mean(), 0.58, 0.61)
        base = at(p)
        mt = base["qt"].mean()
        mtt = base["qtt"].mean()
        yt = (base["et"][0] - base["et"][1]) / DELTA
        ytt = (base["ett"][0] - base["ett"][1]) / DELTA
        ratio = yt / mt
        ratio_t = (ytt - ratio * mtt) / mt
        c = (1 / DELTA, -1 / DELTA)
        mu_h = [2 * c[g] * base["e"][g] - ratio * base["q"][g] for g in range(2)]
        denominator = batch_n.sum() - (
            batch_n[source_omit] if source_omit is not None else 0
        )

        # Global packets establish mu_a and the displacement-specific beta.
        keep_global = np.ones(len(global_b), dtype=bool)
        if source_omit is not None:
            keep_global &= global_b != source_omit
        gb_g = global_g[keep_global]
        gb_d = global_d[keep_global]
        gb_k = global_k[keep_global]
        gb_v = global_values[keep_global]
        flat = gb_g.astype(np.int64) * dcount + gb_d
        f_global = (
            (1 - p) * (p / proposal) ** gb_k
            * ((1 - p) / (1 - proposal)) ** (N - 2 - gb_k) / denominator
        )
        packets = np.zeros((2 * dcount, 6), dtype=np.float64)
        for state, columns in enumerate(((3, 4, 5), (8, 9, 10))):
            w = (1 - p, p)[state]
            score = gb_k + state - N * p
            values = gb_v[:, list(columns)] / (N * 16)
            for field in range(3):
                packets[:, field] += np.bincount(
                    flat, weights=f_global * w * values[:, field], minlength=2 * dcount
                )
                packets[:, field + 3] += np.bincount(
                    flat, weights=f_global * w * score * values[:, field],
                    minlength=2 * dcount,
                )
        packets = packets.reshape(2, dcount, 6)

        mu_a = packets[..., 0]
        covariance = packets[..., 1] - base["q"][:, None] * mu_a
        beta = covariance.sum(axis=0) / (2 * mt)

        # Complete full response, with the same exact C4 missing-NN fill used by
        # the committed full-T scorer.
        retained = packets.sum(axis=1)
        nn = ((-1, 0), (0, -1), (0, 1))
        nn_indices = [displacement_index[d] for d in nn]
        omitted = packets[:, nn_indices].mean(axis=1)
        full = retained + omitted
        a, qa, ea, at_, qat, eat = full.T
        jm = qa - base["q"] * a
        jmt = qat - base["qt"] * a - base["q"] * at_
        jy = (ea - base["e"] * a)
        jyt = eat - base["et"] * a - base["e"] * at_
        t_full = (
            (jyt[0] - jyt[1]) / DELTA
            - ratio * jmt.mean() - ratio_t * jm.mean()
        )

        per_d = np.zeros((dcount, 3), dtype=np.float64)
        v = p * (1 - p)
        keep_carrier = np.ones(len(carrier_b), dtype=bool)
        if source_omit is not None:
            keep_carrier &= carrier_b != source_omit
        cb_g = carrier_g[keep_carrier]
        cb_d = carrier_d[keep_carrier]
        cb_k = carrier_k[keep_carrier]
        cb_v = carrier_values[keep_carrier]
        f_carrier = (
            (1 - p) * (p / proposal) ** cb_k
            * ((1 - p) / (1 - proposal)) ** (N - 2 - cb_k) / denominator
        )
        mu = mu_a[cb_g, cb_d]
        b_d = beta[cb_d]
        s_minus = cb_k - (N - 1) * p
        for state, columns in enumerate(((1, 2, 3, 4, 5), (6, 7, 8, 9, 10))):
            q, e, a16, qa16, ea16 = cb_v[:, list(columns)].T
            count = cb_v[:, 0]
            aa, qaa, eaa = np.array([a16, qa16, ea16]) / (N * 16)
            w = (1 - p, p)[state]
            u = state - p
            score = s_minus + u
            bterm = u * score - v
            sh = 2 * np.take(c, cb_g) * e - ratio * q - np.take(mu_h, cb_g) * count
            sha = (
                2 * np.take(c, cb_g) * (eaa - mu * e)
                - ratio * (qaa - mu * q)
                - np.take(mu_h, cb_g) * (aa - mu * count)
            )
            original = 2 * f_carrier * w * (u * sha - b_d * bterm * sh)
            cc = 2 * f_carrier * w * v * sh
            per_d[:, 0] += np.bincount(cb_d, weights=original, minlength=dcount)
            per_d[:, 1] += np.bincount(cb_d, weights=cc, minlength=dcount)
            per_d[:, 2] += np.bincount(
                cb_d, weights=original - b_d * cc, minlength=dcount
            )

        contact = per_d.sum(axis=0)
        original, cc, horizontal = contact
        beta_cc = original - horizontal
        vector = np.array([
            p, mt, t_full, original, cc, beta_cc, horizontal,
            t_full - original, t_full - horizontal,
            original - horizontal - beta_cc,
        ])
        extra = None
        if detail:
            extra = {
                f"{dx},{dy}": {
                    "beta": beta[index],
                    "T_contact_canonical": per_d[index, 0],
                    "C_contact": per_d[index, 1],
                    "beta_C_contact": beta[index] * per_d[index, 1],
                    "T_contact_horizontal": per_d[index, 2],
                }
                for index, (dx, dy) in enumerate(displacements)
            }
        return vector, extra

    center, per_displacement = evaluate(detail=True)
    source_jk = np.array([evaluate(source_omit=b)[0] for b in range(BATCHES)])
    baseline_jk = np.array([evaluate(base_omit=b)[0] for b in range(BATCHES)])
    factor = (BATCHES - 1) / BATCHES
    source_delta = source_jk - source_jk.mean(axis=0)
    baseline_delta = baseline_jk - baseline_jk.mean(axis=0)
    covariance = factor * (source_delta.T @ source_delta + baseline_delta.T @ baseline_delta)
    se = np.sqrt(np.maximum(0, np.diag(covariance)))
    estimates = {
        name: {
            "value": float(value), "se": float(error),
            "ci95": [float(value - 1.96 * error), float(value + 1.96 * error)],
        }
        for name, value, error in zip(OUTPUTS, center, se)
    }
    payload = {
        "schema": "matching-one/p537-contact-horizontal-sidecar/v1",
        "status": "COMPLETED_EXISTING_SUFFICIENT_STATISTICS",
        "identity": "T_contact_canonical = T_contact_horizontal + sum_d beta_d C_contact,d",
        "gauge_law": "a -> a + c K + d implies T_C -> T_C + c C_C; T_C_horizontal invariant",
        "N65": estimates,
        "per_displacement": per_displacement,
        "joint_covariance": {
            "order": list(OUTPUTS), "matrix": covariance.tolist(),
            "groups": ["new_N65_source_MC_100_batches", "P45_N65_baseline_100_batches"],
        },
        "delete_one": {
            "order": list(OUTPUTS),
            "source": source_jk.tolist(), "baseline": baseline_jk.tolist(),
        },
        "run": run,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    print(json.dumps({name: estimates[name] for name in OUTPUTS[2:9]}, indent=2))


if __name__ == "__main__":
    main()
