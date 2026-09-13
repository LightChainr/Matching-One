#!/usr/bin/env python3
"""#633: restore spin-0 Q/Z from the committed N725 histogram.

Reads the #614/#612 histogram in-tree. Uses ``threshold_quantile_lineage``
(correct ``cos 4theta``, invert-then-weight, pooled + aligned delete-one).
Does not rerun Monte Carlo. Does not overwrite ``n725-zflow.json``.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

import threshold_quantile_lineage as L  # noqa: E402

HIST = ROOT / "results/server-20260907/P612-n725-fullcurve/raw/n725_100m.hist.csv"
OLD_ZFLOW = ROOT / "results/probe-invariant-shape/n725-zflow.json"
P612 = ROOT / "results/p612-n725-score/latest.json"
GLIFT = ROOT / "results/probe-invariant-shape/blindness-and-glift.json"
CENSUS = ROOT / "results/probe-invariant-shape/census-exact.json"
OUT = ROOT / "results/probe-invariant-shape/n725-zflow-corrected.json"

DECILES = L.FROZEN_LEVELS
ANCHORS_Z = (0.25, 0.75)
ANCHORS_G = (0.2, 0.8)
LEVELS = tuple(sorted(set(DECILES + ANCHORS_Z)))
A_Z, B_Z = ANCHORS_Z
A_G, B_G = ANCHORS_G


def z_from(Q: list[float], levels: tuple[float, ...], a: float, b: float) -> list[float]:
    idx = {u: i for i, u in enumerate(levels)}
    qa, qb = Q[idx[a]], Q[idx[b]]
    width = qb - qa
    return [(q - qa) / width for q in Q]


def dz_g(Q: list[float], g: list[float], q_levels: tuple[float, ...],
         g_levels: tuple[float, ...], a: float, b: float) -> list[float]:
    """Lift g (on g_levels) into the Z-chart of Q. g and the comparison
    live only on shared levels; default anchors 0.2 and 0.8."""
    q_idx = {u: i for i, u in enumerate(q_levels)}
    g_idx = {u: i for i, u in enumerate(g_levels)}
    W = Q[q_idx[b]] - Q[q_idx[a]]
    ga, gb = g[g_idx[a]], g[g_idx[b]]
    out = []
    for u in g_levels:
        out.append((g[g_idx[u]] - ga) / W - (Q[q_idx[u]] - Q[q_idx[a]]) * (gb - ga) / (W * W))
    return out


def unoriented_angle_deg(u: list[float], v: list[float]) -> dict:
    nu = math.sqrt(sum(x * x for x in u))
    nv = math.sqrt(sum(x * x for x in v))
    if nu == 0.0 or nv == 0.0:
        return {"undetermined": True, "reason": "zero vector"}
    cos = sum(x * y for x, y in zip(u, v)) / (nu * nv)
    cos = max(-1.0, min(1.0, cos))
    raw = math.degrees(math.acos(abs(cos)))
    return {
        "undetermined": False,
        "abs_cosine": abs(cos),
        "unoriented_angle_deg": raw,
        "oriented_angle_deg": math.degrees(math.acos(cos)),
        "supplement_deg": 180.0 - math.degrees(math.acos(cos)),
    }


def main() -> None:
    raw = L.load_batch_histograms(HIST)
    assert raw["n"] == 725
    assert raw["orientation_representative"]["first"] == [26, 7]
    assert raw["orientation_representative"]["second"] == [23, 14]
    cos4 = raw["orientation_cos4theta"]
    w0 = L.spin_zero_weights(cos4)
    leak = sum(w0[k] * cos4[k] for k in cos4)
    assert abs(sum(w0.values()) - 1.0) < 1e-12
    assert abs(leak) < 1e-12

    buggy_c1 = (26 * 26 - 7 * 7) / 725
    buggy_c2 = (23 * 23 - 14 * 14) / 725
    buggy_w = {
        "first": -buggy_c2 / (buggy_c1 - buggy_c2),
        "second": buggy_c1 / (buggy_c1 - buggy_c2),
    }
    old_leak = buggy_w["first"] * cos4["first"] + buggy_w["second"] * cos4["second"]

    p612 = json.loads(P612.read_text())
    p612_w = p612["weightings"]["spin0"]["label_crossing"]["spin0_weights"]
    p612_c = p612["weightings"]["spin0"]["label_crossing"]["orientation_cos4theta"]
    g = p612["weightings"]["spin0"]["consensus_direction_frozen_from_five_committed_transitions"]

    out = {
        "schema": "matching-one.probe-invariant-shape.n725-zflow-corrected.v1",
        "issue": 633,
        "histogram": str(HIST.relative_to(ROOT)),
        "histogram_git_commit": "8b5f9d1a",
        "n": 725,
        "levels_11": list(LEVELS),
        "deciles": list(DECILES),
        "z_anchors": [A_Z, B_Z],
        "g_anchors": [A_G, B_G],
        "orientation_representative": raw["orientation_representative"],
        "orientation_cos4theta": cos4,
        "spin0_weights": dict(w0),
        "spin0_is_interpolation": L.is_interpolation(w0),
        "weight_checks": {
            "sum_w": sum(w0.values()),
            "sum_w_cos4": leak,
            "agrees_with_p612_weights": (
                abs(w0["first"] - p612_w["first"]) < 1e-15
                and abs(w0["second"] - p612_w["second"]) < 1e-15
            ),
            "agrees_with_p612_cos4": (
                abs(cos4["first"] - p612_c["first"]) < 1e-15
                and abs(cos4["second"] - p612_c["second"]) < 1e-15
            ),
        },
        "withdrawn_n725_zflow": {
            "file": str(OLD_ZFLOW.relative_to(ROOT)),
            "formula": "Re(w^2)/|w|^2 = cos(2 theta), misnamed cos_four_theta",
            "reported_cos4": {"first": buggy_c1, "second": buggy_c2},
            "reported_weights": buggy_w,
            "true_spin4_leakage_of_those_weights": old_leak,
            "cdf_then_invert": True,
            "replaced_0.4_and_0.6_with_quarter_anchors": True,
            "mean_of_batchwise_Z_was_reported_as_shape": True,
        },
        "retractions": {
            "full_rank_batch_noise_does_not_kill_one_parameter_mean": (
                "Counterexample: y_N = m + f(N) g + epsilon with epsilon uniform "
                "on {±e_i}. Mean lies on a line; noise covariance is I/d. "
                "Withdrawn from #628 §W3 / #622."
            ),
            "159.6_deg_is_not_almost_orthogonal": (
                "Affine-removed angle 159.6248° between g and ΔZ_{L=4-L=3} is "
                "20.375° between unoriented lines. Use abs cosine. Tiny-torus "
                "ΔZ cannot prove an asymptotic mismatch."
            ),
        },
        "weightings": {},
    }

    for name, weights in (("spin0", w0), ("equal", L.EQUAL_WEIGHTS)):
        jack = L.jackknife_quantiles(raw, weights, LEVELS)
        qcov = L.jackknife_covariance(jack["full"], jack["deleted"])
        Q = jack["full"]
        Z = z_from(Q, LEVELS, A_Z, B_Z)
        idx = {u: i for i, u in enumerate(LEVELS)}
        Zdel = {}
        for batch, qv in jack["deleted"].items():
            qa, qb = qv[idx[A_Z]], qv[idx[B_Z]]
            Zdel[batch] = [(q - qa) / (qb - qa) for q in qv]
        zcov = L.jackknife_covariance(Z, Zdel)
        Q9 = [Q[idx[u]] for u in DECILES]
        g_lift = dz_g(Q, g, LEVELS, DECILES, A_G, B_G)
        z_on_deciles = [Z[idx[u]] for u in DECILES]
        # Compare g-lift to (Z - Z(a_g)) which is already a chart of Q.
        z_gchart = z_from(Q, LEVELS, A_G, B_G)
        z_g_on_deciles = [z_gchart[idx[u]] for u in DECILES]
        out["weightings"][name] = {
            "weights": dict(weights),
            "Q_pooled": Q,
            "Z_pooled_quarter_anchors": Z,
            "Q_se": [math.sqrt(max(0.0, qcov[i][i])) for i in range(len(LEVELS))],
            "Z_se": [math.sqrt(max(0.0, zcov[i][i])) for i in range(len(LEVELS))],
            "Q_deciles": Q9,
            "g_lift_DZ_Q_anchors_0.2_0.8": g_lift,
            "angle_g_lift_vs_Z_in_0.2_0.8_chart": unoriented_angle_deg(g_lift, z_g_on_deciles),
            "batches": jack["batches"],
            "note": "Z_pooled is from the pooled histogram, not the mean of batchwise Z",
        }
        print(name, "Q_pooled", [round(x, 8) for x in Q])
        print(name, "Z_pooled", [round(x, 6) for x in Z])
        print(name, "g-lift angle", out["weightings"][name]["angle_g_lift_vs_Z_in_0.2_0.8_chart"])

    if GLIFT.exists():
        gl = json.loads(GLIFT.read_text())
        raw_ang = gl["direction6"]["angle_raw_deg"]
        aff_ang = gl["direction6"]["angle_after_affine_removal_deg"]
        out["tiny_torus_g_vs_deltaZ"] = {
            "source": str(GLIFT.relative_to(ROOT)),
            "raw_deg": raw_ang,
            "affine_removed_deg": aff_ang,
            "unoriented_deg": min(aff_ang % 180.0, 180.0 - (aff_ang % 180.0)),
            "claim": "report unoriented 20.4°, not 159.6° as almost-orthogonal",
        }

    if CENSUS.exists():
        census = json.loads(CENSUS.read_text())
        out["exact_lab_census_present"] = True
        out["exact_lab_dual_fail_site"] = census.get("site", {}).get("dual_fail")
        out["exact_lab_note"] = (
            "g vs exact-lab ΔZ uses this branch's repaired census; "
            "do not copy PR #628's broken bond Z."
        )
    else:
        out["exact_lab_census_present"] = False

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
