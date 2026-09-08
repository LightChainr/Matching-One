#!/usr/bin/env python3
"""#661: compatible g-lift vs exact-lab ΔZ, on existing objects only.

Leftover of #633/#655.  Not a new production, not #622.  Reads:

* frozen g — ``results/type582-residual/latest.json`` ``consensus_g``
  (cross-checked against ``results/p612-n725-score/latest.json`` spin0
  consensus; 725 is not in the freeze),
* the repaired exact laboratory — ``results/probe-invariant-shape/census-exact.json``
  (PR #653, bond ``dual_fail = 0``; PR #628's broken bond Z is NOT copied),
* the historical withdrawn reading — ``blindness-and-glift.json``.

Formula (issue #661):

    DZ_Q[g](u) = (g(u) - g(a))/W - (Q(u) - Q(a)) (g(b) - g(a))/W^2
    W = Q(b) - Q(a)

Shared anchors a = 0.2, b = 0.8 (both on g's 9-decile grid).  Angles are
UNORIENTED (|cosine|).  The exact-lab comparison object is the site tangent
ΔZ = Z_{L=4} − Z_{L=3} on the same deciles, the object named in #633 and
recomputed here from the repaired census.  The lift is invariant under
affine reparametrisation of Q (Q ↦ cQ + d scales DZ_Q[g] by 1/c), so
recovering Q on the deciles from the census Z (which fixes Q up to affine)
is exact for the direction.

Covariance: the exact laboratory is deterministic (census), so there is no
sampling covariance for ΔZ_exact; the N725 covariance lives in
``n725-zflow-corrected.json`` and is cited, not recomputed (no rerun).
No N725 Monte Carlo is run.  Equal weighting is a labelled sensitivity.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

CENSUS = ROOT / "results/probe-invariant-shape/census-exact.json"
TYPE582 = ROOT / "results/type582-residual/latest.json"
P612 = ROOT / "results/p612-n725-score/latest.json"
GLIFT_HIST = ROOT / "results/probe-invariant-shape/blindness-and-glift.json"
N725_CORR = ROOT / "results/probe-invariant-shape/n725-zflow-corrected.json"
OUT = ROOT / "results/probe-invariant-shape/issue661-glift-exact-lab.json"

LEVELS = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
A_G, B_G = 0.2, 0.8


def unoriented_angle_deg(u: list[float], v: list[float]) -> dict:
    nu = math.sqrt(sum(x * x for x in u))
    nv = math.sqrt(sum(x * x for x in v))
    if nu == 0.0 or nv == 0.0:
        return {"undetermined": True, "reason": "zero vector"}
    cos = sum(x * y for x, y in zip(u, v)) / (nu * nv)
    cos = max(-1.0, min(1.0, cos))
    oriented = math.degrees(math.acos(cos))
    return {
        "undetermined": False,
        "abs_cosine": abs(cos),
        "unoriented_angle_deg": math.degrees(math.acos(abs(cos))),
        "oriented_angle_deg": oriented,
        "supplement_deg": 180.0 - oriented,
    }


def dz_g_lift(Q: list[float], g: list[float], a: float = A_G, b: float = B_G) -> list[float]:
    """DZ_Q[g](u) = (g(u)-g(a))/W - (Q(u)-Q(a))(g(b)-g(a))/W^2."""
    idx = {u: i for i, u in enumerate(LEVELS)}
    W = Q[idx[b]] - Q[idx[a]]
    if W == 0.0:
        raise ValueError("degenerate anchor width W = Q(0.8) - Q(0.2)")
    ga, gb = g[idx[a]], g[idx[b]]
    return [
        (g[idx[u]] - ga) / W - (Q[idx[u]] - Q[idx[a]]) * (gb - ga) / (W * W)
        for u in LEVELS
    ]


def q_on_deciles_from_census(zdec: list[float], q_quarter: str | float,
                             q_threequarters: str | float) -> list[float]:
    """Recover Q on the deciles from the census Z, exactly up to affine.

    Z(u) = (Q(u) - Q(1/4)) / (Q(3/4) - Q(1/4)), so
    Q(u) = Q(1/4) + W_z Z(u).  The lift direction is invariant under
    Q ↦ cQ + d (verified: DZ_{cQ+d}[g] = DZ_Q[g]/c), so this is exact
    for the angle, not an approximation.
    """
    q4 = float(Fraction(q_quarter)) if isinstance(q_quarter, str) else q_quarter
    q34 = float(Fraction(q_threequarters)) if isinstance(q_threequarters, str) else q_threequarters
    wz = q34 - q4
    return [q4 + wz * z for z in zdec]


def main() -> None:
    census = json.loads(CENSUS.read_text())
    type582 = json.loads(TYPE582.read_text())
    p612 = json.loads(P612.read_text())

    g = type582["consensus_g"]
    assert type582["levels"] == list(LEVELS)
    g612 = p612["weightings"]["spin0"][
        "consensus_direction_frozen_from_five_committed_transitions"]
    assert g == g612, "type582 consensus_g and p612 spin0 consensus disagree"
    g_equal = p612["weightings"]["equal"][
        "consensus_direction_frozen_from_five_committed_transitions"]

    # repaired exact laboratory (PR #653): bond dual_fail must be 0
    assert census["bond"]["3"]["dual_fail"] == 0, (
        "census-exact.json is not the repaired #653 version"
    )

    z3 = census["site"]["3"]["Z_levels_float"]
    z4 = census["site"]["4"]["Z_levels_float"]
    zb = census["bond"]["3"]["Z_levels_float"]
    dz_site = [b - a for a, b in zip(z3, z4)]

    q3 = q_on_deciles_from_census(
        z3, census["site"]["3"]["Q_quarter"], census["site"]["3"]["Q_threequarters"])
    q4 = q_on_deciles_from_census(
        z4, census["site"]["4"]["Q_quarter"], census["site"]["4"]["Q_threequarters"])
    qb = q_on_deciles_from_census(
        zb, census["bond"]["3"]["Q_quarter"], census["bond"]["3"]["Q_threequarters"])

    out = {
        "schema": "matching-one.probe-invariant-shape.issue661-glift-exact-lab.v1",
        "issue": 661,
        "leftover_of": 633,
        "not_a_new_production": True,
        "not_622": True,
        "n725_rerun": False,
        "levels": list(LEVELS),
        "g_anchors": [A_G, B_G],
        "g": g,
        "g_provenance": [
            "results/type582-residual/latest.json consensus_g",
            "results/p612-n725-score/latest.json weightings.spin0."
            "consensus_direction_frozen_from_five_committed_transitions",
        ],
        "g_freeze_excludes_725": True,
        "inputs": {
            "census": "results/probe-invariant-shape/census-exact.json (PR #653 repaired; bond dual_fail = 0)",
            "n725_corrected": "results/probe-invariant-shape/n725-zflow-corrected.json (PR #655, cited not recomputed)",
            "historical_withdrawn": "results/probe-invariant-shape/blindness-and-glift.json",
        },
        "formula": "DZ_Q[g](u) = (g(u)-g(a))/W - (Q(u)-Q(a))(g(b)-g(a))/W^2, W = Q(b)-Q(a), a=0.2 b=0.8",
        "chart_note": (
            "Q on the deciles is recovered from the census Z exactly up to affine; "
            "the lift direction is invariant under Q -> cQ + d "
            "(DZ_{cQ+d}[g] = DZ_Q[g]/c), so the angle is chart-exact."
        ),
        "exact_lab_deltaZ_site_L4_minus_L3": dz_site,
        "angles": {},
        "sensitivity_equal_weighting": {},
        "historical_withdrawn_reading": None,
        "scope_limits": [
            "The exact laboratory is deterministic: no sampling covariance for "
            "Delta Z_exact; the N725 jackknife covariance lives in "
            "n725-zflow-corrected.json and is cited, not recomputed.",
            "Tiny-torus / exact-lab Delta Z cannot prove an asymptotic mismatch.",
            "Equal weighting is a sensitivity, not a second vote.",
        ],
    }

    for gname, gv in (("spin0", g), ("equal", g_equal)):
        lift3 = dz_g_lift(q3, gv)
        lift4 = dz_g_lift(q4, gv)
        liftb = dz_g_lift(qb, gv)
        out["angles"][gname] = {
            "DZ_Q3_vs_deltaZ_site": {
                "DZ_Q[g]": lift3,
                "Q_used": "site L=3 quantile on deciles (from census Z, affine-exact)",
                "angle": unoriented_angle_deg(lift3, dz_site),
            },
            "DZ_Q4_vs_deltaZ_site": {
                "DZ_Q[g]": lift4,
                "Q_used": "site L=4 quantile on deciles (from census Z, affine-exact)",
                "angle": unoriented_angle_deg(lift4, dz_site),
            },
            "DZ_Qbond_vs_deltaZ_site": {
                "DZ_Q[g]": liftb,
                "Q_used": "bond L=3 quantile on deciles (from census Z, affine-exact)",
                "angle": unoriented_angle_deg(liftb, dz_site),
            },
        }

    # labelled sensitivity: spin0 vs equal, primary comparison only
    out["sensitivity_equal_weighting"] = {
        "note": "Equal weighting is a sensitivity, not a second vote.",
        "delta_unoriented_deg_primary": out["angles"]["spin0"]["DZ_Q3_vs_deltaZ_site"]["angle"]["unoriented_angle_deg"]
        - out["angles"]["equal"]["DZ_Q3_vs_deltaZ_site"]["angle"]["unoriented_angle_deg"],
    }

    if GLIFT_HIST.exists():
        gl = json.loads(GLIFT_HIST.read_text())
        d6 = gl["direction6"]
        out["historical_withdrawn_reading"] = {
            "source": "results/probe-invariant-shape/blindness-and-glift.json",
            "raw_deg": d6["angle_raw_deg"],
            "affine_removed_deg": d6["angle_after_affine_removal_deg"],
            "unoriented_deg": min(d6["angle_after_affine_removal_deg"] % 180.0,
                                  180.0 - (d6["angle_after_affine_removal_deg"] % 180.0)),
            "status": "withdrawn in #633 (raw oriented angle and affine-removed reading replaced by unoriented 20.4)",
        }

    if N725_CORR.exists():
        n725 = json.loads(N725_CORR.read_text())
        w = n725["weightings"]["spin0"]
        out["n725_cited_not_recomputed"] = {
            "file": "results/probe-invariant-shape/n725-zflow-corrected.json",
            "g_lift_DZ_Q_anchors_0.2_0.8": w["g_lift_DZ_Q_anchors_0.2_0.8"],
            "angle_g_lift_vs_Z_in_0.2_0.8_chart": w["angle_g_lift_vs_Z_in_0.2_0.8_chart"],
            "tiny_torus_unoriented_deg": n725["tiny_torus_g_vs_deltaZ"]["unoriented_deg"],
            "status": "no second production size; N725 g-vs-DeltaZ remains undetermined",
        }

    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out["angles"]["spin0"]["DZ_Q3_vs_deltaZ_site"]["angle"], indent=1))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
