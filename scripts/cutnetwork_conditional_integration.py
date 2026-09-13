#!/usr/bin/env python3
"""Cut-network conditional integration of the original-U influence function (#578).

Contract A (influence-function diagnostic): freeze the nuisance/root quantities and
compare the empirical variance of phi vs E[phi | S] for a nested conditioning
hierarchy S0..S3. The first-order object is the conditional expectation

    m(S) = E[phi | S],

NOT a ratio. We verify the law of total variance exactly on a small, fully
specified cut network and report the variance decomposition and a wall-clock cost
model, so the issue's decision ("promote one conditioning contract, or stop") has
a reproducible basis.

IMPORTANT HONESTY NOTES (see notes/p578-cutnetwork-conditional-integration-20260913.md):
  * The exact original-U influence-function coefficients (psi_E, psi_q, a_g, R,
    R_xi, D) and the 147 frozen P334 prefixes with their 84-86% suffix-noise
    measurement are NOT in the tree. The phi used here is a faithfully-structured
    PROTOTYPE of the issue's schematic
        phi_g = [c_g psi_E,g - R a_g psi_q,g - R_xi a_g (q - <q>_g)] / D .
    The variance-decomposition PROCEDURE and the tower-identity verification are
    exact and general; the numeric phi values are illustrative of the method.
  * The physical continuation law used is the fixed-cardinality two-terminal
    VERTEX reliability law ("choose uniformly among remaining switchable
    vertices"), i.e. every subset of a fixed occupied-cardinality k is equally
    likely. This is the law #487/#491 specify; we do NOT substitute independent
    edge reliability (no exact transform is asserted).
  * No nonlinear ratio/root is Rao-Blackwellized here; only the first-order
    influence function m(S) = E[phi|S] is integrated.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

OUT = Path("results/cutnetwork-conditional-20260913")
DEROUT = OUT / "derived"
RAWOUT = OUT / "raw"

# ---- A small planar two-terminal vertex network (terminals L,R) -------------
# switchable vertices s1..s5; edges fixed (cut geometry). Planar diamond stack.
SWITCHABLE = ["s1", "s2", "s3", "s4", "s5"]
EDGES = [
    ("L", "s1"), ("L", "s2"),
    ("s1", "s3"), ("s1", "s4"), ("s2", "s3"), ("s2", "s4"),
    ("s3", "s5"), ("s4", "s5"),
    ("s5", "R"),
]


def l_r_connected(occupied: frozenset[str]) -> bool:
    """BFS from L to terminal R over occupied switchable vertices + fixed edges.
    Terminals L and R are always present; only the switchable vertices are
    subject to the future occupation law."""
    available = occupied | {"L", "R"}
    seen = {"L"}
    stack = ["L"]
    while stack:
        u = stack.pop()
        for a, b in EDGES:
            if a == u and b in available and b not in seen:
                seen.add(b); stack.append(b)
            if b == u and a in available and a not in seen:
                seen.add(a); stack.append(a)
    return "R" in seen


def connectivity_indicator(occupied: frozenset[str]) -> int:
    return int(l_r_connected(occupied))


def phi_prototype(occupied: frozenset[str]) -> float:
    """Faithful structural prototype of the original-U influence function.

    phi_g = [c*psi_E - R*a*psi_q - R_xi*a*(q - qbar)] / D
    Here psi_E = L-R connectivity indicator (the "E" estimator),
    psi_q = q-proxy = fraction of switchable vertices occupied,
    a, c, R, R_xi, D are frozen nuisance/root constants.
    """
    c, a, R, R_xi, D = 1.0, 0.7, 0.5, 0.3, 1.2
    qbar = 0.5
    psi_E = connectivity_indicator(occupied)
    q = len(occupied) / len(SWITCHABLE)
    return (c * psi_E - R * a * q - R_xi * a * (q - qbar)) / D


# ---- enumerate the fixed-cardinality future law ------------------------------
def all_futures():
    """Every (k, subset-of-size-k) under the fixed-cardinality vertex law;
    uniform over subsets of each k (equal weight per cardinality bucket)."""
    futures = []
    for k in range(len(SWITCHABLE) + 1):
        for combo in itertools.combinations(SWITCHABLE, k):
            futures.append(frozenset(combo))
    return futures


def empirical_variance(values: list[float]) -> float:
    a = np.array(values, dtype=float)
    return float(a.var(ddof=0))


# ---- conditioning states S0..S3 ---------------------------------------------
def state_of(level: int, occ: frozenset[str]) -> tuple:
    if level == 0:
        return ("ALL",)                       # S0: nothing conditioned
    if level == 1:
        return ("k", len(occ))                # S1: survival cardinality
    if level == 2:
        # S2: low-cost cut-network invariant = (cardinality, s1 state, s5 state)
        return ("inv", len(occ), int("s1" in occ), int("s5" in occ))
    # S3: full typed cut network = exact occupied subset
    return ("full", tuple(sorted(occ)))


def variance_decomposition(level: int, futures, phi):
    groups: dict = {}
    for occ in futures:
        s = state_of(level, occ)
        groups.setdefault(s, []).append(phi(occ))
    # weight each conditioning state by its probability (group size / total)
    total_n = sum(len(v) for v in groups.values())
    means = {s: float(np.mean(v)) for s, v in groups.items()}
    vars_ = {s: empirical_variance(v) for s, v in groups.items()}
    p = {s: len(v) / total_n for s, v in groups.items()}
    grand = sum(p[s] * means[s] for s in groups)
    between = sum(p[s] * (means[s] - grand) ** 2 for s in groups)
    within = sum(p[s] * vars_[s] for s in groups)
    total = between + within
    n_states = len(groups)
    inner_cost = sum(len(v) for v in groups.values())  # total future evals
    return {
        "level": level,
        "n_states": n_states,
        "V_between": between,
        "V_within": within,
        "V_total": total,
        "var_reduction_fraction": (between / total) if total > 0 else 0.0,
        "inner_cost": inner_cost,
        "outer_cost_states": n_states,
        "tower_residual": abs(total - empirical_variance([phi(f) for f in futures])),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    DEROUT.mkdir(parents=True, exist_ok=True)
    RAWOUT.mkdir(parents=True, exist_ok=True)

    futures = all_futures()
    phi_vals = [phi_prototype(f) for f in futures]
    total_var = empirical_variance(phi_vals)

    decomps = [variance_decomposition(lvl, futures, phi_prototype) for lvl in range(4)]
    # wall-clock proxy: cost = outer_states + inner_future_evals (each eval is one
    # configuration; conditional integration cost also includes the S-identification
    # which we proxy by outer_states). information per wall-clock = var_reduction/cost.
    for d in decomps:
        cost = d["outer_cost_states"] + d["inner_cost"]
        d["wallclock_cost_proxy"] = cost
        d["info_per_wallclock"] = d["var_reduction_fraction"] / cost

    # tiny exact control: tower identity must hold to machine precision
    tower_ok = all(d["tower_residual"] < 1e-9 for d in decomps)

    result = {
        "issue": 578,
        "contract": "A (influence-function diagnostic)",
        "network": {
            "terminals": ["L", "R"],
            "switchable": SWITCHABLE,
            "n_futures": len(futures),
            "fixed_cardinality_law": "uniform over subsets of each k",
        },
        "total_var_phi": total_var,
        "tower_identity_holds_exactly": tower_ok,
        "decomposition_S0_S3": decomps,
        "decision_input": {
            "best_var_reduction_level": max(decomps, key=lambda d: d["var_reduction_fraction"])["level"],
            "note": "S3 (full typed network) gives V_within=0 by construction; coarser "
                    "levels trade variance reduction against inner cost.",
        },
        "buy_back": {
            "real_original_U_phi_coefficients": "NOT in tree (psi_E, psi_q, a_g, R, R_xi, D) -> declared buy-back",
            "P334_147_saved_prefixes": "NOT in tree -> saved-prefix asset audit BLOCKED",
            "suffix_noise_84_86pct": "cannot be re-attributed to true phi without the archive -> buy-back",
        },
    }
    DEROUT.joinpath("cutnetwork_conditional.json").write_text(json.dumps(result, indent=2))
    RAWOUT.joinpath("phi_per_future.json").write_text(
        json.dumps([{"future": sorted(f), "phi": phi_prototype(f)} for f in futures], indent=2))
    meta = {
        "issue": 578,
        "contract": "A",
        "tower_identity_holds_exactly": tower_ok,
        "total_var_phi": total_var,
    }
    OUT.joinpath("metadata.json").write_text(json.dumps(meta, indent=2))
    OUT.joinpath("commands.txt").write_text(
        "PY=/Users/lc/.workbuddy/binaries/python/envs/default/bin/python\n"
        "git checkout analysis/p578-cutnetwork-conditional-integration-20260913\n"
        "$PY scripts/cutnetwork_conditional_integration.py\n")

    print(f"n_futures={len(futures)} total_var_phi={total_var:.6f} tower_ok={tower_ok}")
    for d in decomps:
        print(f"  S{d['level']}: n_states={d['n_states']:>3d} V_between={d['V_between']:.4f} "
              f"V_within={d['V_within']:.4f} var_red={d['var_reduction_fraction']:.3f} "
              f"cost={d['wallclock_cost_proxy']} info/wall={d['info_per_wallclock']:.2e}")


if __name__ == "__main__":
    main()
