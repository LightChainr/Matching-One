#!/usr/bin/env python3
"""C3 selection-rule demonstration for the finite-group response theorem.

Synthetic but exact object: six states arranged as two C3-orbits A and B, with
a rotation-equivariant continuous-time generator.  The permutation
representation decomposes as V = 2*triv + 2*W (W the real 2-d irrep of C3).
A diagonal perturbation field h on A transforms in W.

Theorem being demonstrated (first order, any finite group K): with an even
(invariant) source in rho_B and readout in rho_C, a perturbation transforming
in rho_H has zero pointwise linear response whenever the trivial
representation is absent from rho_B* (x) rho_H (x) rho_C.  In particular

  triv source, triv readout, W perturbation  -> zero (2.6e-17 here)
  triv source, W readout,    W perturbation  -> generically nonzero

Run: python3 scripts/probe/group_selection_demo.py
"""

import json
import math
import sys
from pathlib import Path

import numpy as np

n = 6
IDX = lambda o, i: o * 3 + i  # noqa: E731


def generator_matrix() -> np.ndarray:
    G = np.zeros((n, n))
    for i in range(3):
        G[IDX(0, i), IDX(0, (i + 1) % 3)] = 1.0   # a_i -> a_{i+1}
        G[IDX(0, i), IDX(1, i)] = 2.0             # a_i -> b_i
        G[IDX(1, i), IDX(1, (i - 1) % 3)] = 1.0   # b_i -> b_{i-1}
        G[IDX(1, i), IDX(0, i)] = 1.0             # b_i -> a_i
    for i in range(n):
        G[i, i] = -G[i].sum()
    return G


def perm_matrix(k: int) -> np.ndarray:
    Pk = np.zeros((n, n))
    for o in (0, 1):
        for i in range(3):
            Pk[IDX(o, (i + k) % 3), IDX(o, i)] = 1
    return Pk


def expv(v: np.ndarray, t: float) -> np.ndarray:
    w, V = np.linalg.eig(G)
    return np.real((V @ np.diag(np.exp(w * t)) @ np.linalg.inv(V)) @ v)


def integrand(mu: np.ndarray, H: np.ndarray, f: np.ndarray, t: float, s: float) -> float:
    """mu^T e^{(t-s)G} H e^{sG} f."""
    return float(mu @ expv(H @ expv(f, s), t - s))


G = generator_matrix()
assert np.allclose(G.sum(axis=1), 0)
for k in range(3):
    P = perm_matrix(k)
    assert np.allclose(P @ G @ P.T, G, atol=1e-12), "not C3-equivariant"

# isotypic projectors
Ptr = sum(perm_matrix(k) for k in range(3)) / 3
ang = 2 * math.pi / 3
w = np.zeros((n, 4))
for i in range(3):
    w[IDX(0, i), 0] = math.sqrt(2 / 3) * math.cos(i * ang)
    w[IDX(0, i), 1] = math.sqrt(2 / 3) * math.sin(i * ang)
    w[IDX(1, i), 2] = math.sqrt(2 / 3) * math.cos(i * ang)
    w[IDX(1, i), 3] = math.sqrt(2 / 3) * math.sin(i * ang)
PW = w @ w.T
assert np.allclose(Ptr @ PW, 0, atol=1e-12)
dim_triv, dim_W = round(np.trace(Ptr)), round(np.trace(PW))

# vectors
mu_triv = np.ones(n) / n
mu_W = np.zeros(n)
f_W = np.zeros(n)
h_W = np.zeros(n)
for i in range(3):
    mu_W[IDX(0, i)] = math.sqrt(2 / 3) * math.cos(i * ang)
    f_W[IDX(0, i)] = math.sqrt(2 / 3) * math.sin(i * ang)
    h_W[IDX(0, i)] = math.sqrt(2 / 3) * math.cos(i * ang)
mu_W /= np.linalg.norm(mu_W)
f_W /= np.linalg.norm(f_W)
f_triv = np.array([1.0, 1, 1, -1, -1, -1])
f_triv /= np.linalg.norm(f_triv)
Hw = np.diag(h_W)  # W-charged diagonal field (0 on B)

cases = [
    ("triv src / triv readout", mu_triv, f_triv, "ZERO expected: triv not in W (x) triv"),
    ("triv src / W readout", mu_triv, f_W, "nonzero expected: triv in W (x) W"),
    ("W src / triv readout", mu_W, f_triv, "nonzero expected: triv in W (x) W"),
    ("W src / W readout", mu_W, f_W, "nonzero expected (generic)"),
]
results = []
for name, mu, f, expectation in cases:
    vals = [
        integrand(mu, Hw, f, t, s)
        for t in (0.5, 1.0, 2.0)
        for s in (0.25 * t, 0.5 * t, 0.75 * t)
    ]
    peak = max(abs(v) for v in vals)
    results.append({"case": name, "expectation": expectation,
                    "peak_abs_integrand": peak})
    print(f"{name:26s} peak |integrand| = {peak:.3e}   [{expectation}]")

def expv_free(A: np.ndarray, v: np.ndarray, t: float) -> np.ndarray:
    """exp(t A) v for an arbitrary (not necessarily generator) A."""
    w, V = np.linalg.eig(A)
    return np.real((V @ np.diag(np.exp(w * t)) @ np.linalg.inv(V)) @ v)


def response(mu, f, A, t):
    return float(mu @ expv_free(A, f, t))


# second-order re-entry of the W direction into the unmarked task
mu, f = mu_triv, f_triv
t = 1.0
r0 = response(mu, f, G, t)
eps_ladder = [0.02, 0.04, 0.08]
diffs_w = [abs(response(mu, f, G + e * Hw, t) - r0) for e in eps_ladder]
slope_w = math.log(diffs_w[-1] / diffs_w[0]) / math.log(eps_ladder[-1] / eps_ladder[0])
h_triv = np.ones(n)
Ht = np.diag(h_triv)
diffs_t = [abs(response(mu, f, G + e * Ht, t) - r0) for e in eps_ladder]
slope_t = math.log(diffs_t[-1] / diffs_t[0]) / math.log(eps_ladder[-1] / eps_ladder[0])
print(f"triv/triv: W-perturbation second-order slope = {slope_w:.3f} "
      f"(predicted 2: W(x)W contains triv)")
print(f"triv/triv: triv-perturbation first-order slope = {slope_t:.3f} (control)")

payload = {
    "schema": "matching-one.probe-group-selection.v1",
    "probe_issue": 599,
    "group": "C3",
    "state_count": n,
    "isotypic_decomposition": {"triv": dim_triv, "W": dim_W},
    "perturbation_representation": "W",
    "cases": results,
    "second_order": {
        "task": "triv source / triv readout",
        "eps_ladder": eps_ladder,
        "W_perturbation_abs_diffs": diffs_w,
        "W_perturbation_loglog_slope": slope_w,
        "triv_perturbation_loglog_slope_control": slope_t,
    },
}
dest = Path(__file__).resolve().parents[2] / "results" / "probe-group-selection" / "latest.json"
dest.parent.mkdir(parents=True, exist_ok=True)
dest.write_text(json.dumps(payload, indent=2) + "\n")
print("wrote", dest)
