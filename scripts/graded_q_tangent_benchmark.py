#!/usr/bin/env python3
"""Graded Q->1 tangent benchmark for issue #586 (positive-control-first).

Reuses the exact finite-volume gate of #581 (scripts/qtangent_scale_decomposition)
for the LATTICE control, and adds the Phase A CONTINUUM internal check:
differentiate the BPZ/fusion (2,1)-degenerate ODE in the coupling and verify the
Q-tangent satisfies the resulting inhomogeneous linear ODE.

This is a *graded tangent table*, not a Jordan yes/no verdict. The two
measure-score contributions -- B_even (duality-even Betti) and X (ambient
homology) -- are reported SEPARATELY and only summed afterwards.

KEY HONESTY NOTES (see notes/p586-graded-q-tangent-20260913.md):
  * The continuum check uses a representative (2,1)-degenerate BPZ ODE
    (hypergeometric form) to demonstrate the *method* of differentiating the
    equation and verifying the tangent against the inhomogeneous ODE. The exact
    closed-form Cai arXiv:2603.28161 boundary four-point is a buy-back; this
    script does NOT claim to reproduce that specific closed form.
  * Phase B's precise generic-Q bulk-logarithmic pair numeric (VJS / Camia-Feng)
    is not in the tree; it is declared a buy-back. The graded *method* and the
    exact L=2,3 lattice graded table are delivered.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

try:
    from scripts.qtangent_scale_decomposition import (
        OBSERVABLES,
        covariance_split,
        enumerate_torus,
        exact_identity_gate,
    )
except ModuleNotFoundError:  # pragma: no cover
    from qtangent_scale_decomposition import (  # type: ignore
        OBSERVABLES,
        covariance_split,
        enumerate_torus,
        exact_identity_gate,
    )

OUT = Path("results/graded-q-tangent-20260913")
DEROUT = OUT / "derived"
RAWOUT = OUT / "raw"


# --------------------------------------------------------------------------
# Phase A: continuum internal check -- differentiate the BPZ/(2,1) ODE in h.
# --------------------------------------------------------------------------
def bpz_solve(h: float, lam: np.ndarray) -> np.ndarray:
    """Solve the representative (2,1)-degenerate BPZ hypergeometric ODE

        L(h) G = lam*(1-lam) G'' + (1-2*lam) G' - h*(1-h) G = 0

    with boundary conditions G(0)=0, G(1)=1 (a boundary 4-point normalisation).
    Returns G evaluated on the grid `lam` (interior finite-difference solve).
    """
    n = len(lam)
    dlam = lam[1] - lam[0]
    # Build tridiagonal system A G = rhs for interior points 1..n-2.
    A = np.zeros((n, n))
    rhs = np.zeros(n)
    for i in range(1, n - 1):
        lp = lam[i]
        c = lp * (1 - lp)
        a = c / dlam**2 - (1 - 2 * lp) / (2 * dlam)       # coeff of G[i-1]
        b = c / dlam**2 + (1 - 2 * lp) / (2 * dlam)       # coeff of G[i+1]
        d = -2 * c / dlam**2 - h * (1 - h)                # coeff of G[i]
        A[i, i - 1] = a
        A[i, i] = d
        A[i, i + 1] = b
    # boundary conditions
    A[0, 0] = 1.0
    rhs[0] = 0.0
    A[n - 1, n - 1] = 1.0
    rhs[n - 1] = 1.0
    G = np.linalg.solve(A, rhs)
    return G


def bpz_continuum_check(h0: float = 0.0, eps: float = 1e-4, n: int = 4001) -> dict:
    """Verify the Q-tangent satisfies the inhomogeneous ODE.

    L(h) G = 0  =>  differentiating in h:
        L(h) T + (dL/dh) G = 0,   T = dG/dh,
    where dL/dh multiplies G by -(1-2h).  We compute G(h0+/-eps) by solving the
    homogeneous ODE, finite-difference T, evaluate the residual
        R = L(h0) T + (1-2h0) G(h0)
    on the interior grid, and report max|R|.  R ~ 0 confirms the tangent obeys
    the differentiated equation -- the issue's internal continuum check.
    """
    lam = np.linspace(0.0, 1.0, n)
    dlam = lam[1] - lam[0]
    Gp = bpz_solve(h0 + eps, lam)
    Gm = bpz_solve(h0 - eps, lam)
    G0 = bpz_solve(h0, lam)
    T = (Gp - Gm) / (2 * eps)
    residual = np.zeros(n)
    for i in range(1, n - 1):
        lp = lam[i]
        c = lp * (1 - lp)
        LhT = (c * (T[i + 1] - 2 * T[i] + T[i - 1]) / dlam**2
               + (1 - 2 * lp) * (T[i + 1] - T[i - 1]) / (2 * dlam)
               - h0 * (1 - h0) * T[i])
        residual[i] = LhT - (1 - 2 * h0) * G0[i]
    # normalise residual by typical G magnitude
    scale = max(abs(G0.max()), abs(G0.min()), 1e-12)
    return {
        "h0": h0,
        "eps": eps,
        "n_grid": n,
        "max_residual_l2_normed": float(np.max(np.abs(residual)) / scale),
        "max_residual_raw": float(np.max(np.abs(residual))),
        "passes": bool(np.max(np.abs(residual)) / scale < 1e-3),
    }


# --------------------------------------------------------------------------
# Phase A (lattice) + Phase B: exact graded tangent table via #581 gate.
# --------------------------------------------------------------------------
def graded_tangent_table() -> dict:
    out = {"lattice_control": {}, "identity_gate": {}, "observables": list(OBSERVABLES)}
    for L in (2, 3):
        table = enumerate_torus(L)
        gate = exact_identity_gate(table)
        split = covariance_split(table)
        # keep B_even and X SEPARATE (never recombined before reporting)
        rows = []
        for row in split:
            rows.append({
                "observable": row["observable"],
                "cov_with_q_score": row["cov_with_q_score"],
                "betti_even_piece": row["betti_even_piece"],
                "ambient_homology_piece": row["ambient_homology_piece"],
                "split_is_exact": row["split_is_exact"],
                "topological_fraction": row["topological_fraction"],
            })
        out["lattice_control"][f"L{L}"] = {
            "configurations": gate["configurations"],
            "split_rows": rows,
        }
        out["identity_gate"][f"L{L}"] = gate
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    DEROUT.mkdir(parents=True, exist_ok=True)
    RAWOUT.mkdir(parents=True, exist_ok=True)

    continuum = bpz_continuum_check()
    table = graded_tangent_table()

    # Phase B proxy: apply the same graded split to a chosen pair of observables
    # whose difference is a candidate bulk-logarithmic signature.  We report the
    # graded signature (B_even vs X loading) exactly; the precise generic-Q
    # LCFT pair numeric is declared a buy-back (see notes).
    result = {
        "issue": 586,
        "phase_A_continuum_check": continuum,
        "phase_A_and_B_lattice_graded_table": table,
        "notes": {
            "continuum_model": "representative (2,1)-degenerate BPZ hypergeometric ODE; demonstrates the differentiate-the-equation method. Exact Cai closed form is a buy-back.",
            "phase_B_bulk_log_numeric": "NOT in tree (VJS / Camia-Feng generic-Q data absent); declared buy-back. Graded method + exact L=2,3 table delivered.",
        },
    }
    DEROUT.joinpath("graded_q_tangent.json").write_text(json.dumps(result, indent=2))

    # raw copy of the lattice split for transparency
    RAWOUT.joinpath("lattice_graded_table.json").write_text(json.dumps(table, indent=2))

    meta = {
        "issue": 586,
        "method": "reuse #581 exact gate (enumerate_torus + covariance_split) + representative BPZ ODE tangent check",
        "continuum_check_passes": continuum["passes"],
        "max_residual_normed": continuum["max_residual_l2_normed"],
    }
    OUT.joinpath("metadata.json").write_text(json.dumps(meta, indent=2))
    commands = [
        "PY=/Users/lc/.workbuddy/binaries/python/envs/default/bin/python",
        "git checkout theory/p586-graded-q-tangent-20260913",
        "$PY scripts/graded_q_tangent_benchmark.py",
    ]
    OUT.joinpath("commands.txt").write_text("\n".join(commands) + "\n")

    print("Phase A continuum check passes:", continuum["passes"],
          "max|R|/||G||:", continuum["max_residual_l2_normed"])
    for L in (2, 3):
        print(f"L{L}: {table['lattice_control'][f'L{L}']['configurations']} configs, "
              f"identity failures: {table['identity_gate'][f'L{L}']['failures']}")


if __name__ == "__main__":
    main()
