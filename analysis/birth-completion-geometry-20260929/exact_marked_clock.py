#!/usr/bin/env python3
"""Exact iid-label kernels of the persistent rank-one direction sectors.

Transform the new 512-subset/prefix-DP output, without sampling labels or
rerunning permutation enumeration. SymPy performs rational polynomial work.
"""
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import sympy as S


HERE = Path(__file__).resolve().parent


def main():
    started = time.perf_counter()
    raw = (HERE / "result.json").read_bytes()
    source = json.loads(raw)
    p, q = S.symbols("p q")
    a, b, c = S.symbols("a b c")
    A = lambda x: 2*x**3-6*x**2+3*x
    B = lambda x: -2*x**3+3*x+1
    cases = []
    for case in source["cases"]:
        n = case["N"]
        assert n == 9
        total = math.factorial(n)
        pooled = S.Integer(0)
        sectors = []
        for sector in case["direction_kernels"]:
            hist = sector["paired_histogram"]
            direction = sector["direction"]
            def h(k, ell):
                return sum(count for j1, j2, count in hist if j1 <= k and j2 > ell)
            kernel = S.expand(sum(
                S.Rational(h(k, ell), total)*math.comb(n, k)*math.comb(n-k, ell-k)
                * p**k*(q-p)**(ell-k)*(1-q)**(n-ell)
                for k in range(n+1) for ell in range(k, n+1)))
            pooled += kernel
            if case["lattice"] == "square":
                if 0 in direction:
                    orbit = "axial"
                    closed = 3*p**3*(1-q)**3*((1+q)**3-p**3)
                    det_closed = -9*a**3*b**3*(1-b)**3*(1-c)**3*(b**3-a**3)*((1+c)**3-(1+b)**3)
                    mean_nu = 3-3*(1-q)*(1+q)**2/((1+q)**3-p**3)
                else:
                    orbit = "diagonal"
                    closed = 3*p**6*(1-q)**3
                    det_closed = S.Integer(0)
                    mean_nu = S.Integer(3)
            else:
                assert case["lattice"] == "triangular"
                orbit = "three_direction_orbit"
                closed = 3*p**3*(1-q)**3*(A(p)+B(q))
                det_closed = 9*a**3*b**3*(1-b)**3*(1-c)**3*(A(b)-A(a))*(B(c)-B(b))
                mean_nu = 3-(1-q)*(3-6*q**2)/(A(p)+B(q))
            assert S.expand(kernel-closed) == 0
            H = lambda x, y: kernel.subs({p: x, q: y}, simultaneous=True)
            determinant = H(a, c)*H(b, b)-H(a, b)*H(b, c)
            assert S.expand(determinant-det_closed) == 0
            # This uses the general finite completion-hazard identity, not
            # a second independent enumeration of the marked geometry.
            assert S.cancel(-(1-q)*S.diff(kernel, q)/kernel-mean_nu) == 0
            times = [S.Rational(1, 3), S.Rational(1, 2), S.Rational(2, 3)]
            x, y, z = times
            Hab, Hac, Hbb, Hbc = H(x, y), H(x, z), H(y, y), H(y, z)
            assert Hab > 0 and Hbb > Hab
            difference = S.cancel(Hac/Hab-(Hbc-Hac)/(Hbb-Hab))
            coefficients = S.Matrix(n+1, n+1, lambda i, j: kernel.coeff(p, i).coeff(q, j))
            sectors.append({
                "direction": direction, "symmetry_orbit": orbit,
                "kernel_closed_form": str(closed),
                "coefficient_matrix_rank": int(coefficients.rank()),
                "markov_determinant": str(det_closed),
                "markov_identity_for_all_label_triples": det_closed == 0,
                "cohort_mean_completion_count": str(mean_nu),
                "specified_label_triple": [str(t) for t in times],
                "early_minus_late_survival_at_triple": str(difference),
                "all_closed_form_identities_exact": True,
            })
        pooled_closed = (6*p**3*(1-q**2)**3 if case["lattice"] == "square"
                         else 9*p**3*(1-q)**3*(A(p)+B(q)))
        assert S.expand(pooled-pooled_closed) == 0
        cases.append({"lattice": case["lattice"], "L": 3, "N": n,
                      "unmarked_kernel": str(pooled_closed),
                      "marked_kernels_sum_to_unmarked": True, "sectors": sectors})
    output = {
        "schema": "matching-one.birth-marked-clock.v1",
        "as_of": "2026-09-29",
        "source": "analysis/birth-completion-geometry-20260929/result.json",
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "definition": "H_d(p,q)=P(T1<=p,T2>q,D=d), unconditional sector mass; 0<=p<=q<=1",
        "normalization": "full 9! permutation count, not conditional on visiting the sector",
        "clock_transform": "multinomial thinning of the same uniform permutation, iid Uniform site labels",
        "determinant": "H_d(a,c)H_d(b,b)-H_d(a,b)H_d(b,c)",
        "cohort_mean_completion_count": "-(1-q) partial_q log H_d(p,q), fixed entry cutoff p",
        "cases": cases,
        "boundaries": [
            "Exact L3 statement about a specified observer and clock, not a continuum result.",
            "Revealing a persistent direction can destroy a weak Markov projection; no contradiction with unmarked rank Markovness.",
            "Polynomial coefficient rank is not a positive hidden-state dimension.",
            "Four-state count-clock realization does not establish four-state label-clock closure with hidden count.",
        ],
        "python": platform.python_version(), "sympy": S.__version__,
        "wall_seconds": time.perf_counter()-started,
    }
    (HERE / "marked-clock.json").write_text(json.dumps(output, indent=2)+"\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
