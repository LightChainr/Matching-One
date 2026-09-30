#!/usr/bin/env python3
"""Exact L3 two-time uniform-label kernels from archived 9! permutation tables.

Uses SymPy for polynomial arithmetic, not new enumeration or simulations.
"""
from pathlib import Path
import gzip
import hashlib
import json
import math
import platform
import time
import sympy as S


HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "birth-gap-20260929/data"


def main():
    started = time.perf_counter()
    source = json.loads((DATA / "run.json").read_text())
    p, q = S.symbols("p q")
    a, b, c = S.symbols("a b c")
    A = lambda x: 2*x**3-6*x**2+3*x
    B = lambda x: -2*x**3+3*x+1
    records = []
    for ref in source["controls"]:
        raw = (DATA / ref["file"]).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == ref["sha256"]
        batch = json.loads(gzip.decompress(raw))
        n, total = batch["N"], batch["samples"]
        assert n == 9 and total == math.factorial(n) and batch["mode"] == "exact"
        assert sum(row[2] for row in batch["histogram"]) == total
        def h(k, ell):
            return sum(count for i, j, count in batch["histogram"] if i <= k and j > ell)
        kernel = S.expand(sum(
            S.Rational(h(k, ell), total)*math.comb(n, k)*math.comb(n-k, ell-k) *
            p**k*(q-p)**(ell-k)*(1-q)**(n-ell)
            for k in range(n+1) for ell in range(k, n+1)))
        if batch["lattice"] == "square":
            closed_form = 6*p**3*(1-q**2)**3
            determinant_form = S.Integer(0)
        else:
            assert batch["lattice"] == "triangular"
            closed_form = 9*p**3*(1-q)**3*(A(p)+B(q))
            determinant_form = 81*a**3*b**3*(1-b)**3*(1-c)**3*(A(b)-A(a))*(B(c)-B(b))
        assert S.expand(kernel-closed_form) == 0
        substitute = lambda x, y: kernel.subs({p: x, q: y}, simultaneous=True)
        determinant = substitute(a, c)*substitute(b, b)-substitute(a, b)*substitute(b, c)
        assert S.expand(determinant-determinant_form) == 0
        coefficients = S.Matrix(n+1, n+1, lambda i, j: kernel.coeff(p, i).coeff(q, j))
        times = (S.Rational(1, 3), S.Rational(1, 2), S.Rational(2, 3))
        x, y, z = times
        Hab, Hac, Hbc, Hbb = substitute(x, y), substitute(x, z), substitute(y, z), substitute(y, y)
        difference = S.cancel(Hac/Hab-(Hbc-Hac)/(Hbb-Hab))
        records.append({
            "lattice": batch["lattice"], "L": 3, "N": n,
            "source_path": "analysis/birth-gap-20260929/data/"+ref["file"],
            "source_sha256": ref["sha256"], "permutations_in_archive": total,
            "kernel_closed_form": str(closed_form),
            "kernel_polynomial": str(kernel),
            "coefficient_matrix_rank": int(coefficients.rank()),
            "markov_determinant_closed_form": str(determinant_form),
            "markov_identity_holds_for_all_label_triples": determinant_form == 0,
            "specified_label_triple": [str(x) for x in times],
            "conditional_survival_difference_at_triple": str(difference),
            "exact_polynomial_identity": True,
        })
    result = {
        "schema": "matching-one.birth-kernel-exact-clock.v1",
        "as_of": "2026-09-29",
        "definition": "H(p,q)=P(T1<=p,T2>q), 0<=p<=q<=1, iid Uniform site labels",
        "clock_transform": "sum_(k<=ell) Multinomial(N;k,ell-k,N-ell;p,q-p,1-q) H_J(k,ell)",
        "determinant": "H(a,c)H(b,b)-H(a,b)H(b,c), a<b<c",
        "triangular_A": str(A(p)), "triangular_B": str(B(q)),
        "cases": records,
        "scope": "Exact archived L3 engine convention; rank natural filtration, ordinary deterministic-time Markov property; no continuum extrapolation or hidden-positive-state dimension claim.",
        "no_new_permutations_or_random_labels": True,
        "python": platform.python_version(), "sympy": S.__version__,
        "wall_seconds": time.perf_counter()-started,
    }
    (HERE / "exact-clock.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
