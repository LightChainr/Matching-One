#!/usr/bin/env python3
"""Program E: an explicit low-rank / high-nonnegative-rank response family.

Construction (abstract experiment algebra, the shape a native branching
realisation must eventually instantiate):

  * hidden states x_0..x_{n-1} = the vertices of a regular n-gon (points in
    the plane), read as "activated micro-configurations";
  * experiments e_0..e_{n-1} = the n facets (edges) of the polygon; the
    success "probability" (a nonnegative response) of experiment e_i on state
    x_j is the slack  b_i - a_i^T x_j  >= 0  (distance of the vertex from the
    facet, scaled), i.e. a nonnegative affine readout.

The response matrix M_{i,j} = b_i - a_i^T x_j is a slack matrix.  It has

  rank(M) <= 1 + 2 = 3        (M = 1 b^T - A^T X, A and X are 2 x n)

while its nonnegative rank grows with n.  For the square (n = 4) this note
gives a fully self-contained certificate: rank = 3, nonnegative rank = 4
(explicit fooling set of size 4).  For general even n the nonnegative rank is
n (Fiorini--Rothvoss--Tiwari, extension complexity of the regular n-gon);
what matters here is only that it is unbounded.

Contrast with #549: there the hidden coordinate is a SCALAR a, so every
nonnegative affine readout is a combination of {1, a} and the nonnegative
rank is 2 — no separation.  The exact feature that creates the separation is a
hidden state of dimension >= 2 whose nonnegative affine readouts form a slack
(facet) structure, not merely a one-parameter family.

Deterministic; no sampling.
"""

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def polygon_slack(n: int):
    ang = 2 * np.pi * np.arange(n) / n
    X = np.stack([np.cos(ang), np.sin(ang)])            # 2 x n vertices
    A = np.zeros((2, n))
    b = np.zeros(n)
    for i in range(n):
        mid = X[:, i] + X[:, (i + 1) % n]
        mid = mid / np.linalg.norm(mid)                 # outward unit normal
        A[:, i] = mid
        b[i] = float(mid @ X[:, i])
    M = b[None, :] - (A.T @ X).T                        # careful indexing
    # build M row-by-row for clarity:
    M = np.zeros((n, n))
    for i in range(n):
        M[i, :] = b[i] - A[:, i] @ X
    assert np.all(M >= -1e-12), "slack matrix must be nonnegative"
    M[M < 0] = 0.0
    return M, X


def rank(M):
    return int(np.linalg.matrix_rank(M, tol=1e-9))


def is_fooling_set(M, F):
    for (i, j) in F:
        if M[i, j] <= 1e-12:
            return False
    for a in range(len(F)):
        for c in range(a + 1, len(F)):
            i, j = F[a]
            ip, jp = F[c]
            if i == ip or j == jp:
                return False
            if M[i, jp] > 1e-12 and M[ip, j] > 1e-12:
                return False
    return True


def main() -> None:
    out = {}
    for n in (4, 6, 8, 10):
        M, X = polygon_slack(n)
        out[str(n)] = {"rank": rank(M), "nonneg_rank_lower": None}
        print(f"n={n}: rank = {rank(M)} (<= 3)")

    M4, _ = polygon_slack(4)
    F = [(0, 3), (1, 0), (2, 1), (3, 2)]
    out["4"]["fooling_set"] = F
    out["4"]["fooling_set_valid"] = is_fooling_set(M4, F)
    out["4"]["nonneg_rank"] = 4
    print(f"\nsquare slack matrix =\n{np.round(M4, 3).astype(int)}")
    print(f"fooling set {F} valid: {is_fooling_set(M4, F)}")
    print(f"square: ordinary rank = {rank(M4)}, nonnegative rank = 4 (exact)")
    out["summary"] = {
        "ordinary_rank": "<= 3 for all n",
        "nonnegative_rank": "4 at n=4 (self-contained); unbounded / =n for even n "
                             "(Fiorini-Rothvoss-Tiwari, extension complexity of the n-gon)",
        "the_necessary_feature": "hidden state of dim >= 2 with nonnegative affine "
                                 "(slack) readouts; #549's scalar coordinate cannot "
                                 "separate (nonneg rank = rank = 2 there)",
    }
    dest = ROOT / "results" / "cutnetwork-positive-vs-signed-slack-family" / "latest.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1, default=float))
    print("wrote", dest)


if __name__ == "__main__":
    main()
