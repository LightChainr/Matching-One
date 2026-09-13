#!/usr/bin/env python3
"""Issue #593 Phase B Deliverable 2 (r_linear): true ``dim span{G0^k f}`` at
widths 8, 9, 10 by exact elimination over the prime 2147483647.

Same mathematics as the repository's ``observable_reachable_dimension``
(block Arnoldi with deflation, integer arithmetic mod the Mersenne prime),
with the echelon stored as numpy int64 rows so width 10 is affordable.  The
Krylov depth runs until the dimension is unchanged for 3 consecutive levels,
which is the same saturation criterion the repository's frontier loop uses.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))

from p593_noncrossing_fast import WideGenerator, noncrossing_states_fast

P = 2147483647
OUTPUT = Path(__file__).resolve().parents[1] / "results" / "p593-width9-10-memory-degree" / "raw"


def _blocks(state):
    out = {}
    for point, label in enumerate(state):
        out.setdefault(label, []).append(point)
    return out.values()


def observable_vectors(states):
    size = len(states)
    blocks = [max(s) + 1 for s in states]
    singletons = [sum(1 for b in _blocks(s) if len(b) == 1) for s in states]
    wrap = [1 if s[0] == s[-1] else 0 for s in states]
    one = [1] * size
    return one, [blocks, singletons, wrap]


def sparse_mod_matrix(generator):
    """G0 as a scipy sparse int64 matrix acting on functions (row = source)."""
    data, indices, indptr = [0], [0], [0]
    rows = generator.rows(generator.baseline_rates())
    for row in rows:
        for column, value in row:
            indices.append(column)
            data.append(int(round(value)) % P)
        indptr.append(len(indices))
    return sp.csr_matrix(
        (np.array(data, dtype=np.int64), np.array(indices, dtype=np.int64),
         np.array(indptr, dtype=np.int64)),
        shape=(generator.size, generator.size),
    )


def krylov_dimension_numpy(sparse, size, seeds, max_depth=100000, budget=None):
    """dim span{G^k f} mod p; numpy echelon, frontier loop, 3-plateau stop."""

    limit = size if budget is None else min(budget, size)

    def matvec(vec):
        return sparse.dot(vec) % P

    pivots = {}  # pivot column -> numpy int64 vector (row echelon, not reduced)

    def insert(vec):
        working = np.array(vec, dtype=np.int64)
        for col in sorted(pivots):
            f = int(working[col])
            if f:
                working = (working - f * pivots[col]) % P
        nz = np.nonzero(working)[0]
        if nz.size == 0:
            return False
        col = int(nz[0])
        inv = pow(int(working[col]), P - 2, P)
        pivots[col] = (working * inv) % P
        return True

    frontier = [np.array([int(x) % P for x in v], dtype=np.int64) for v in seeds]
    history = []
    while frontier:
        added = []
        for v in frontier:
            if len(pivots) >= limit:
                break
            if insert(v):
                added.append(v)
        history.append(len(pivots))
        if not added or len(pivots) >= limit:
            break
        frontier = [matvec(v) for v in added]
        if len(history) >= 3 and history[-1] == history[-2] == history[-3]:
            break
    return len(pivots), len(history)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    widths = [int(a) for a in sys.argv[1:]] or [8, 9, 10]
    for width in widths:
        start = time.time()
        generator = WideGenerator(width, noncrossing_states_fast(width))
        sparse = sparse_mod_matrix(generator)
        one, seeds = observable_vectors(generator.states)
        dim, depth = krylov_dimension_numpy(sparse, generator.size, seeds)
        result = {
            "width": width,
            "states": generator.size,
            "r_linear_D0_true": dim,
            "krylov_depth": depth,
            "arithmetic": "exact modulo the prime 2147483647",
            "is_a_certified_lower_bound_on_the_rational_rank": True,
            "seeds": ["blocks", "singletons", "wrap"],
            "seconds": round(time.time() - start, 1),
        }
        (OUTPUT / f"rlinear_width{width}.json").write_text(json.dumps(result, indent=2))
        print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
