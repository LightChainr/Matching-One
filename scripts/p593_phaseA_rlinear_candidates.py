#!/usr/bin/env python3
"""Issue #593 Phase A: which candidate definition reproduces r_linear = 10/26/72?

The #599 probe measured block-Hankel / McMillan degrees of 9/8, 4/5, 2/2 at
widths 4..6, while #593/#588 report r_linear = 10, 26, 72.  This script
computes every candidate "linear dimension" notion at widths 4..8 and checks
which one (if any) reproduces the reported triple:

  candidate 1  joint Krylov span  dim span{G0^k f} over the readout seeds
               (variants: D0 seeds, D0+constant, D2 seeds, k>=1 only, H0 powers)
  candidate 2  linear lumping degrees of freedom
               (coarsest exact strong lumping vs G0 only / G0+H jointly,
                distinct readout value tuples)
  candidate 3  pencil ranks and combinatorial structures
               (rank G0, rank H, rank J, rank D mod p; distinct rows of G0)

All modular arithmetic is exact over the prime 2147483647, matching the
repository's RANK_PRIME convention (integer vectors, no pivot tolerance).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from p398_intervention_transport import (  # noqa: E402
    DICTIONARIES,
    Generator,
    RANK_PRIME,
    exact_lumping,
    observable_reachable_dimension,
)
from noncrossing_connectivity_codec import noncrossing_states  # noqa: E402
from planar_state_operations import rgs_to_blocks  # noqa: E402

P = RANK_PRIME

OBSERVABLES = {
    "blocks": lambda s: float(max(s) + 1),
    "singletons": lambda s: float(sum(1 for b in rgs_to_blocks(s) if len(b) == 1)),
    "wrap": lambda s: 1.0 if s[0] == s[-1] else 0.0,
    "max_block": lambda s: float(max(len(b) for b in rgs_to_blocks(s))),
    "linked_pairs": lambda s: float(
        sum(len(b) * (len(b) - 1) // 2 for b in rgs_to_blocks(s))
    ),
    "boundary_span": lambda s: float(sum(max(b) - min(b) for b in rgs_to_blocks(s))),
    "halves_linked": lambda s: _halves_linked(s),
    "covering_depth": lambda s: _covering_depth(s),
}


def _halves_linked(state):
    width = len(state)
    half = width // 2
    for block in rgs_to_blocks(state):
        if any(p < half for p in block) and any(p >= half for p in block):
            return 1.0
    return 0.0


def _covering_depth(state):
    blocks = [sorted(b) for b in rgs_to_blocks(state)]
    best = 0
    for point in range(len(state)):
        depth = sum(
            1 for b in blocks if b[0] < point < b[-1] and point not in b
        )
        best = max(best, depth)
    return float(best)


def sparse_mod_rows(gen, rates):
    """Integer mod-p rows of the generator (row = source, acting on functions)."""
    rows = []
    for row in gen.rows(rates):
        rows.append([(c, int(round(v)) % P) for c, v in row])
    return rows


def matvec_mod(rows_sparse, size, vec):
    out = [0] * size
    for i, row in enumerate(rows_sparse):
        acc = 0
        for c, v in row:
            acc += v * vec[c]
        out[i] = acc % P
    return out


def rank_mod(vectors, size):
    """Exact rank of integer vectors in Z_p^size via echelon insertion."""
    pivots = {}  # pivot column -> normalized vector (numpy int64)

    def insert(vec):
        working = np.array(vec, dtype=np.int64)
        for col in sorted(pivots):
            pv = pivots[col]
            f = int(working[col])
            if f:
                working = (working - f * pv) % P
        nz = np.nonzero(working)[0]
        if nz.size == 0:
            return False
        col = int(nz[0])
        inv = pow(int(working[col]), P - 2, P)
        working = (working * inv) % P
        # back-substitute for reduced echelon (keeps future inserts cheap)
        for other, pv in list(pivots.items()):
            if int(pv[col]):
                pivots[other] = (pv - int(pv[col]) * working) % P
        pivots[col] = working
        return True

    for v in vectors:
        insert(v)
    return len(pivots)


def krylov_dimension(gen, rows_sparse, size, seed_vectors, max_depth=None):
    """dim span{G^k f} mod p, with depth tracked until 2 consecutive plateaus."""
    frontier = [[int(round(x)) % P for x in v] for v in seed_vectors]
    echelon = {}

    def insert(vec):
        working = list(vec)
        for col, pv in echelon.items():
            f = working[col]
            if f:
                working = [(a - f * b) % P for a, b in zip(working, pv)]
        col = next((i for i, a in enumerate(working) if a), None)
        if col is None:
            return None
        inv = pow(working[col], P - 2, P)
        working = [(a * inv) % P for a in working]
        for other, pv in list(echelon.items()):
            if pv[col]:
                echelon[other] = [(a - pv[col] * b) % P for a, b in zip(pv, working)]
        echelon[col] = working
        return working

    depth = 0
    history = []
    added_any = True
    while frontier and added_any:
        added = []
        for v in frontier:
            w = insert(v)
            if w is not None:
                added.append(w)
        history.append(len(echelon))
        depth += 1
        if max_depth is not None and depth >= max_depth:
            break
        if not added:
            added_any = False
            break
        frontier = [matvec_mod(rows_sparse, size, v) for v in added]
        if len(history) >= 3 and history[-1] == history[-2] == history[-3]:
            break
    return len(echelon), depth


def float_krylov_dimension(gen, rows_float, size, seed_vectors, tol=1e-8):
    """Floating-point Krylov dimension (the documented failure mode)."""
    basis = []
    frontier = [list(map(float, v)) for v in seed_vectors]
    while frontier:
        added = []
        for v in frontier:
            w = list(v)
            for b in basis:
                dot = sum(a * x for a, x in zip(b, w))
                w = [a - dot * x for a, x in zip(w, b)]
            nrm = sum(a * a for a in w) ** 0.5
            if nrm > tol:
                basis.append([a / nrm for a in w])
                added.append(w)
        if not added:
            break
        frontier = [
            [sum(v * vec[c] for c, v in row_i) for row_i in rows_float]
            for vec in added
        ]
    return len(basis)


def lumping_block_count(gen, colours, operators):
    result = exact_lumping(gen, colours, operators)
    return result["blocks"]


def analyse_width(width):
    gen = Generator(width)
    size = gen.size
    states = gen.states
    baseline = gen.baseline_rates()
    rows_sparse = sparse_mod_rows(gen, baseline)
    h_rows = sparse_mod_rows(gen, gen.coefficients("uniform_join_minus_detach"))

    # integer observation vectors
    obs = {name: [int(OBSERVABLES[name](s)) for s in states] for name in OBSERVABLES}
    D0 = ["blocks", "singletons", "wrap"]
    D2 = DICTIONARIES["D2_plus_nonlocal_topology"]

    seeds_D0 = [obs[n] for n in D0]
    seeds_D0_one = [[1] * size] + seeds_D0
    seeds_D2 = [obs[n] for n in D2]
    colour_D0 = list(zip(obs["blocks"], obs["singletons"], obs["wrap"]))

    out = {"width": width, "size": size}

    # --- candidate 1: joint Krylov spans -----------------------------------
    dim, depth = krylov_dimension(gen, rows_sparse, size, seeds_D0)
    out["c1_krylov_D0"] = dim
    out["c1_krylov_D0_depth"] = depth
    dim, _ = krylov_dimension(gen, rows_sparse, size, seeds_D0_one)
    out["c1_krylov_D0_plus_constant"] = dim
    dim, _ = krylov_dimension(gen, rows_sparse, size, seeds_D2)
    out["c1_krylov_D2"] = dim

    # k >= 1 only: apply G once to the seeds, then iterate
    shifted = [matvec_mod(rows_sparse, size, [int(x) % P for x in v]) for v in seeds_D0]
    dim, _ = krylov_dimension(gen, rows_sparse, size, shifted)
    out["c1_krylov_D0_kge1"] = dim

    # H-powered Krylov
    dim, _ = krylov_dimension(gen, h_rows, size, seeds_D0)
    out["c1_krylov_D0_with_H"] = dim

    # float failure-mode check (cheap, widths <= 6 only)
    if width <= 6:
        rows_float = gen.rows(baseline)
        out["c1_float_krylov_D0"] = float_krylov_dimension(
            gen, rows_float, size, seeds_D0
        )

    # cross-check against the repository's own implementation (width <= 7)
    if width <= 7:
        out["repo_observable_reachable_dimension"] = (
            observable_reachable_dimension(gen, seeds_D0, size)["dimension"]
        )

    # --- candidate 2: linear lumping degrees of freedom ---------------------
    out["c2_lumping_G0_only"] = lumping_block_count(gen, colour_D0, [baseline])
    out["c2_lumping_G0_and_H_joint"] = lumping_block_count(
        gen, colour_D0, [baseline, gen.coefficients("uniform_join_minus_detach")]
    )
    out["c2_distinct_readout_tuples"] = len(set(colour_D0))

    # --- candidate 3: pencil ranks / combinatorial structures ---------------
    out["c3_rank_G0_mod_p"] = _matrix_rank_mod(rows_sparse, size)
    out["c3_rank_H_mod_p"] = _matrix_rank_mod(h_rows, size)

    join_only = {m: (1.0 if m[0] == "join" else 0.0) for m in gen.moves}
    detach_only = {m: (1.0 if m[0] == "detach" else 0.0) for m in gen.moves}
    out["c3_rank_J_mod_p"] = _matrix_rank_mod(sparse_mod_rows(gen, join_only), size)
    out["c3_rank_D_mod_p"] = _matrix_rank_mod(sparse_mod_rows(gen, detach_only), size)

    row_patterns = set()
    for row in gen.rows(baseline):
        row_patterns.add(tuple((c, round(v)) for c, v in row))
    out["c3_distinct_rows_G0"] = len(row_patterns)

    return out


def _matrix_rank_mod(rows_sparse, size):
    """Rank of the generator matrix itself over Z_p (rows = sources)."""
    dense_rows = []
    for row in rows_sparse:
        v = [0] * size
        for c, val in row:
            v[c] = val
        dense_rows.append(v)
    return rank_mod(dense_rows, size)


def main():
    results = []
    for width in (4, 5, 6, 7, 8):
        t0 = time.time()
        row = analyse_width(width)
        row["seconds"] = round(time.time() - t0, 1)
        results.append(row)
        print(json.dumps(row), flush=True)
    Path("/tmp/p593-phaseA-results.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
