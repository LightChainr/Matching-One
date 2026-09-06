#!/usr/bin/env python3
"""Monte Carlo empirical control for #581: do the two typed Q-score pieces keep
distinct scale / noise profiles beyond the tiny enumerated tori?

#581's exact gate (``qtangent_scale_decomposition.py``) ran to exhaustion on the
L=2 and L=3 square-bond tori and found that the duality-even Betti tangent
``B_even`` and the ambient-homology source ``X = r - 1`` separate at Boolean
degree one by parity -- every degree-one Walsh coefficient of ``B_even`` is
exactly zero and every degree-one coefficient of ``X`` is one value.  That gate
is deterministic and adds no evidence.  The owner's roadmap note asks for the
next step: **one modest exact-critical square-bond block** asking only whether
``B_even`` and ``X`` retain measurably different finite-size scale/noise
profiles beyond L=2,3, with raw covariance vectors primary and one frozen
supported readout.

This module is that block.  It samples square-bond configurations at ``p = 1/2``
on L=3 (a correctness control against the enumerated gate), L=4 and L=8, and
decomposes the readout's covariance with ``B_even`` and with ``X`` over a small
number of predeclared Euclidean scales.  L=4 has 2**32 configurations and L=8
2**128, so the exact conditional-replica machinery cannot be used; the module
uses an **unbiased paired estimator** instead, described below.

## The paired estimator (why it is unbiased)

#256's increment matrix is ``Gamma_j = E[D_j D_j^T]`` with ``D_j = m_j - m_{j-1}``
and ``m_j = E[Y | F_j]``.  For a martingale the increments are orthogonal, so

    Gamma_j = E[m_j m_j^T] - E[m_{j-1} m_{j-1}^T].

Each ``E[m_j m_j^T]`` is a conditional expectation squared, which has the exact
U-statistic form ``E[Y(c) Y(c'_j)^T]`` where ``c`` and ``c'_j`` agree on the
bonds revealed by ``F_j`` and are independent elsewhere.  At ``p = 1/2`` that
partner is ``c'_j = (c & S_j) | (r & ~S_j)`` for a fresh random ``r``.  So

    Gamma_j = E[Y(c) Y(c'_j)^T] - E[Y(c) Y(c'_{j-1})^T],

estimated by plain paired averages.  Unlike a conditional-replica estimate of
``m_j`` itself, this has **no replica-noise bias**: the paired product already
integrates the replica out in expectation.  Its consistency check is
telescoping -- the sum over scales must equal the total covariance -- and the
L=3 run must reproduce the enumerated gate's ``Gamma_j`` to within sampling
error.  Both are asserted and tested.

## Claim boundary

Nothing here names a mechanism.  ``B_even`` is a duality-even measure tangent,
not automatically a local energy field; ``X`` is topological, but a covariance
with it does not demonstrate a defect theory.  The only question answered is
whether the two typed channels keep a distinguishable scale organisation on
tori larger than the enumerated ones.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path
from typing import Any, Sequence

try:  # pragma: no cover - import shape depends on how the script is invoked
    from scripts.square_bond_kappa3 import square_bond_pairs
    from scripts.qtangent_scale_decomposition import spatial_filtration_levels
    from scripts.torus_homology import HomologyUnionFind, _extend_basis
except ModuleNotFoundError:  # pragma: no cover
    from square_bond_kappa3 import square_bond_pairs
    from qtangent_scale_decomposition import spatial_filtration_levels
    from torus_homology import HomologyUnionFind, _extend_basis

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "qtangent-empirical" / "latest.json"
SCHEMA = "matching-one.qtangent-empirical.v1"
ISSUE = 581

#: The frozen readout, the same primary channel the exact gate carried.
READOUT = "wrap_either"

#: Sizes.  L=3 is the control: its enumerated gate numbers are known exactly,
#: so a Monte Carlo run that does not reproduce them within error is a bug.
SIZES = (3, 4, 8)

#: Number of predeclared Euclidean scales.
N_SCALES = 4

#: Configurations sampled per size.
SAMPLES = 100_000

#: Seed that fixes the whole random stream.
SEED = 581


def _primal_statistics(length: int, mask: int, pairs: Sequence[Any]) -> dict[str, int]:
    """One union-find pass over the primal open edges.

    Returns the same quantities as ``configuration_statistics`` plus the two
    wrapping flags, so a configuration is read exactly once rather than three
    times.  ``span_rank`` is the rank of the image of ``H1(A)`` in ``H1(T^2)``.
    """
    vertices = length * length
    union_find = HomologyUnionFind(vertices, (length, length))
    edges = 0
    for index, pair in enumerate(pairs):
        if (mask >> index) & 1:
            union_find.add_edge(*pair.primal)
            edges += 1
    seen: set[int] = set()
    basis: list[Any] = []
    max_rank = 0
    top_rank = 0
    for vertex in range(vertices):
        root, _, _ = union_find.find(vertex)
        if root in seen:
            continue
        seen.add(root)
        rank = union_find.component(vertex).rank
        max_rank = max(max_rank, rank)
        top_rank = max(top_rank, rank)
        for winding in union_find.basis[root]:
            _extend_basis(basis, winding)
    components = len(seen)
    return {
        "components": components,
        "open_edges": edges,
        "cycle_rank": edges - vertices + components,
        "span_rank": len(basis),
        "max_component_rank": max_rank,
        "wrap_either": int(top_rank > 0),
        "wrap_cross": int(top_rank == 2),
    }


def _dual_mask(mask: int, pairs: Sequence[Any], primal_index: dict[Any, int]) -> int:
    """The geometric dual transport, with the primal->index map precomputed."""
    transported = 0
    for bond_index, pair in enumerate(pairs):
        if not (mask >> bond_index) & 1:
            transported |= 1 << primal_index[pair.dual]
    return transported


def evaluate(length: int, mask: int, pairs: Sequence[Any],
             primal_index: dict[Any, int]) -> tuple[int, int, int]:
    """The joint vector for one configuration, as integers.

    Returns ``(wrap_either, twice_b_even, x_source)`` where ``twice_b_even`` is
    ``2 * B_even`` -- kept doubled so the whole pipeline is integer arithmetic.
    ``x_source = span_rank - 1`` is the ambient-homology source.
    """
    primal = _primal_statistics(length, mask, pairs)
    dual = _primal_statistics(length, _dual_mask(mask, pairs, primal_index), pairs)
    betti = primal["components"] + primal["cycle_rank"]
    betti_dual = dual["components"] + dual["cycle_rank"]
    twice_b_even = betti + betti_dual
    x_source = primal["span_rank"] - 1
    return primal["wrap_either"], twice_b_even, x_source


def coarse_scales(length: int, n_scales: int) -> list[tuple[int, int]]:
    """``n_scales`` nested Euclidean scales, as cumulative bond selectors.

    The base filtration is #581's radius-ordered nesting
    ``spatial_filtration_levels``, the same one the enumerated gate used.  When
    that nesting already has ``<= n_scales`` levels (L=2, L=3) it is returned
    unchanged, so the L=3 control is comparable level-by-level with the exact
    gate.  Otherwise the levels are coarsened to ``n_scales`` anchors spread
    evenly over the level index.  Each entry is ``(selector, revealed_bond_count)``
    with the last entry revealing every bond.
    """
    levels = spatial_filtration_levels(length)
    cumulative: list[tuple[int, int]] = []
    acc = 0
    for level in levels:
        for index in level:
            acc |= 1 << index
        cumulative.append((acc, acc.bit_count()))
    if len(cumulative) <= n_scales:
        return cumulative
    # Coarsen to n_scales anchors evenly spread over the level index.
    anchors = sorted({round((len(cumulative) - 1) * s / (n_scales - 1))
                      for s in range(n_scales)})
    return [cumulative[a] for a in anchors]


def estimate_gammas(length: int, scales: Sequence[tuple[int, int]],
                    samples: int, seed: int, pairs: Sequence[Any],
                    primal_index: dict[Any, int]) -> dict[str, Any]:
    """Paired Monte Carlo estimate of the per-scale increment matrix.

    Returns the raw 3x3 ``Gamma_j`` (over ``wrap_either, twice_b_even,
    x_source``), the telescoping residual, and the total covariance.
    """
    rng = random.Random(seed)
    bonds = 2 * length * length
    n = 3  # wrap_either, twice_b_even, x_source
    # cross[j][a][b] accumulates sum_s Y_a(c_s) Y_b(partner_j(c_s)).
    cross = [[[0.0 for _ in range(n)] for _ in range(n)] for _ in scales]
    # mean[a] accumulates sum_s Y_a(c_s), for the j = -1 base E[Y] E[Y]^T.
    mean = [0.0 for _ in range(n)]
    for _ in range(samples):
        mask = rng.getrandbits(bonds)
        y = evaluate(length, mask, pairs, primal_index)
        for a in range(n):
            mean[a] += y[a]
        for j, (selector, _) in enumerate(scales):
            partner_mask = (mask & selector) | (rng.getrandbits(bonds) & ~selector)
            yp = evaluate(length, partner_mask, pairs, primal_index)
            for a in range(n):
                for b in range(n):
                    cross[j][a][b] += y[a] * yp[b]
    # Base: E[Y] E[Y]^T at the empty filtration (j = -1).
    base = [[mean[a] * mean[b] / samples**2 for b in range(n)] for a in range(n)]
    gammas = []
    running = [[0.0 for _ in range(n)] for _ in range(n)]
    for j in range(len(scales)):
        gamma = [[cross[j][a][b] / samples - base[a][b] for b in range(n)]
                 for a in range(n)]
        if j > 0:
            for a in range(n):
                for b in range(n):
                    gamma[a][b] -= cross[j - 1][a][b] / samples - base[a][b]
        for a in range(n):
            for b in range(n):
                running[a][b] += gamma[a][b]
        gammas.append(gamma)
    # Total covariance by the paired convention: cross[last] - base, where the
    # last partner equals the configuration itself.
    total = [[cross[-1][a][b] / samples - base[a][b] for b in range(n)]
             for a in range(n)]
    residual = [[running[a][b] - total[a][b] for b in range(n)] for a in range(n)]
    return {
        "gammas": gammas,
        "total_covariance": total,
        "telescoping_residual": residual,
    }


def degree_one(length: int, samples: int, seed: int, pairs: Sequence[Any],
               primal_index: dict[Any, int]) -> dict[str, Any]:
    """Monte Carlo degree-one Walsh coefficients of ``B_even`` and ``X``.

    ``f({i}) = E[f chi_i]`` with ``chi_i = (-1)^{bit_i}``, so
    ``f({i}) = 1/2 (E[f | bit_i=0] - E[f | bit_i=1])``.  The exact gate's parity
    argument says ``B_even`` has all-zero degree-one coefficients and ``X`` has
    one constant value; this measures both rather than assuming them.
    """
    rng = random.Random(seed)
    bonds = 2 * length * length
    # sums conditioned on bit i: (sum_when_1, sum_when_0, count_1, count_0).
    be_s1 = [0.0] * bonds
    be_s0 = [0.0] * bonds
    x_s1 = [0.0] * bonds
    x_s0 = [0.0] * bonds
    c1 = [0] * bonds
    c0 = [0] * bonds
    for _ in range(samples):
        mask = rng.getrandbits(bonds)
        _, twice_be, x = evaluate(length, mask, pairs, primal_index)
        for i in range(bonds):
            if (mask >> i) & 1:
                be_s1[i] += twice_be
                x_s1[i] += x
                c1[i] += 1
            else:
                be_s0[i] += twice_be
                x_s0[i] += x
                c0[i] += 1
    twice_be_coeff = [0.5 * (be_s0[i] / c0[i] - be_s1[i] / c1[i])
                      for i in range(bonds)]
    x_coeff = [0.5 * (x_s0[i] / c0[i] - x_s1[i] / c1[i]) for i in range(bonds)]
    # B_even itself is half of twice_b_even, so halve its degree-one coefficient.
    be_coeff = [v / 2 for v in twice_be_coeff]
    return {
        "twice_b_even_degree_one_max_abs": max(abs(v) for v in twice_be_coeff),
        "betti_even_degree_one_max_abs": max(abs(v) for v in be_coeff),
        "betti_even_degree_one_coefficients": be_coeff,
        "x_degree_one_value": float(sum(x_coeff) / len(x_coeff)),
        "x_degree_one_spread": max(abs(v - sum(x_coeff) / len(x_coeff))
                                   for v in x_coeff),
        "x_degree_one_coefficients": x_coeff,
    }


def analyse_size(length: int, samples: int, seed: int) -> dict[str, Any]:
    pairs = square_bond_pairs(length)
    primal_index = {pair.primal: i for i, pair in enumerate(pairs)}
    scales = coarse_scales(length, N_SCALES)
    est = estimate_gammas(length, scales, samples, seed, pairs, primal_index)
    deg1 = degree_one(length, samples, seed + length, pairs, primal_index)

    # Cumulative coupling curves for the readout against B_even and against X,
    # in natural units: Gamma[wrap, B_even] = Gamma[wrap, twice_b_even]/2.
    gammas = est["gammas"]
    cumulative_betti = []
    cumulative_x = []
    run_b = 0.0
    run_x = 0.0
    scale_rows = []
    for j, gamma in enumerate(gammas):
        run_b += gamma[0][1] / 2
        run_x += gamma[0][2]
        cumulative_betti.append(run_b)
        cumulative_x.append(run_x)
        scale_rows.append({
            "level": j,
            "revealed_bonds": scales[j][1],
            "gamma_wrap_b_even": gamma[0][1] / 2,
            "gamma_wrap_x": gamma[0][2],
            "gamma_raw": gamma,
        })

    total_betti = est["total_covariance"][0][1] / 2
    total_x = est["total_covariance"][0][2]
    # Sign coherence: is the X profile single-signed and the B_even profile
    # alternating, as the enumerated gate found?
    x_inc = [r["gamma_wrap_x"] for r in scale_rows]
    b_inc = [r["gamma_wrap_b_even"] for r in scale_rows]
    x_single_signed = all(v >= 0 for v in x_inc) or all(v <= 0 for v in x_inc)
    b_alternates = any(a * b < 0 for a, b in zip(b_inc, b_inc[1:]))

    return {
        "length": length,
        "bonds": 2 * length * length,
        "samples": samples,
        "n_scales": N_SCALES,
        "readout": READOUT,
        "scale_rows": scale_rows,
        "cumulative_coupling_betti_even": cumulative_betti,
        "cumulative_coupling_x": cumulative_x,
        "total_cov_wrap_b_even": total_betti,
        "total_cov_wrap_x": total_x,
        "x_profile_single_signed": x_single_signed,
        "betti_profile_alternates": b_alternates,
        "telescoping_residual": est["telescoping_residual"],
        "total_covariance_raw": est["total_covariance"],
        "degree_one": deg1,
    }


def assemble(sizes: Sequence[int] = SIZES, samples: int = SAMPLES,
             seed: int = SEED) -> dict[str, Any]:
    per_size = [analyse_size(length, samples, seed) for length in sizes]
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "method": {
            "estimator": "unbiased paired U-statistic E[Y(c)Y(c'_j)] - E[Y(c)Y(c'_{j-1})]",
            "samples_per_size": samples,
            "n_scales": N_SCALES,
            "seed": seed,
            "readout": READOUT,
            "coordinates": ["wrap_either", "twice_b_even", "x_source"],
            "partner_rule": "partner = (c & selector_j) | (fresh & ~selector_j)",
        },
        "sizes": per_size,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", nargs="+", type=int, default=list(SIZES))
    parser.add_argument("--samples", type=int, default=SAMPLES)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    payload = assemble(tuple(args.sizes), args.samples, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {args.output}")
    for row in payload["sizes"]:
        print(f"  L={row['length']:2d}  total Cov(wrap,B_even)={row['total_cov_wrap_b_even']:+.5f}  "
              f"Cov(wrap,X)={row['total_cov_wrap_x']:+.5f}  "
              f"x_single_signed={row['x_profile_single_signed']}  "
              f"b_alternates={row['betti_profile_alternates']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
