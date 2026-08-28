# C00: general 2x2 integer-period torus homology engine

Queue item C00 (issue 22 / coordination issue 23). Started from PR #21
homology on `agent/c00-general-homology` at
`fad82a0fc044c377993c9a0d59aabd772b92fa7c`.  The existing rank-0/1/2
classifier was extended to an arbitrary integer period matrix.  A second
topology implementation was not added.  C01 production data was not started.

Wall time: 3.30 s (18:29:59–18:30:03 Asia/Shanghai, 2026-08-28).
Single-thread CPU, no GPU.

## What changed

`scripts/torus_homology.py` now accepts a 2x2 integer period matrix `P`
with nonzero determinant.  Closed cover displacements convert to generator
windings by the exact identity `w = adj(P) d / det(P)` (integer adjugate,
integer division, divisibility check).  Scalar `(px, py)` input remains
the diagonal special case, so axis and diamond embeddings are unchanged.

`scripts/matched_torus_reference.py` records `period_matrix` on every
geometry and adds `integer_period_geometry` / `gaussian_geometry` /
`diamond_xy_geometry` on the same homology path.

Conventions: `notes/integer-period-torus-conventions.md` (copied here).

## Stop-gate tests

| test | result | notes |
|---|---|---|
| exhaustive axis L=3 | PASS | counts equal PR #21 `{rank0:259, rank1:162, rank2:91, d0:175, d1:175}` |
| exhaustive diamond L=2 | PASS | counts equal PR #21 `{rank0:143, rank1:68, rank2:45, d0:81, d1:81}` |
| axis L=3 via general `P=diag(3,3)` | PASS | same PR #21 counts |
| diamond L=2 in xy coordinates | PASS | same PR #21 counts |
| Gaussian `(2,1)` N=5 exhaustive | PASS | `{rank0:16, rank1:10, rank2:6, d0:11, d1:11}` |
| random SL(2,Z)/det +/-1 basis change, axis L=3 | PASS | seed 20260828, 8 matrices |
| random basis change, diamond L=2 | PASS | seed 20260828, 8 matrices |
| random basis change, Gaussian `(2,1)` | PASS | seed 20260828, 12 matrices |
| full unittest suite | PASS | 25 tests, 2.08 s |

Boolean `either` wrapping still matches `cluster_stats` on every enumerated
configuration.  Gaussian cyclic labels agree with
`scripts/gaussian_circulant_geometry.py`.

## Preserved negative / design result

Generator-relative flags `direction_0`, `direction_1` and `both` are **not**
invariant under unimodular shears.  A rank-1 wrap along the first generator
becomes a spiral after a mixing shear.  Rank, `either` and `cross` stay
invariant.  This is required generator-relative behaviour, recorded rather
than “fixed”.

## Stop gate

All exhaustive and basis-invariance tests passed.  C01 was not started.

Machine-readable counts: `regression.json`.
