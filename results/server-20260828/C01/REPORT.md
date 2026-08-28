# C01: same-N Gaussian orientation discovery

Queue item C01 (issue 22 / coordination issue 23). Started from C00
`agent/c01-gaussian-orientation` at
`08c4079798cc37913c6a7de3c97bc8d1c6bda16a`.  The C00 2x2 homology engine is
the single DSU; `src/pell_matching_mc.cpp` now uses
`src/homology_union_find.hpp` (`either` = PR #21 Boolean wrap).  A second
topology implementation was not added.  `scripts/torus_homology.py` math was
not rewritten.

Wall time: 25.49 s (18:46:24–18:46:49 Asia/Shanghai, 2026-08-28).
OpenMP 8 threads, no GPU.

## Frozen choices (pilot, not reused)

p_ref = 0.592746050790 (discovery coordinate, not a pc claim); complement =
0.407253949210. Bond control at exact pc = 1/2.

Graph-only structural multipliers, then 200000 independent site-pilot
replicas per discovery pair (Philox4x32-10 seed 20260828, replicas
[0, 200000)). Evaluation uses the same seed with replicas
[1000000, 3000000). Bond uses seed 20260829 and the same eval replica range.

| N | t_struct (graph-only) | frozen t | reason |
|---|---|---|---|
| 65 | 4 | 1 | min replica var of ΔS_either; t=4 not better |
| 85 | 3 | 1 | min replica var of ΔS_either; t=3 not better |
| 145 | 8 | 1 | canonical t=1 within 5% of min var |

Primary effect channel frozen from the N=65 matching-even pilot ranking:
`S_both`.  All five homology channels are still reported.  Default coupling
is same_U_j (t=1).

## Bond control: PASS

Square-bond percolation at p=1/2 on the same Gaussian engine.  Matching-graph
bond wrapping is essentially always true (degree-8 circulant at p=1/2), so
the angular signal lives in primal wrapping.

| run | channel | Δ primal | batch SE | z | N Δ / Δcos4 |
|---|---|---|---|---|---|
| bond N=65 | both | -0.006048 | 0.000290 | -20.82 | -0.2884 ± 0.0138 |
| bond N=85 | both | -0.0048705 | 0.000428 | -11.38 | -0.2596 ± 0.0228 |
| bond N=65 | either | +0.0022405 | 0.000401 | +5.58 | +0.1068 ± 0.0191 |
| bond N=85 | either | +0.0022435 | 0.000317 | +7.09 | +0.1196 ± 0.0169 |

The frozen channel `both` is the largest bond signal.  Scaled amplitudes at
N=65 and N=85 agree at the ~1–2 SE level, with the sign of `both` opposite
to Δcos4 (so a negative spin-4-like bond amplitude).  Geometry, period-matrix
homology, cyclic coupling, and the angular estimator are functioning.

## Site discovery: matching-even NOT resolved; matching-odd NOT resolved

Evaluation: 2000000 independent site replicas at p_ref, t=1, N=65 (8,1) vs
(7,4) and N=85 (9,2) vs (7,6).  Uncertainties are batch SE from 20 equal
batches.

Matching-even ΔS = S(θ1)−S(θ2), S=(R_primal+R_matching)/2 on the frozen
channel `both`.  Matching-odd D_N uses M_either = primal_either −
matching_either (identical for all five channels on these quotients, as in
the C00 tiny-torus identity).

| run | ΔS_both | SE | z_S | D_N (M_either) | SE | z_M | A4_M = N^{13/8} D_N / Δcos4 |
|---|---|---|---|---|---|---|---|
| site N=65 | -0.00021975 | 0.000276 | -0.80 | +0.0007695 | 0.000460 | +1.67 | +0.498 ± 0.298 |
| site N=85 | -0.00055775 | 0.000342 | -1.63 | +0.0002145 | 0.000437 | +0.49 | +0.184 ± 0.374 |

No site homology channel reaches |z|≥3 in ΔS, Δ primal, or D_N.  The largest
site |z| is N=85 ΔS_either z=+1.98 (and Δ matching_cross z=-2.47), still
below the predeclared resolution threshold.

CRN correlation of primal_either across the two orientations is 0.333 (N=65)
and 0.304 (N=85).  Relabel t=t_struct did not reduce replica variance of ΔS
relative to t=1 on the 200k pilots.

## Power / sensitivity (negative site result)

A naive O(1)×L^{-2} matching-even amplitude would give
|ΔS| ~ Δcos4 / N ≈ 0.021 at N=65.  The measured |ΔS_both| = 2.2e-4 is smaller
by two orders of magnitude; that O(1) L^{-2} hypothesis is excluded by tens
of sigma.  3σ upper bounds at 2e6 samples:

- N=65: |ΔS_both| < 8.3e-4
- N=85: |ΔS_both| < 1.03e-3

Matching-odd D_N is likewise consistent with zero (z=1.67 and 0.49).  The
A4_M point estimates are positive (same sign as Δcos4) but the batch SEs are
larger than the estimates.  Detecting |D_N|=0.001 at z=3 would need roughly
the present sample count; the observed |D_N| is smaller.

Increasing N makes an L^{-2} (or faster) signal smaller, so N=145 / 205 /
425 / 1105 evaluation was not started.

## Stop gate

- Bond control angular signal: **PASS**
- Site matching-even (larger sector) at N=65 or N=85: **FAIL to resolve**
- Site matching-odd D_N: **FAIL to resolve**

Decision: **stop**.  Do not increase N.  Do not start C02 production analysis.
This is not a geometry-debug stop: the exact bond control on the same engine
resolves a reproducible orientation harmonic, so the site null is a preserved
negative result with a power statement.

## RNG

Philox4x32-10 (Random123 official KATs in `--self-test`).
key=(seed_lo, seed_hi); ctr=(index, replica_lo, replica_hi, stream);
stream 0 = site occupation by cyclic vertex j; stream 1 = bond occupation
by packed (src, dx, dy). Uniforms are the top 53 bits of the first two
32-bit words, in [0,1).

Site seed 20260828, bond seed 20260829.
Pilot replicas [0, 200000); evaluation replicas [1000000, 3000000).
Coupling: U_j^{(1)}=U_j, U_j^{(2)}=U_{t j mod N} with frozen t=1.

## Tests

- Gaussian (2,1) exhaustive: rank0=16, rank1=10, rank2=6, d0=11, d1=11 (C00)
- Gaussian (3,2) exhaustive: rank0=4629, rank1=2340, rank2=1223, d0=2471, d1=2471
- PR #21 axis L=2,3 and diamond L=2 matching polynomials still pass after the DSU swap
- Unittest (homology + Pell MC + Gaussian MC): PASS (21 tests)

## Negative / preserved results

All channels, both signs, and the site null are retained in
`long_form_channel_means.csv` and `covariance_matrices.json`.  Pilot t
shortlists are in `coupling_pilot.json`.  Confirmation sizes 205/425/1105
and C02 sector fits were not run.

C00 result files under `results/server-20260828/C00/` were not modified.
