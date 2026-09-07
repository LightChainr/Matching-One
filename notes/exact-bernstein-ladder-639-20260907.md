# Exact Bernstein ladder to L=5 (and honest L=6 attempt) — #639, 2026-09-07

Ticket: #639. Artifact: `results/exact-bernstein-ladder-639/latest.json`.
Kernel: `scripts/exact_matching_ladder_bruteforce.cpp` (two kernels, K1 and K2).
Checks: `scripts/verify_power_sums.py`, `scripts/exact_roots_physical.py`.
Branch: `exact-bernstein-ladder-639` (from `claude/matching-one-workspace-pwr5pv`).

## What was run

Two independent brute-force kernels over the full 2^N configuration space,
pure integer arithmetic end to end (g++ -O2, ARM 16-vCPU Huawei fleet):

- **K1** — per-configuration rebuild of the displacement-aware union-find,
  a faithful C++ port of `scripts/exact_matching_polynomial.py` +
  `scripts/matched_torus_reference.py`. Guard: refuses N > 40 (the Python
  reference still refuses N > 26).
- **K2** — DFS over site colors with a fully journaled incremental
  union-find (undo on backtrack). No transfer matrix, no symmetry
  reduction, no symmetry assumptions. Adjacency keeps every directed edge
  (no dedup), matching the reference's behavior on small-period tori.

**Regression before any new rung was reported:** both kernels reproduced all
five committed rungs (axis L=2/3/4, diamond L=2/3) bit-for-bit against
`results/exact_small_matching_polynomials.md`.

## New rungs

| rung | N | configs | counts agree | K2 wall | K1 wall |
|---|---:|---:|---|---:|---:|
| axis L=5 | 25 | 2^25 | K1 = K2, 26 coefficients | 0.26 s | 2.06 s |
| diamond L=4 | 32 | 2^32 | K1 = K2, 33 coefficients | 36.2 s | 961.7 s (7 thr, fleet) |
| axis L=6 | 36 | 2^36 | K1 = K2, 37 coefficients | 526.2 s | 5416.5 s (14 thr, fleet) |

diamond L=4 was additionally confirmed by a **third independent
environment**: a locally built K1 on the author's Mac (ARM64, 8 threads,
2753.7 s) reproducing the identical 33 integers. Three agreeing
computations of the same 2^32-space count, from two kernels and two
machines.

The three count vectors (see the artifact) pass every structural check in
`scripts/verify_power_sums.py`: `a_0 = M(0) = -1` and `a_N = M(1) = +1`
exactly, integer coefficients, degree exactly N — on the three new rungs
and, as regression, on the five committed ones.

## Roots and factorization (exact)

`scripts/exact_roots_physical.py`: Bernstein → power basis over
`fractions.Fraction`, `sympy.factor_list` over ZZ, Sturm `count_roots(0,1)`
for uniqueness, exact rational isolation intervals, then 170 rounds of
exact `Fraction` bisection. The exact artifact is the bracket and the
factorization; the decimal is labelled diagnostic. Regression: the script
reproduces the committed axis L=2 root
`0.54119610014619698439972320536638942006107206337802` from the committed
counts.

- **axis L=5**: power-basis poly irreducible over ZZ; exactly one root in
  (0,1); diagnostic 0.59198825651833384461096868021192887904787477719722.
- **diamond L=4**: irreducible over ZZ; one root in (0,1); diagnostic
  0.59318323467361685147129923181974517395503790167638.
- **axis L=6**: irreducible over ZZ; one root in (0,1); diagnostic
  0.59239507081770423769385580764250543411218819923508.

The orientation asymmetry recorded for smaller rungs persists: axis roots
still climb toward the threshold from below, diamond roots still sit above
it, and at L=5/6 the axis estimator has not crossed the diamond estimator.

## diamond L=5 is out of independent reach — reported, not worked around

N=50 means 2^50 configurations. Measured K2 rate on one 16-vCPU box:
~1.19e8 configs/s. One machine would need ~2638 hours. Splitting across the
ten-box fleet by precomputed prefix colors (the `--split` decomposition) is
bounded below by the same total CPU, still far out of reach for this rung.

The known symmetry reduction (torus automorphism group, order 288) would
bring this to roughly 9 machine-hours, but that harness is not built here,
and the ticket's hard constraint — independent, non-circular verification —
means a symmetry-reduced K2 would itself need a brute-force cross-check on
at least one smaller rung before it could be trusted, plus a second
independent symmetry-reduced kernel for the actual L=5 number. None of that
exists yet, so **no diamond L=5 number is claimed**. Building the symmetry
harness with its own K1/K2 double-cover is the obvious next ticket.

## What K1 confirmed

All three new rungs are double-kernel confirmed:

- axis L=5: K1 = K2 bit-for-bit.
- diamond L=4: K1 = K2 bit-for-bit, plus a third independent environment
  (local Mac K1 build) agreeing again.
- axis L=6: K1 = K2 bit-for-bit (2^36, 5416.5 s at 14 threads).

## What we could not do

- diamond L=5: no number, reason above. This is the only rung requested by
  the ticket that remains uncomputed, and the honest limit is the reported
  result.

## What finishing would need

- The symmetry-reduction harness (order-288 quotient + its own independent
  double kernel + a brute-force regression rung) to reach diamond L=5.
