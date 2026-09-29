# Independent L512 completion-hazard experiment

2026-09-29. Draft #838; continuation of the
[finite completion-geometry mechanism](birth-completion-geometry-20260929.md).

**Outcome: completed, 140,000 new filtrations.** The new block resolves
differences in completion geometry more clearly than the raw finite-lag
survival contrast. After exact direction conditioning, or direction plus
integer completion count, the predeclared weighted survival differences
are not resolved. This shifts the next candidate toward a direction-mixture
model; it does not establish direction sufficiency or license adding an
unlimited sequence of descriptors to the same block.

| Fixed readout | Estimate, percentage-point scale | Aligned batch SE, same scale |
|---|---:|---:|
| Early-minus-late survival | +0.686 | 0.444 |
| Within exact direction, overlap-weighted survival | -0.415 | 0.451 |
| Within exact (direction, nu_b), overlap-weighted survival | +0.122 | 0.476 |
| Integrated-hazard prediction of pooled survival difference | +1.550 | 0.351 |
| Actual minus integrated prediction | -0.865 | 0.639 |
| Instantaneous completion contrast, window-scaled, hazard sign | -1.504 | 0.330 |

The last row is dimensionless with a percentage-point display convention,
**not** an observed finite-window exit probability. The two conditional
comparisons have different target weights from the pooled row. All outputs
share this block and its full covariance. The uncertainty is a 14-batch
delete-one estimate, not an exact finite-sample confidence guarantee.

There are 24,218 early and 30,505 late rank-one histories. Their mean
completion counts are 99.402 and 101.587. The early cohort is 92.295% axial,
the late cohort 86.094%; axial directions also have greater subsequent
survival than diagonals in the descriptive table. Composition is therefore
a concrete candidate explanation for pooled memory, not a measured causal
or mediated fraction. Exact-D comparison retains 99.973% of the risk set;
exact-(D,nu_b) comparison retains 98.898%, across 1,202 shared cells.

The random-time prediction is an unbiased estimator of the **same
population pooled contrast**, but not the same sample statistic as raw
survival. Its stronger estimate-to-SE ratio is not another independent
replication. Actual-minus-predicted closure is about 1.35 batch SE from
zero; groupwise residuals are -0.009411 +/- 0.006315 (early) and
-0.000766 +/- 0.004129 (late). No resolved inconsistency with the exact
hazard map appears here; this is not proof of correctness or state closure.

The old 70k square-L512 contrast was +1.422 +/- .511 percentage points,
with cutoffs re-estimated on its own batch deletions. The new raw estimate
is smaller and retains only about 1.5 estimate/SE. We do not hide that
weakening, pool the blocks, or redefine the primary outcome around the
stronger ancillary hazard statistic. The old exact L3 results are unchanged.

Full [numerical readout](../analysis/completion-hazard-production-20260929/results/RESULT.md),
[covariance and cells](../analysis/completion-hazard-production-20260929/results/result.json),
and [raw batch manifest](../analysis/completion-hazard-production-20260929/data/run.json)
are retained.

## Question and before-data choices

The L3 calculation supplied an exact microscopic hazard and, at that size,
a binary rank-one phase realization. The existing L512 archive showed a
small birth-history association but stored no geometry. The new block asks
whether current direction, or current direction plus the **exact integer**
completion count, removes a specified finite-lag history contrast at L512.
It is new acquisition, not a reconstruction of unavailable old snapshots.

The [contract](../analysis/completion-hazard-production-20260929/contract.json),
engine and runner were committed and pushed at
`9520c9bdf702a62435f1054e64184b1df4fc0ab3` **before** new sampling.
The contract specifies:

- Square NN torus, L=512, N=262144, count clock only.
- Fixed `(a,b,c)=(154646,155385,156120)`, taken from the previous independent
  block's birth-mixture quartiles. No re-estimation on this new block.
- Fourteen independent batches of 10,000 permutations. Separate namespace
  from all earlier runs; 32 benchmark permutations excluded.
- Early cohort `J1<=a,J2>b`; late cohort `a<J1<=b,J2>b`.
- Current direction D and completion count nu_b at b. An independent
  uniform integer audit time tau in `[b,c-1]` supplies nu_tau when rank one.
- Three primary contrasts: pooled survival difference, direction-stratified
  difference, and exact `(D,nu_b)`-stratified difference. Common-support
  strata use the predeclared overlap weight `n_E*n_L/(n_E+n_L)`.
- Full aligned delete-one-batch covariance; no cutoff scan, feature search,
  significance stopping, old/new pooling or added size/geometry.

Each predictive-sufficiency null implies zero for its corresponding
contrast. Zero is not a proof of full Markovness. The stratified targets
use different weights from the pooled target; their ratio is **not** a
fraction of memory explained or a causal mediation effect.

## Efficient exact completion geometry

The new [engine](../analysis/completion-hazard-production-20260929/engine.cpp)
retains the existing lifted union-find convention. It scans vacant sites
only at b and tau, without changing the permutation or topology.

To see why a trial site's completion count can be evaluated without a
full rebuild, contract each old occupied component. A newly added vertex
and its incident edges form a star with possibly repeated component ends.
Every new independent cycle is generated by a pair of neighbours in the
**same** old component. For each occupied neighbour, its union-find root
and lifted potential specify the new vertex's implied displacement from
that root. The difference between two such displacements is a torus-period
cycle. Existing cycles span the current rank-one line; rank two appears
iff one new cycle is not parallel to it. Distinct old components linked
once each do not create an additional cycle.

This gives an `O(N d^2 alpha(N))` snapshot scan for fixed degree d, using
the existing union-find arrays. The direction is primitive and unoriented.
Path compression changes representation, not geometry or future sampling.

One narrow preproduction control used the 512 L3 subsets of each lattice:
648 square and 810 triangular explicit copied-state insertion probes
matched the fast count. All 162 square and 180 triangular rank-one
direction/count rows also matched the previously saved independent
lifted-DFS census. This checks the new measurement convention; it is not
another 9! enumeration or a test of the scientific null.

## A random-time integral, not a midpoint proxy

For a trajectory in either cohort at b, define

    I=(c-b) 1{J2>tau} nu_tau/(N-tau).

The audit time is uniform and independent of the permutation. Conditioning
on a cohort g and averaging over tau, then using the exact insertion hazard,

    E[I | g]
      = sum_(k=b)^(c-1) E[1{J2>k}nu_2(A_k)/(N-k) | g]
      = P(b<J2<=c | g).

Thus `E[I|late]-E[I|early]` predicts the early-minus-late survival difference
without assuming that the midpoint hazard persists over the whole window.
An individual I can exceed one; it is an unbiased integral estimator, not
a pathwise probability. The actual/predicted closure residual uses their
**paired** covariance. Agreement checks the implemented microscopic map,
not the claim that nu_b is recursively sufficient: only the conditional
future-history contrasts address that latter question.

The scaled midpoint difference
`(c-b)/(N-b)*(mean(nu_b|early)-mean(nu_b|late))` is also retained, explicitly
as an instantaneous readout, not an alternative finite-lag prediction.

## Execution provenance

Only DevEnvC_TV2N0X, UUID `4a8d1d443419434889e49148ed0a7ba6`, is used.
The target was Ready, then started for this job; one authorized SSH-key
rotation after Running succeeded and the key mode was 600. No other
environment was started, reset or inspected for jobs.

Observed container limits: aarch64, nproc=16, CPU quota 14.5 cores, memory
25 GiB. At entry there were no other research processes. The task uses
14 workers in the new, previously absent directory
`/workspace/Matching-One-TV2N0X/completion-hazard-20260929`.
After the container restart, g++ was absent; gcc-c++ 10.3.1 and its package
dependencies were installed on this target. Python is 3.9.9; the runner
uses only its standard library. No SSH key or credentials are transferred.

The three uploaded source files match their committed SHA256 hashes.
The 32-sample excluded benchmark took .661 seconds of engine time and
predicted about 209 seconds under ideal 14-worker scaling. The runner
records actual command, source/binary hashes, seeds, platform, cgroup,
batch completion and timings. Scientific results and terminal status are
reported only after the fixed block is complete: all 14 batches succeeded,
with runner wall time **325.512 seconds** (individual engines 308.9--323.2 s).
All 140,000 rows were scored once locally using the existing Python 3.11.15
research environment. Every compressed batch matched its recorded SHA256.
No production was extended, discarded or repeated after reading its result.

The manifest SHA256 is
`7c63af9d9c7dcf5eff4b0f14644fe83560d73190532d8528f3e0cf0c5f0af2df`;
the result JSON SHA256 is
`35cfa9ae64e97b4ef81e4fb339e4c37fa1f6615c8968510e32711c79c3f88abf`.
After download and scoring, the target had no remaining research process;
its task tunnel was stopped and power-off requested. No remote data were
deleted. Terminal Ready status is recorded in the execution receipt.

## Next scientific choice, not a new production queue

The most economical next hypothesis is a **persistent-direction mixture**
of rank-one survival kernels. In an ideal marked-Markov model each sector
would satisfy `H_d(a,c)H_d(b,b)=H_d(a,b)H_d(b,c)`, while their sum need not.
The present finite-lag conditional contrasts are compatible with this
candidate; they neither prove sectorwise equality nor determine all times.
Opposite cell effects or changing hazard differences can cancel at a
single triple. The exact L3 square label-time counterexample remains an
explicit warning against promoting count-clock compatibility to label time.

Next theory should derive a directional exit kernel or a second specific
prediction that could contradict that mixture, under a declared clock.
If additional data are later justified, reserve that prediction for a new
block; do not turn this completed experiment into an adaptive age/feature
scan. No new larger size or extra sampling is automatically dispatched.
The original-U forward-map question remains a separate unresolved lane.

## Reproduce the analysis without new acquisition

```bash
python3 analysis/completion-hazard-production-20260929/analyze.py \
  --data-dir analysis/completion-hazard-production-20260929/data \
  --output-dir analysis/completion-hazard-production-20260929/results
```

The [estimand note](../analysis/completion-hazard-production-20260929/ESTIMANDS.md)
gives the exact signs, conditional targets and covariance formula. The
scorer requires the complete fixed block. Its source was developed in
parallel with acquisition from the precommitted contract, before scoring
the data; the acquisition source and contract, not this later scorer,
were in the pre-run commit. No full-repository or repeated test suite ran.
