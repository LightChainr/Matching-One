# Completion-hazard production: fixed estimands and interpretation

This scorer uses the newly generated independent square block only. It does
not pool older random data or count derived readouts as independent evidence.
The initial design is L=512, N=262144, 14 equal batches of 10000 filtrations,
with fixed insertion cutoffs a=154646, b=155385, c=156120. The scorer reads
these parameters from the generating contract rather than embedding them in
the estimator. The clock is uniform site-permutation **insertion count**,
not an iid-label time or a fitted near-critical clock.

## Inputs and invocation

```sh
python3 analysis/completion-hazard-production-20260929/analyze.py \
  --data-dir analysis/completion-hazard-production-20260929/data \
  --output-dir analysis/completion-hazard-production-20260929/results
```

`contract.json` must be in the parent of `--data-dir`. Named integer fields
`a,b,c,L` are required, either at the top level or nested inside objects;
repeated occurrences must agree. N is read if present, otherwise inferred
as L squared. The square geometry and `0<=a<b<c<=N` are checked. Matching
named fields in `data/run.json` are cross-checked, and absent metadata names
are reported, not presented as verified matches.

`data/run.json["batches"]` is a nonempty list of objects with `file`,
`samples`, `seed`, `sha256`. File names are relative to data-dir (or absolute
inside it). Distinct files/seeds/content hashes are required; duplicated
blocks are not independent replicates. Each listed CSV or `.csv.gz` is read
**once**; the same raw bytes are hashed, decompressed and parsed. Its
expected number of samples and every required row field are checked.

**The fixed production block must be complete before any scoring.** Require
`run.status == "completed"`, `len(run.batches) == contract.batches`, and
each `batch.samples == contract.samples_per_batch`. The latter two contract
fields are required top-level integers. For the current design this means
all 14 batches, each with 10000 rows; the CSV row counts are also checked
against those declarations. An in-progress run or its subset of finished
batches is rejected, not described as a complete analysis. All contracted
batches are included; there is no data-dependent batch selection.

Required CSV columns are

    J1,J2,dx_b,dy_b,nu_b,tau,nu_tau.

For rank one at b, D=(dx_b,dy_b) is a nonzero primitive unoriented direction.
Only its sign is canonicalized: the first nonzero coordinate is positive.
Nonprimitive directions are errors, not silently rounded or merged.
`nu_b=-1` outside rank one. At tau the analogous rule holds for `nu_tau`.
Require `b<=tau<c`, and `nu_tau=nu_b` if tau=b. These are semantic checks,
not a new validation of the engine's topological count. Hashes and range
checks cannot establish that tau is independent of the filtration or that
the seeds generate independent blocks; those remain generation assumptions.

## Fixed initial cohorts

The two initial cohorts partition all rank-one risk at b:

    E: J1<=a, J2>b;
    L: a<J1<=b, J2>b.

Here L used as a cohort label means **late**, not the linear size. These
cohorts are fixed using the history by b. They are not redefined by rank
at tau or by survival to c. Rows outside rank one at b are excluded from
these cohort estimands, with their count reported.

For each member of either cohort define

    Z = 1{J2>c},
    I = (c-b) * 1{J2>tau} * nu_tau/(N-tau).

When J2<=tau, I is zero; the sentinel -1 is not multiplied into it.
Otherwise the initial cohort condition ensures rank one at tau. I is
**not** truncated to [0,1]; its individual value can exceed one.

## Why the random-time compensator is unbiased

Let G_k be the revealed permutation history. A rank-one set A_k completes
at the next insertion with conditional probability nu_2(A_k)/(N-k).
For either initial cohort g in G_b and every k>=b,

    E[1_g 1{J2>k} nu_2(A_k)/(N-k)] = P(g,J2=k+1).

Independence and uniformity of tau on `{b,...,c-1}` therefore give

    E[I | g]
      = sum_(k=b)^(c-1) E[1{J2>k}nu_2(A_k)/(N-k) | g]
      = P(b<J2<=c | g),
    E[Z+I | g] = 1.                                  (1)

Conditioning on the entire sampled filtration instead of g gives
`E_tau[I | filtration]` equal to its sum of conditional one-step hazards,
not generally to its realized completion indicator. There is no pointwise
identity `I=1-Z`. Completion changes the risk set along the path, which is
why retaining `1{J2>tau}` is essential. Fresh entrants after b are not part
of this compensator cohort.

## Primary scalar outputs and signs

The JSON uses these ordered primary scalars; all enter the joint aligned
batch covariance:

| Key | Definition |
|---|---|
| `delta_survival` | mean Z in E minus mean Z in L |
| `delta_integrated_hazard` | minus mean I in E plus mean I in L |
| `closure_residual` | delta_survival minus delta_integrated_hazard |
| `delta_nu_instantaneous_scaled` | (c-b)/(N-b) times (mean nu_b in E minus mean nu_b in L) |
| `delta_overlap_D` | overlap-weighted survival difference within exact D |
| `delta_overlap_D_nu` | overlap-weighted survival difference within exact (D, integer nu_b) |
| `closure_early` | mean Z in E plus mean I in E minus one |
| `closure_late` | mean Z in L plus mean I in L minus one |

Equation (1) makes the expected pooled survival and integrated-hazard
differences equal. The groupwise residuals additionally expose common
closure bias that could cancel in their difference. The per-cohort means
of Z, I and nu_b are also reported as descriptive components.

The instantaneous scaled nu difference has the **hazard sign**: a larger
early nu_b predicts lower *next-step* early survival. Multiplication by
c-b is only an explicitly labelled window scaling. It is not integration
over the evolving completion intensity and does not assume geometry or
hazard stays constant during the interval.

Closing the random-time identity on the new block supports the implemented
geometry-to-exit map and compensator. It does **not** demonstrate that nu_b
is a sufficient recursive state: I observes geometry at a later random
time, whereas a model based on the state at b must predict its evolution.

## Two fixed common-support conditional comparisons

For either fixed stratification S=D or S=(D,nu_b), let a cell contain early
count e_s and late count l_s, with survival counts z_Es and z_Ls. Only
cells with both e_s>0 and l_s>0 contribute. Define

    w_s = e_s*l_s/(e_s+l_s),
    Delta_S = sum_s w_s*(z_Es/e_s-z_Ls/l_s) / sum_s w_s.       (2)

D is the exact primitive unoriented integer direction; nu_b is its exact
integer value. There are no bins, regressions, fitted cutoffs, new features,
or post-hoc subgroup searches. For each stratification report the number
of shared and one-sided cells, supported counts, and excluded fractions
in E, L and the combined initial risk set. Empty or entirely non-overlapping
groups are `not_scoreable`, not zero. The overlap-weight sum is a weight
sum, not an asserted effective sample size. Exact cell counts/survivors
are retained in JSON; every direction has descriptive cohort means.

These overlap-weighted targets differ from the pooled risk difference and
from each other. Therefore their reduction relative to the pooled estimate
is **not** a mediated fraction, explained fraction, or causal adjustment
effect. No such ratio is reported.

If the state (time b, rank one, D,nu_b) suffices to predict future survival,
then early and late entry must have equal survival probability within
every such population cell on common support. A nonzero population
Delta_(D,nu) is a genuine violation of this necessary prediction. A zero
weighted difference can hide opposite cell differences; it establishes
neither cellwise equality nor a recursively Markov compressed process.
Sampling variation must be distinguished from a population violation.
The D-only statistic addresses the analogous weaker conditioning claim.

## Aligned whole-batch delete-one uncertainty

The point estimate pools all rows through per-batch sufficient statistics.
For each of B batches, remove that entire batch simultaneously for all
cohorts, all directions and all scalar outputs, then recompute (2), including
which cells still share support and their weights. No direction is deleted
independently and no original full-data stratum weights are held fixed.

For the vector theta_(-j) and its delete-one mean theta_bar,

    Cov_JK = (B-1)/B * sum_j
               (theta_(-j)-theta_bar)(theta_(-j)-theta_bar)^T.

The JSON contains the complete ordered covariance matrix, all B aligned
replicate vectors and marginal SEs. The estimates shown are pooled plug-in
estimates, not jackknife bias-corrected estimates. Exact identities such
as `closure_residual=closure_early-closure_late` make a singular covariance
expected; the scorer neither inverts it nor invents independent tests.

This formula assumes the planned equal-size independent batches. The CLI
rejects partial or unequal-size production blocks using the fixed contract
above. Internally the covariance routine additionally returns `not_scoreable`
for fewer than two batches or unequal batch sizes. If one scalar loses its necessary groups
or common support in any delete-one replicate, its SE and covariance row/
column are null with an explicit reason. Other fully scoreable scalars still
use all the same B deletions; no pairwise subset covariance is substituted.

## Files produced and boundaries

The command writes `result.json` and `RESULT.md` under --output-dir. Structural
input errors (hash mismatch, invalid rows, missing files, inconsistent
parameters, duplicated batch identities) terminate the run rather than
silently omitting data. A failed scientific comparison due to absent groups
instead remains in the report as `not_scoreable`.

This is one fixed finite-lag comparison from one new random block. Survival,
compensator, directions and completion counts share its randomness. Old
archives are not independently pooled. The scorer makes no statement about
causal intervention, original-U operator identity, iid-label-clock closure,
or a limiting non-Markov theorem. The finite model distinction here is
whether the specified b-state leaves predictive entry-history dependence,
not whether adding ever more descriptive coordinates yields a smaller
retrospective residual.
