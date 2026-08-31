# P439: real same-stream crosswalk has an unresolved unmarked loading

## Result first

The existing four-generation archive does **not reject a common ray**, but it
also does not resolve a nonzero unmarked `M` loading on that ray. Thus it does
not yet connect the marked-current radial law to the original Matching-One
observable. This is an additional outcome missing from #439's A–D matrix:
**the common-ray model can survive because one observer is unresolved**.

No new samples were generated. This is a real-archive score, not another
synthetic control. It reuses the already reconstructed canonical observables
and paired-orientation delete-one vectors in immutable source `6123955`.

| N | angular M, estimate ± SE | angular K_A, estimate ± SE | correlation |
|---|---|---|---|
| 85 | -0.00014931 ± 0.00100547 | -0.01605092 ± 0.00497150 | -0.7261 |
| 170 | 0.00021928 ± 0.00016207 | -0.01111149 ± 0.00097563 | -0.6914 |
| 340 | 0.00005803 ± 0.00013684 | -0.00485726 ± 0.00124889 | -0.7122 |
| 680 | -0.00003976 ± 0.00004371 | -0.00216756 ± 0.00055693 | -0.7130 |

The paired covariance matters: the two estimates are anticorrelated by about
0.7 in every generation. No unmarked coordinate exceeds 1.36 standard errors;
the marked current is resolved in each generation.

## Primary fixed-p score

At `p=0.592746050790`, use the declared exact orientation covectors and
`Y_N=(A_M,A_K)`. For the three adjacent wedges

`D_N = A_M(2N) A_K(N) - A_K(2N) A_M(N)`,

the result is

| N | D_N | SE |
|---|---|---|
| 85 | -5.1787132e-6 | 1.0753219e-5 |
| 170 | 4.2032234e-7 | 1.5345959e-6 |
| 340 | 3.1891219e-7 | 3.8439493e-7 |

Full-covariance quadratic form: **2.445811 / 3 df, p=0.485164**.
The correlation-matrix eigenvalues are approximately `0.2554, 1, 1.7446`;
the wedge covariance is not rank deficient.

The covariance-aware common-ray fit `Y_N=a_N(r,1)` gives:

- `r=u_M/u_K=-0.00466903`;
- residual **2.648739 / 3 df, p=0.449009**;
- 95% profile interval for r: **[-0.0222220, 0.0173954]**;
- 99% profile interval: **[-0.0271720, 0.0257492]**.

The submodel `u_M=0` gives **2.859833 / 4 df, p=0.581549**.
Allowing nonzero M loading improves the quadratic form by only **0.211094**
(`p=0.645911` in the asymptotic one-parameter comparison). Non-rejection of
the larger ray model therefore cannot count as positive coupling evidence.

These are measurement-only asymptotic chi-square/profile diagnostics using
estimated batch covariances. They are not exact finite-sample tests or
intervals that include source/model uncertainty.

## Reconstruction and dependence

The raw `tau1,tau2` fields are integer activation ranks, not sampled continuous
Bernoulli threshold times. The source reconstructs `F1,F2` using exact
binomial-tail mixing of each rank histogram, then `M=F1+F2-1`. It retains
direct rank-two births naturally. A step function at `k/N` would be a different
observable and is not used here.

The four archive blocks have **20,80,80,80 batches**, respectively, and are
independent across generations. Within each block the two orientations and
all observer coordinates are paired. The scorer removes one paired batch
from one generation, holds other generations at their full estimates, and
recomputes the entire nonlinear wedge vector. Its covariance is the sum of
the four independent generation contributions. It never invents alignment
between N85 batch 1 and N170 batch 1, nor treats wedge products as independent.

The output retains the selected delete-one vectors, archive hashes, sample
counts, dependency groups, full two-coordinate covariances, and each
generation's contribution to the wedge covariance. All other P337 conclusions
from these same blocks remain correlated reuse, not new independent evidence.

## Sensitivity and exploratory extension

Using the source's per-generation pooled matching roots, each already
recomputed in its delete-one replicates, makes essentially no difference:
wedge p=0.48516 and common-ray p=0.44901.

After the primary result, two structurally motivated exploratory observer
substitutions were examined: `d_eta M` and `E_top`. Neither resolves a nonzero
loading either. Their ray p-values are 0.6753 and 0.9305; these are **not**
substitute endpoints or additional confirmations. All three unmarked rows
share the same archive. The result file preserves both exploratory scores.

## Scientific consequence and next move

This result does not reject a common radial state, and does not show that the
unmarked response is exactly zero. It does show why the marked N680 recurrence
cannot currently be promoted to a statement about M merely because the
common-ray fit passes. A fixed-transfer rejection in this situation could
still be driven entirely by K_A.

The next informative acquisition must resolve the **unmarked loading or its
drift**, not only improve K_A. In particular, another high-N current-only block
does not repair this missing observation. The present data supply no positive
lower bound on the loading, so they do not justify a finite universal sample
budget for identifying it. A geometry with a separately resolved M response,
or a justified variance-reduced estimator of that same M, is a sharper target.

No transfer-model order was frozen before this reconstruction was viewed;
therefore this report does not relabel a post-reveal recurrence fit as a
prospective test. Direct02/plateau conditional-M decomposition and N1360
remain unperformed. #439 should stay open.

## Scientific card

- Mechanism space changed: common-ray survival now has an explicit
  unresolved-loading alternative; no automatic transfer of marked dynamics.
- Not proved: distinct radial states, nonzero shared loading, any CFT identity.
- Observer / sector / source / geometry: canonical unmarked M and natural
  axis activity K_A; exact H4 orientation contrast; square-site birth archives;
  N85/170/340/680 doubling lineage.
- Dependency group: the four P337 source blocks, named in the result JSON;
  no new randomness and no independent replication claim.
- Lifecycle: real-data retrospective score, measurement covariance retained,
  source artifact SHA-pinned, no prospective recurrence claim.
- Next promotion observation: resolved nonzero M loading followed by a
  covariance-aware cross-generation relation or loading drift.

## Reproduce

With NumPy and SciPy, from a clone containing source commit `6123955`:

```sh
python scripts/score_p439_same_stream_production.py
```

Alternatively pass `--input` pointing to the immutable source JSON; the same
SHA-256 is required. The main computation checks the reused covariance against
its delete-one vectors and the F1/F2 identity. No raw Monte Carlo, synthetic
control suite, or broad regression suite was rerun.
