# P537 N65 contact-stage held-out result

## Decision

`CONTACT_FUSION_COMPLETION_TRANSMITS`.

The frozen N25 contact-fusion by topological-completion sign pattern survives
in one held-out N65, 20,000,000-sample production scored as a complete
pooled-root Schur/original-`U` allocation.  No sample extension, descriptor,
distance, size or mixture was added after reading the result.

The collapsed `stage x contact` matrix is

| birth stage | single contact | double contact |
|---|---:|---:|
| `0->1` | `-1.578572908e-7 +/- 2.909771556e-8` | `-9.227657800e-8 +/- 2.233238230e-8` |
| `1->2` | `-3.091298461e-7 +/- 2.143886591e-8` | `+3.696802086e-7 +/- 3.045279085e-8` |

All four marginal 95% intervals exclude zero with the frozen row-major signs
`-- -+`.  The primary determinant is

```text
Delta_cs    = -8.688216055121765e-14
SE          =  1.4139776863031884e-14
95% CI      = [-1.1459612320276013e-13, -5.916819789967515e-14]
Delta / SE  = -6.144521331052191
theta_cs    = -1.0
```

The complete uncollapsed mask matrix, marginal intervals, displacement-level
source means and `beta_y`, and full covariance matrices are in
[`result.json`](result.json).  Covariance is retained separately for the old
P45 100-batch baseline and the new 100-batch MC, then combined only after each
delete-one statistic is recomputed.

## Scientific read

The N25 sign rotation was not a finite-size exact-enumeration curiosity.  At
N65, single and double contact are both negative at first homology birth, while
double contact alone reverses strongly positive at second-direction completion.
Thus the selected transmission map cannot be represented by a scalar contact
counterterm times one scalar birth amplitude.  It contains a reproducible
`contact fusion x topological completion` interaction that reaches the full
original-`U` score.

This result rejects the frozen sign-rotation-null model for this carrier.  It
does not establish that the carrier exhausts the full norm-4 residual, identify
a continuum field or exponent, or turn the nested N25 analyses into independent
evidence.  The exact all-`z` N25 census remains useful as a finite total-closure
decomposition, but it is no longer a prerequisite for transmission: the
prospective N65 gate has already answered that question.

A secondary [full-covariance shape comparison](../p537-contact-shape-transport/latest.md)
asks whether the entire N65 tensor is only one scalar times the exact N25
tensor.  It finds `a=0.0538876+/-0.0025709` and `chi2=10.3371/3`
(`p=.01591`): the two matrices have the same signs and cosine `.9889`, but the
one-amplitude shape is nominally disfavoured.  This question was formulated
after N65 was opened, so it is an adaptive mechanism-compression diagnostic,
not a second prospective decision or independent evidence.

## Frozen production and provenance

- freeze/base commit: `76e2d82e4a1faa76cd71d377fa924a50cb0b6033`
- contract: `analysis/p537_contact_stage_n65_contract.json`
- geometries: `(8,1)` and `(7,4)` with shared counter-derived randomness
- proposal root: `0.5927311266364432`
- new production: 20,000,000 total samples, four deterministic shards,
  100 batches, seed `20260901537`
- external baseline: frozen P45 N65 threshold-rank histogram, 100 batches
- production wall time: 67.61 s; full covariance scoring wall time: 173.63 s
- host: Darwin arm64; Python 3.13.7; Apple clang 21.0.0

Raw shard tables are deterministic, about 115 MB in total, and are not added to
Git history.  Their hashes, row counts and all code/input/output hashes are
recorded in [`run-manifest.json`](run-manifest.json); the committed producer,
seed and frozen inputs regenerate them exactly.
