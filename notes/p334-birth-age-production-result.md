# P334/P337 result: production-scale birth age predicts the next rank birth

## Main result

At the frozen intrinsic-center layers, the next-step `rank 1 -> rank 2`
hazard retains a strong continuous dependence on `K1` after primitive-line
fixed effects.

The coefficient uses density-normalized age

```text
age=(k0-K1)/N.
```

| archive | orientation | beta_age | SE | two-sided p |
|---|---|---:|---:|---:|
| N325 | first | -0.06725 | 0.00952 | `2.27e-10` |
| N325 | second | -0.08140 | 0.00896 | `1.07e-14` |
| N425 | first | -0.06944 | 0.00848 | `9.67e-13` |
| N425 | second | -0.06692 | 0.00781 | `1.44e-13` |

The shared-randomness, two-orientation joint tests are

```text
N325: chi2=127.856/2, p=1.72e-28
N425: chi2=151.675/2, p=1.16e-33.
```

There are two production evidence blocks, not four: each size's orientations
share one counter-keyed permutation archive.  Negative beta means that, at
fixed current layer and primitive line, older rank-one plateaux are **less**
likely to exit on the next insertion.  Equivalently, recent first births carry
the larger immediate second-birth hazard.  This agrees in sign with the exact
N10 `1/57` witness without reusing it as evidence.

## Recorded birth geometry does not absorb the signal

The secondary diagnostic refines each line by its strict `0->1` birth-local
axis, diagonal, landing and H4 marks.  Across all four primary rows, the largest
absolute coefficient change is only `7.80e-4`, about 1.17% of the corresponding
slope.  The strong signal is therefore not explained by those recorded
birth-local marks.

This is not an intrinsic-memory identification.  The archive does not store
the full microscopic configuration at `k0`; the age coefficient may proxy
unrecorded current geometry.  The correct conclusion is narrower:

> `(k, rank, primitive line)` is not a lumpable predictive state for these
> production streams.

## Complement pair

The exact transform

```text
K1c=N+1-K2,
K2c=N+1-K1,
k0c=N-k0
```

preserves the rank-one risk set and direct-birth diagonal.  All six complement
audit failure sums are zero in both archives.  The transformed age scores are
also resolved jointly:

```text
N325: p=2.60e-30
N425: p=4.94e-24.
```

These are paired views of the same paths, not new evidence rows.

## Direct-birth collision channel

| archive | first D_N | second D_N | orientation mean Dbar +/- SE |
|---|---:|---:|---:|
| N325 | 0.0048895 | 0.0049395 | 0.0049145 +/- 0.0000359 |
| N425 | 0.0039075 | 0.0037980 | 0.0038528 +/- 0.0000315 |

The largest absolute batchwise correlation between `beta_age` and `D_N` is
only 0.0855.  Collision mass is kept as a separate path channel; it is not
subtracted from the birth-age result.

The pre-score conditional adversary

```text
Dbar_425/Dbar_325=(425/325)^(-5/6)=0.7996723
```

is not rejected.  The observed ratio is `0.783956`, with
`chi2=3.285/1`, `p=0.0699`.  This is compatibility, not confirmation: the two
production designs are not a matched fixed-shape lineage, so the test also
contains possible shape-amplitude drift.  Its outcome says nothing for or
against thermal Q4.

## Provenance and boundary

No new simulation was run.  The scorer streamed the two original 2M sparse
tables from the local external-artifact cache after verifying their Huawei
SHA256 hashes.  The raw files remain outside Git; the compact score, protocol,
hashes and complete `10 x 10` delete-one covariance are committed.

This result does not identify a temporal latent variable, Jordan block, CFT
operator, or non-Markov scaling limit.  It establishes a production-scale
failure of coarse-state birth-age independence and points to the next useful
acquisition: record a small current-geometry summary at `k0`, rather than
adding another marginal K1/K2 histogram.

## Scientific card

- Mechanism changed: next-step rank-two hazard is not determined by current
  layer plus primitive line in either production archive.
- Not proved: intrinsic temporal memory or scaling persistence; unrecorded
  current geometry remains a live explanation.
- Observer/source/geometry: P267 N325/N425 two-orientation projective birth
  streams at their frozen intrinsic centers.
- Dependency groups: one paired-orientation archive per size; forward and
  complement views remain paired.
- Next lift: co-record a minimal current-geometry covariate at the observation
  layer and test whether it closes the age coefficient.
