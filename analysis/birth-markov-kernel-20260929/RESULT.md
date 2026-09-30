# Does rank one retain birth-history memory?

2026-09-29; retrospective reuse of the completed 560k cloud block. No new samples.

One contrast, selected before this scoring: cutoffs a,b,c are the .25/.50/.75 lower quantiles of the pooled J1/J2 births.
Among paths with J1<=b<J2, compare P(J2>c) for early (J1<=a) versus late (a<J1<=b) first birth.
Every current-rank-only time-inhomogeneous Markov process predicts difference zero.

## Result

L512 is the primary available-size readout; smaller sizes supply context. Errors are 14-batch delete-one SE, not confidence intervals.
All thresholds are re-estimated in every deletion; the JSON retains their values and the complete metric covariance.

| Lattice | L | Rank-one risk count | Early survival | Late survival | Difference (percentage points) |
|---|---:|---:|---:|---:|---:|
| square | 64 | 29042/70000 | 0.4649 | 0.4524 | +1.248 ± 0.553 |
| square | 128 | 28218/70000 | 0.4516 | 0.4437 | +0.785 ± 0.480 |
| square | 256 | 27775/70000 | 0.4578 | 0.4438 | +1.395 ± 0.623 |
| square | 512 | 27537/70000 | 0.4491 | 0.4349 | +1.422 ± 0.511 |
| triangular | 64 | 28622/70000 | 0.4493 | 0.4430 | +0.623 ± 0.521 |
| triangular | 128 | 27429/70000 | 0.4427 | 0.4391 | +0.352 ± 0.712 |
| triangular | 256 | 26832/70000 | 0.4428 | 0.4323 | +1.054 ± 0.476 |
| triangular | 512 | 26681/70000 | 0.4325 | 0.4255 | +0.702 ± 0.515 |

## Existing exact L3 archive, a different evidential object

No 9! enumeration rerun: the stored complete tables are evaluated with integer/Fraction arithmetic.
All 120 count-time triples 0<=a<b<c<=9 are checked here. This finite exhaustive identity check is separate from the ONE fixed production triple; it is not a production cutoff search.

- square L3: 0 violating triples out of 120; rank-only insertion process Markov under the full-kernel criterion: True.
  - Exact iid-label bridge at p=(1/3,1/2,2/3): contrast 0 (+0.0000 percentage points); Markov identity violated: False.
- triangular L3: 2 violating triples out of 120; rank-only insertion process Markov under the full-kernel criterion: False.
  - [3, 4, 5]: early/late survival 3/5 / 13/25; difference 2/25; conditional covariance 1/90.
  - [4, 5, 6]: early/late survival 3/16 / 1/7; difference 5/112; conditional covariance 1/90.
  - Exact iid-label bridge at p=(1/3,1/2,2/3): contrast -1156/716639 (-0.1613 percentage points); Markov identity violated: True.

The exact clock bridge is a multinomial mixture, not a replacement of count k by its expectation Np. Insertion and label clocks have separately evaluated properties; neither finite result is extrapolated to larger sizes or a continuum limit.

## Scope

- Not prospective or independent model-elimination certification; this block already supplied gap readouts.
- Finite-size conditional dependence does not prove non-Markovianity of a continuum limit.
- Conditioning on rank one tests the Markov property; it is not a causal effect of earlier birth.
- A zero contrast at one triple would not prove Markovianity; the theorem requires all triples.
- Square and triangular primitive-coordinate tori have different moduli; no universality comparison.
- No fitted exponent, pooled cross-lattice significance, adaptive cutoff search or hidden-state identification.

Source: `../birth-gap-20260929/cloud-data-5k/run.json`; all 112 production batches only.
No local-pilot pooling or benchmark rows. Grain: one uniform-permutation filtration, represented by weighted (J1,J2) histogram rows.
One run took 2.74 s on Python 3.11.15 (arm64).
One exact four-cell control is recorded in result.json; old physical enumeration was not rerun.

Reproduce from the repository root:

```sh
python3 analysis/birth-markov-kernel-20260929/analyze.py
```
