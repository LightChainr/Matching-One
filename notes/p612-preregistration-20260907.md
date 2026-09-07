# #612 pre-registration, recorded before the N725 histograms are read

**Recorded (UTC):** 2026-09-07T03:41Z
**Recipient of the record:** this repository; the N725 artifact is written but not
yet scored at the time of writing.

## 1. Was the scorer fixed before its target was exposed?

**No. This is a reanalysis, not a blind prediction test.**

The target was exposed before the executing team ever touched the data:

- `a0 = -0.0004680842487166532` and its 5% band
  `[-0.0004914884611524859, -0.0004446800362808205]` are written in #612's own
  body, and are the same numbers #609 pre-registered as
  `-4.68084e-04` with band `-4.91489e-04 .. -4.44680e-04`.
- The curvature numbers `+6.70567e-04` (one exponent) and `+1.03745e-03`
  (the 1.547 factor) are in #609.
- The chart transport that revises the second of those to `+0.00103047842`
  is in `docs/astra/ANSWER-610-20260907.md`, already read.

So the executing team has seen every target. Labelling the exercise a blind
forecast would be false. What *is* fixed, and what is being recorded here
before the read, is the **estimator**, not the target:

## 2. What is frozen before the read

Declared now, in the coordinate that will actually be used:

1. **The direction `g`.** The published #582 consensus direction, rebuilt by
   power iteration on the five unit affine residuals exactly as
   `p582_amplitude_law.consensus_direction` does. It is not refit after 725
   arrives, and 725 is not allowed into the consensus.
2. **The quantile grid.** The nine frozen deciles, unchanged.
3. **The weighting.** `spin0` is primary; `equal` is a declared sensitivity
   from the same block, not a second attempt.
4. **The estimator for the 290 -> 725 first difference.** Basis
   `[1, Q_290, g]`, observation `(Q_725 - Q_290)/log(2.5)`, covariance
   `(S_290 + S_725)/log(2.5)^2` — the delete-one jackknife covariance of each
   size, summed, divided by the squared log-step. The number reported is the
   **third coefficient of that fit**, i.e. the `g` coefficient. It is not an
   amplitude read off a differently normalized curve, and it is not the
   un-basis-projected norm of the displacement.
5. **Cross-size covariance.** None. The three productions use distinct seeds
   (`2026105003` for 145, `2026105004` for 290, `2026105011` for 725). Matching
   batch indices do **not** establish coupling: the streams are counter-derived
   SplitMix64 with different seeds, and the systems have different `N`, so the
   permutations diverge. `Var(Q_725 - Q_290) = S_725 + S_290` is therefore the
   declared convention. Where a shared-randomness relationship *is* claimed
   anywhere in this repository it must be demonstrated by a batch-level
   correlation, not by batch IDs.
6. **The decision rule** (from #612, quoted unchanged): a three-standard-error
   measurement interval wholly outside the tolerance band **stops** the
   practical first-amplitude forecast; wholly inside **supports** it; overlap
   with the band boundary is **unresolved**.
7. **The stop rule.** One block, one decision. No larger N, no free second
   exponent, no quantile window chosen after seeing the score, no additional
   production to rescue a failed forecast.

## 3. What is *not* claimed

This is a declared tolerance check with a nominal jackknife error. It is not a
theorem about an exponent, not an exact-coverage test, and not a test of a
second physical scale. Per the chart result, a raw curvature near `1.03e-3` is
already predicted by transporting one exponent through the width contraction
and therefore cannot separate the old "universal 1.547 factor" from a new
scale; the curvature is reported but is **not** the discriminating statistic.

## 4. Production declaration (fixed before the run)

```text
design       N = 725, orientations (26,7) and (23,14)
matrices     [[26,-7],[7,26]] and [[23,-14],[14,23]]
samples      100,000,000 paired configurations  (100 batches x 1,000,000)
             = 100,000,000 per orientation
seed         2026105011   (fresh; 145 used 2026105003, 290 used 2026105004)
replica      7,000,000,000 .. 7,100,000,000
threads      8
```

The corrected sample count matters: `#609` writes "100 batches x 100M", which
would be 100 times this. `n290`'s own metadata says
`samples_per_pair = 100000000`, i.e. 100 batches x 1M. We match `n290`.
