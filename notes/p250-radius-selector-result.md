# P250 common-counter insertion-radius selector

Status: the preregistered `R={1,2,3,4}` scan is complete.  No candidate passes
the production selector for a two-point-normalized charged cubic.

## Frozen comparison

Every radius used seed `25011312220260901`, counters `[0,4000)`, 40 aligned
batches, the same translated triangles and `p=0.59274605079`.  The field and
translation digests agree batch by batch across all four candidates.  All
annuli are injective in both N325 hands.

The primary gate requires all four `G1/G2 x plus/minus` denominators to have
`|z|>=2` and the aligned eight-real local-variance-normalized cubic covariance
to be invertible at at least two separations.  Phase scores are excluded from
the selector.

| R | d=1 minimum z | d=2 minimum z | d=3 minimum z | usable separations |
|---:|---:|---:|---:|---|
| 1 | 4.507 | 0.480 | 0.0039 | `{1}` |
| 2 | 4.599 | 0.642 | 0.122 | `{1}` |
| 3 | 4.522 | 0.756 | 0.430 | `{1}` |
| 4 | 6.086 | 1.106 | 0.214 | `{1}` |

Increasing the landing radius strengthens some individual pair coordinates,
but never produces a second separation at which all four charged denominators
are resolved.  This is not a sample-count problem at the compact triangle: the
same-counter scan shows that changing only the annular cutoff does not create
the missing mesoscopic two-point tail.

## The apparent R4 exception is not a selector pass

R4 has descriptive local-variance-normalized cubic zero scores

```text
d=1: chi2=13.257/8, p=0.103
d=2: chi2=17.663/8, p=0.0239
d=3: chi2=703.668/8, p=1.16e-146
```

but its weakest two-point z scores at `d=2,3` are only `1.106` and `0.214`.
The separated `G` normalization is therefore unstable exactly where the raw
cubic looks largest.  Selecting R4 from that cubic p-value would change the
question after seeing the data and would not yield the requested normalized
OPE ratio.

This does leave a different, clearly labeled lead: R4 may expose a
local-variance-normalized three-body response whose neutral two-point channel
cancels.  That is not this production target; it would need a fresh frozen
counter and a score formulated without division by `G`.

## Decision

Do not run production for R1--R4 of the current landing-H4 insertion.  The next
useful selector must change the operator, not merely its radius: for example a
leg-defect row or a mesoscopic charged row with a demonstrably resolved
two-point tail at two separations.  Reuse the same support-first score once such
an insertion is defined.

This result does not reject charged OPEs, and it does not claim that the R4 raw
cubic is asymptotic.  It closes one inexpensive geometry family before a large
sample allocation.
