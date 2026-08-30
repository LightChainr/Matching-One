# Exact direct-rank2 mechanism decomposition

This control addresses the mechanism diagnostic in Issue 439. At one exact
threshold, every orientation's unmarked Matching-One numerator is partitioned
by the immutable birth-path tag:

```text
P2-P0 = (P2_direct-P0_direct) + (P2_plateau-P0_plateau).
```

Both terms retain the original total-sample denominator, so their orientation
contrasts add exactly to `A_M`. The certificate freezes

```text
A_M_total    = -1/12
A_M_direct   = -1/6
A_M_plateau  =  1/12.
```

It separately reports the cumulative direct `0->2` birth fraction and
`M_without_direct02`, whose denominator contains only plateau rows. That last
quantity is descriptive conditioning and is deliberately not used in the
additive closure or named as a new Matching-One observable.

The parser rejects malformed birth rows, equal covectors, unequal orientation
sample counts, missing direct or plateau sectors, and inexact floating inputs.

## Boundary

This exact synthetic certificate does not estimate production direct-birth
rates, import raw archives, propagate covariance, assign scaling exponents, or
identify a continuum mechanism. Issue 439 remains open.
