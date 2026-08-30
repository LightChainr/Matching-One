# P321 fresh multiscale decision freeze

Status: **frozen before production; no listed fresh counter has been run**.

## Why 100k per shape

The exploratory 20k-per-shape multiscale pilot completed the already fixed
`N^-2/N^-3` covariance scorer:

```text
scale fit:       chi2 = 1.490179138584335 / 4
E4 residual:     chi2 = 7.531750170483802 / 3
E4 pilot p:      0.05674824191100037
```

The E4 residual vector in the frozen order `(16/9,9/4,4)` is

```text
(83.23718955539628, 60.65219308159848, 100.27195086427044).
```

These numbers are exploratory.  They choose precision, not the curve,
exponents or decision threshold.  With the conservative planning estimate

```text
lambda_20k = max(7.531750170483802 - 3, 0)
            = 4.531750170483802,
```

five times the samples gives `lambda_100k=22.65875085241901`.  At
`alpha=.01`, the fixed three-degree-of-freedom critical value is
`11.34486673014437`, so the noncentral-chi-square planning power is
`0.9498366103759819` and the expected score is `25.65875085241901`.  A 60k
design would have power only `0.7279526360658312`.  This freezes 100k per
shape, not an optional or sequential target.

The pilot artifact used for this decision has SHA-256

```text
d425d2f66025c54e5f7d9faa310df004d899d07bca9d2459ee3c02a8ecf1cbe6
```

and is never pooled with the fresh score.

## Fresh streams

RNG domains follow `matching-one-rng-domain-v1`.  Shapes share one common
field within each `N`; distinct `N` values are independently derived from the
tag `P321-fresh-multiscale-100k-v1` and base seed `20260830321`.

| N | effective seed | counter interval | batches |
|---:|---:|---:|---:|
| 144 | 6512046139770053406 | `[100000,200000)` | 50 |
| 576 | 9111365054197340087 | `[200000,300000)` | 50 |
| 1296 | 15228459834092891587 | `[300000,400000)` | 50 |

Every batch contains 2,000 replicas.  At fixed `N`, the four square/rectangle
pair invocations use the identical seed, interval and batch boundaries, so the
repeated square rows must be byte-identical and the complete aligned `5x5`
root covariance is retained.  Cross-size covariance is block diagonal.

All seeds differ from the 20k pilot seeds and every listed counter interval is
new and disjoint.  There is no continuation, overlap or pilot augmentation.

## Frozen score and decision

No exponent or shape is selected on the fresh data.  The only scale model is

```text
p(N,rho)=pc+C_N(rho) N^-2+D_N(rho) N^-3.
```

Convert `C_N` to the transverse-width convention before scoring:

```text
C_width(rho)=C_N(rho)/rho^2.
```

The ordinary thermal-Q4 curve remains the `65b3830` E4 oracle.  The three
parameter-independent primary ratios are frozen at

| rho | `C_width(rho)/C_width(1)` |
|---:|---:|
| 16/9 | 0.6892481080041541 |
| 9/4 | 0.6870445687187882 |
| 4 | 0.6869250536811942 |

`rho=9` is an endpoint diagnostic and never enters a primary decision.

Two scores are reported in order:

1. The four-degree-of-freedom fixed-scale fit is adequate at `alpha=.01` iff
   its chi-square is at most `13.276704135987622`.  If it fails, archive the E4
   score but classify the result as scale-law failure rather than E4-shape
   rejection.
2. Conditional on the scale gate, reject the frozen E4 shape iff its
   three-degree-of-freedom chi-square exceeds `11.34486673014437`; otherwise
   do not reject it.

Rejecting E4 does not validate the pilot residual vector as a new curve.  The
pilot vector is only the explicit planning alternative.  There is no E4
amplitude fit, exponent fit, curve refit, interim look or optional stopping.

## Launch boundary

This branch freezes the decision only.  It does not execute the engine.  The
root agent may launch after the exact commit containing this note and manifest
has been pushed and posted publicly to Issue 321.  Every engine invocation
must use that public commit in `--git-commit`, 100,000 samples, 50 batches and
the exact domain row above.
