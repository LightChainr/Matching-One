# Minimal even rank-2 transfer after the norm-4 reveal

## Result used

The frozen Issue #154 score rejects the ordinary q=2 scalar and full-jet
closures, while the rank-2 Jordan versions land at `p=0.0673` and `p=0.0543`.
This note asks the smallest post-reveal question that the submitted sufficient
statistics can answer: are the two lineages' remaining Jordan curvatures
compatible with **one shared conjugation-even mode direction**?

Write the width-normalized thermal jet at generation `k` as

```text
x_k = a + b k + c lambda^k v.
```

The fixed `a+b k` term is the inherited rank-2 Jordan direction. Its second
difference vanishes. The one extra even mode gives

```text
s_0 = x_2 - 2 x_1 + x_0 = c (lambda-1)^2 v.
```

With only three generations, `lambda` is not identifiable: it is absorbed
into `c`. What *is* identifiable is the rank-one prediction that the two
lineage residual vectors are proportional across thermal orders 2 through 6.

## Covariance-aware rank-one fit

The model uses six fitted quantities for ten correlated residuals: five
entries of `v` (anchored to the 65-lineage amplitude) and one 85/65 amplitude
ratio. The full frozen 10x10 Jordan-residual covariance is used.

```text
amplitude ratio rho_85/65 = 0.177756
68% profile interval       = [-0.4467, 1.0623]
95% profile interval       = [-1.1866, 2.7847]
rank-one fit               = chi2(4) = 4.7563, p = 0.3132
improvement over zero mode = Delta chi2(6) = 13.286, exploratory p = 0.0387
```

The fitted shared direction over orders `(2,3,4,5,6)` is

```text
(0.13673, 0.25944, -0.27415, -2.95931, -0.79924).
```

This is a viable compression of the residual, not a precise amplitude law:
the amplitude-ratio interval is wide and includes zero. The high-information
result is therefore structural: the existing covariance does not require a
free matrix or one correction per jet order. A single shared residual
direction is enough at current precision.

## Frozen next-generation fork

To make the unidentifiable eigenvalue falsifiable, freeze the most economical
existing value rather than fit it: `lambda=1/2`, inherited from the ordinary
analytic q=2 irrelevant mode. This does **not** revive the rejected pure-q=2
model; it places that mode below the Jordan leading block.

For each lineage the future null is

```text
x_3 - 2 x_2 + x_1
  - lambda (x_2 - 2 x_1 + x_0) = 0.
```

The machine-readable artifact freezes points for N520 and N680. Its scalar-U
side predictions are

```text
U520 = 1.9279025735
U680 = 2.5725311898
```

The corresponding thermal-jet predictions are recorded order by order in
`predictions/norm4_even_rank2_next_generation_20260829.json`.

The fourth generation separates mechanisms without a free exponent fit:

| next/current second difference | mechanism |
|---:|---|
| 0 | pure rank-2 Jordan has been reached |
| 1/2 | Jordan plus analytic even irrelevant mode |
| 1 | persistent curvature; rank-3/log-squared candidate |
| another common ratio | a different coherent even eigenmode |
| no common ratio | one-secondary-mode transfer is false |

## Evidence boundary

The rank-one fit is a post-reveal mechanism inference. It is not new primary
evidence and its nested improvement p-value is descriptive. The N520/N680
recurrence is a real frozen prediction for any future block. Scalar U is not
added to the jet chi-square because the submitted scalar and jet JSON files do
not contain their cross-covariance.

The bolder continuation is clear: if the half-eigenvalue recurrence survives,
carry the same two-dimensional even state into the norm-10 phase matrix and
test Gaussian-semigroup composition. If it fails, do not add more scalar
exponents; promote the minimal source-stable matrix under Issue #180.
