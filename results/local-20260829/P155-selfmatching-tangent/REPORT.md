# P155 exact self-matching odd-tangent analysis

## Result

The exact score-function engine and three frozen tests all point to the same
finite-size conclusion: the five wrapping channels resolve one matching-odd
direction at the available precision. Neither a second wrapping channel nor an
orientation/H4 contrast produced a statistically resolved second response
direction. The source-frozen scalar orthogonalization then transferred from
`N=130` to the previously unseen `N=170` with a residual of only `-0.0244`
standard errors.

This is outcome A in Issue #155: for this topological readout family, the
staggered microscopic tangent is asymptotically parallel to the thermal
response. It does **not** prove that the full microscopic tangent plane is
rank one, and it does not exclude a separate odd irrelevant direction in a
different local or sublattice-sensitive readout.

## Exact engine

The Monte Carlo engine samples only the exact center `p=1/2`. For every
configuration it records the Bernoulli scores

```text
S_t      = 4[(K_e-N_e/2) + (K_o-N_o/2)]
S_lambda = 4[(K_e-N_e/2) - (K_o-N_o/2)]
```

and their covariance with five wrapping observables. Its exhaustive `N=10`
self-test reproduces `d_t R=15/8`, `d_lambda R=5/4`, and the exact Fisher
matrix `4N I_2` for every channel.

## Frozen rank gates at N=130

The first representation was `(a,b)=(11,3)`, with 10,000,000 samples in 100
aligned batches. All ten pairs of wrapping-response rows failed the frozen
gate `condition number <= 50` and `|determinant z| >= 3`. The best determinant
significance was about `1.45 sigma`, while raw condition numbers were of order
thousands or worse. No second generalized eigenvalue was extracted or
regularized.

An orientation rescue used the equal-area representation `(9,7)` with the
same random counters. Its frozen cross matrix combined the orientation mean
and `N` times the `cos(4 theta)` contrast. It had row-angular condition number
`69.94` and determinant significance `-0.098 sigma`. No diagnostic channel
exceeded `0.74 sigma` in absolute determinant significance.

## Prospective scalar holdout

After both matrix gates failed, the scalar protocol was committed before the
`N=170` block was generated. The `N=130` cross channel fixed

```text
c_130 = d_lambda R_cross / d_t R_cross
      = 0.6740950243296661.
```

The held-out `(13,1)`, `N=170` block used the same seed and aligned counter
interval as the source, so every delete-one replicate recomputed the source
ratio and retained cross-size covariance. The result was

| score | estimate | jackknife SE | signed z |
| --- | ---: | ---: | ---: |
| `d_lambda(170) - c_130 d_t(170)` | -0.000153706 | 0.00628772 | -0.02445 |
| `d_t(170) - (170/130)^(3/8) d_t(130)` | -0.00125417 | 0.00440492 | -0.28472 |

The direct coupling ratio was `0.6740655896003787` at `N=170`; its drift from
the frozen source ratio was `-2.94347e-05`, again `-0.02445 sigma`.

## Interpretation and next discriminator

Within wrapping/topological projections, the data favor a single stable odd
response ray rather than a visible thermal-plus-irrelevant two-dimensional
block. More samples of the same rows are therefore low information value. The
next useful experiment is a genuinely local, sublattice-odd readout whose
continuum overlap is not forced to be parallel to wrapping probability. A new
matrix-pencil campaign should be attempted only after that readout passes the
same rank and determinant gates at one source size.

## Reproduction

All production blocks use seed `2026105501`, 100 batches, and the counter
interval `[1000000000,1010000000)`. The source and held-out target contain
10,000,000 samples each. Runtime was 31.0 seconds for `N=130` and 37.9 seconds
for `N=170` on the local Apple Silicon host using one thread. Raw response
tables, metadata, covariance analyses, the prospective score, and checksums
are stored beside this report.
