# Frozen intrinsic quantile-center score, N=145 -> 290

## Outcome

The Issue #101 doubling prediction survives its first true-doubling score.
The coordinate and target were committed before the N290 result:

- freeze commit `3762b342b8b376e587df0073044b2c7f6452aa8e`,
  `2026-08-29T06:41:21Z`;
- first target-result commit `9675bce`, `2026-08-29T06:47:18Z`.

For

```text
Q_N = c_0.05 - c_0.025,
Q_290 / Q_145 = 2^(-3/4) = 0.5946035575013605,
```

the observed ratio is `0.5958454898404986` with delta-method standard error
`0.0012034365742889918`.  Equivalently,

```text
Q_145 = -3.7526150840072603e-06
Q_290 = -2.2359787729131497e-06
Q_290 - 2^(-3/4) Q_145 = -4.660494029166079e-09
SE = 4.5130098417110095e-09
z = -1.0326797841413868
chi-square = 1.0664275365743012 / 1 df
```

This is a secondary score.  It does not alter the already frozen primary
P50 scoring order or any primary result.

## Width diagnostic

The dimensionless widths `w_u N^(3/8)` show a small but precisely resolved
positive drift:

| u | N145 | N290 | N290-N145 | SE | z |
|---:|---:|---:|---:|---:|---:|
| 0.025 | 0.0143137596362 | 0.0143295282797 | 1.57686436e-05 | 6.93163854e-07 | 22.7488 |
| 0.050 | 0.0286418746996 | 0.0286736549678 | 3.17802682e-05 | 1.38613036e-06 | 22.9273 |

Thus the leading `N^(-3/8)` width law is accurate at roughly the per-mille
level here but is not correction-free.  This diagnostic is not part of the
one-degree-of-freedom primary quantile-center score.

## Covariance protocol

All four crossings (`-u,+u` at both frozen levels) are solved again inside
every size-local delete-one replicate.  The resulting pseudovalue covariance
therefore carries the nonlinear crossing and shared-curve correlations.
N145 and N290 use independent RNG domains, so their covariance blocks are
combined as independent blocks.  The full size-local and three-residual
covariance matrices are stored in `intrinsic_quantile_center_score.json`.
