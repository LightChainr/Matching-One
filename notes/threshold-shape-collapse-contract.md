# Threshold-distribution shape-collapse contract

Status: methodology-only slice of Issue 122.

The comparison separates three effects before any target curve is scored:

```text
location = q50(target)-q50(reference),
scale    = IQR(target)/IQR(reference),
shape    = standardized quantile residuals.
```

Quantiles use the weighted generalized inverse `min{x:F(x)>=q}`. Standardization uses median center
and interquartile scale. The frozen grid is

```text
[0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95].
```

This makes every positive affine transformation exactly shape-equivalent while leaving non-affine tail
or skew deformations visible. Inputs with a zero IQR fail closed rather than receiving an arbitrary
fallback scale.

The committed examples are synthetic fixtures only. They validate the contract but provide no evidence
for collapse across sizes, geometries, or microscopic models. A production analysis still needs a
frozen covariance/bootstrap plan and independently held-out curves before assigning significance.
