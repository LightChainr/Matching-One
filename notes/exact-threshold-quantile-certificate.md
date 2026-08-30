# Exact threshold quantile certificate

This certificate applies the merged rational Sturm isolator to the frozen N=4 synthetic threshold CDF

```text
F(p)=p/2+3p^2/2-p^3.
```

For targets `1/4`, `1/2`, and `3/4`, it isolates the unique solution of `F(p)=q` on `[0,1]` with exact dyadic endpoints. Every bracket stores exact CDF endpoint values and verifies the required signs. Non-rational quartiles remain brackets rather than being mislabeled as exact decimals.

The median is the exact rational root `1/2`. The lower and upper quartile brackets reflect exactly about `1/2`; their endpoints therefore give a rigorous rational interval containing the interquartile range. At 24 isolation bits, each nontrivial quantile bracket has width `1/16777216`, and the IQR bracket has width `1/8388608`.

Targets outside `(0,1)`, equations with zero or multiple roots, and nonmonotone synthetic CDF candidates fail closed.

## Boundary

This is an exact synthetic polynomial certificate. It does not read production histograms, estimate empirical quantiles, bootstrap uncertainty, fit tails, or establish universality. Issue #28 remains open.

