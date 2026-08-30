# Exact polynomial root certificate

This utility supplies a gauge-free, optimizer-free root certificate for rational threshold-profile polynomials.

It uses exact `Fraction` arithmetic for polynomial division, gcd/square-free reduction, Sturm sequences, root counts, and dyadic isolation on `[0,1]`. Endpoint roots and rational roots encountered during bisection are extracted exactly. Repeated roots are reported once, while derivative sign changes classify isolated stationary points.

For the frozen N=4 synthetic density,

```text
rho(p)=1/2+3p-3p^2,
rho'(p)=3-6p,
```

the only stationary point is the exact root `p=1/2`, classified as a strict maximum.

Synthetic controls cover repeated roots, endpoint roots, no-root polynomials, and a density with three alternating maxima/minima. The isolator fails closed on invalid precision or interval contracts.

## Boundary

This is exact rational-polynomial infrastructure and one synthetic certificate. It does not read production histograms, estimate an empirical mode, fit a tail, or establish a universal scaling function. Issue #28 remains open.

