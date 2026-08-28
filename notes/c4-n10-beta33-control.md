# Exact Beta(3,3) threshold distribution on the N=10 C4 self-matching control

The exhaustive `(a,b)=(3,1)` checkerboard-triangulation control in PR #82 has

```text
M Bernstein sums by occupation k=0..10:
[-1,-10,-45,-100,-100,0,100,100,45,10,1].
```

Converting these exact Bernstein sums to the power basis gives

\[
M_{10}(p)=12p^5-30p^4+20p^3-1.
\]

Equivalently,

\[
\boxed{M_{10}(p)=2I_p(3,3)-1},
\]

where `I_p(3,3)` is the regularized incomplete beta function, i.e. the CDF of a `Beta(3,3)` random variable.

Thus the threshold-distribution interpretation

\[
F(p)=\frac{1+M(p)}2
\]

is exactly Beta(3,3) on this finite quotient.

## Persistent self-matching factor

The same polynomial factors over `Q` as

\[
\boxed{
M_{10}(p)=(2p-1)
\left(6p^4-12p^3+4p^2+2p+1\right).
}
\]

The linear factor is structural, not an N=10 accident. For any finite self-matching quotient whose matching observable obeys the exact complement antisymmetry

\[
M(1-p)=-M(p),
\]

one has `M(1/2)=0`; hence every rational-coefficient polynomial representation contains the factor `2p-1`.

This provides the positive algebraic control required by the finite-factor/GCD program: an exact self-matching mechanism produces a persistent physical linear factor, whereas the square-site axis/diamond target polynomials certified so far are irreducible over `Q` and do not contain a bounded-degree physical factor.

Immediate exact consequences at the self-matching center `p=1/2` are

\[
M'(1/2)=\frac{15}{4},\qquad
M'''(1/2)=-60,\qquad
M^{(5)}(1/2)=1440,
\]

and the dimensionless derivative ratios are

\[
\kappa_3=-\frac{256}{225},\qquad
\kappa_5=\frac{32768}{16875}.
\]

These values are finite-size control constants, not universal limits and not evidence for the `-5/3` conjecture. Their value is methodological: they give a compact exact regression for

- Bernstein-to-canonical reconstruction;
- threshold-rank CDF semantics;
- analytic derivatives at the center;
- `kappa_3` / `kappa_5` code paths;
- self-matching central parity;
- recovery of the exact persistent self-matching factor `2p-1`.

`tests/test_c4_self_matching_exact.py` recomputes the polynomial identity from the exhaustive integer coefficients without using floating-point special functions.
