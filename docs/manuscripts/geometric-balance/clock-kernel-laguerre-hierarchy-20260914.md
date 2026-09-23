# Laguerre hierarchy for the intensity-clock gap process

Date: 2026-09-14

Status: exact consequence of the abstract Poisson intensity-clock semigroup in `clock-kernel-separation-20260914.md`. This is an unmarked fixed-location process result, not an additional SITE convergence theorem.

## 1. One-sided stationary process

Let `U_x` be the normalized nearest-barrier distance on one side at log intensity `x=log Lambda`. For `h>=0`, put

```text
q = exp(-h),
r = exp(h)=1/q.
```

The two-time Laplace transform is

```text
E exp[-s U_x - t U_{x+h}]
 = (r+s)/[(1+s)(s+r(1+t))].                              (1.1)
```

Let `L_n` be the ordinary Laguerre polynomial, orthonormal under the stationary `Exp(1)` measure:

```text
int_0^infinity e^-u L_n(u)L_m(u) du = delta_nm.
```

Using the generating function

```text
sum_{n>=0} L_n(u) z^n
 = (1-z)^-1 exp[-u z/(1-z)],                              (1.2)
```

substitution of (1.1) gives

```text
E [sum_n L_n(U_x) z^n] [sum_m L_m(U_{x+h}) w^m]
 = [1-(1-q)z] / [1-(1-q)z-q z w].                        (1.3)
```

Therefore, for `n,m>=1`,

```text
E[L_n(U_x)L_m(U_{x+h})]
 = 0,                                                      m>n,

 = binom(n-1,m-1) q^m (1-q)^(n-m),                        n>=m.   (1.4)
```

In particular the diagonal mode correlations are exactly

```text
E[L_n(U_x)L_n(U_{x+h})] = exp(-n h).                      (1.5)
```

The process is not reversible: off-diagonal entries occur only on one side of the Laguerre matrix. Time reversal flips this triangular direction.

## 2. Two-sided fixed-location gap

Let

```text
S_x = U_x^- + U_x^+,
```

so the stationary law is `Gamma(shape=2,rate=1)`. Let `L_n^(1)` be the generalized Laguerre polynomial. Its norm is

```text
int_0^infinity u e^-u [L_n^(1)(u)]^2 du = n+1.
```

Because the left/right processes are independent, the bivariate generating function is the square of (1.3). Hence for `n,m>=1`,

```text
E[L_n^(1)(S_x)L_m^(1)(S_{x+h})]
 = 0,                                                      m>n,

 = (m+1) binom(n-1,m-1)
   q^m (1-q)^(n-m),                                       n>=m.  (2.1)
```

The normalized diagonal correlation of the n-th Gamma-Laguerre mode is again

```text
exp(-n h).                                                 (2.2)
```

The familiar ordinary length correlation `exp(-h)` is only the `n=1` member of this hierarchy.

## 3. Why this is useful for #780

If the actual common-label SITE barrier process converges to the clock-kernel limit, then after calibrating the intensity ratio no free parameter remains in (1.4) or (2.1). A two- or three-parameter finite-width check can therefore test substantially more than the first covariance:

```text
mode n=1  -> exp(-h),
mode n=2  -> exp(-2h),
mode n=3  -> exp(-3h),
```

plus the asymmetric triangular cross-mode coefficients.

This is especially useful for distinguishing:

- a true record/splitting process;
- a reversible Markov surrogate with the same one-time Gamma law;
- a process with barrier mergers, which should create forbidden/reweighted cross-mode entries.

No new sampling is requested by this note; the formulas are regression targets if a paired common-label calculation is already performed.

## 4. Claim boundary

- Exact: (1.3)--(2.2) for the abstract Poisson clock process.
- Conditional: transfer of the hierarchy to square SITE under the #780 process convergence assumptions.
