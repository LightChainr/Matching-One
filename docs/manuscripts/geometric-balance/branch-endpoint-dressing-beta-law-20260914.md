# Branch endpoint dressing, singularity order, and the split-fraction Beta law

Date: 2026-09-14

Status: exact convolution/Tauberian algebra plus a mechanism diagnostic for #758/#762/#800. It generalizes the equal-rate exponential/Uniform-split picture in `span-double-pole-branch-convolution-20260914.md`.

## 1. One-sided branch ansatz

Suppose a one-sided excess branch of integer length h has asymptotic weight

```text
b_h
 ~ C/Gamma(alpha) * h^(alpha-1) rho^h,                   (1.1)
```

with `alpha>0` and `0<rho<1`.

Equivalently, under the usual transfer/Tauberian regularity, its generating function has singular part

```text
B(z) ~ C (1-rho z)^(-alpha).                             (1.2)
```

The exponent `alpha` packages endpoint dressing/rooting information in the longitudinal variable. It is not asserted here to be a universal exponent of the SITE model.

## 2. Two-ended convolution

For two asymptotically independent ends with the same `alpha,rho`, the total excess weight is the convolution

```text
d_h = sum_{k=0}^h b_k b_(h-k).
```

Therefore

```text
B(z)^2
 ~ C^2 (1-rho z)^(-2 alpha),                             (2.1)
```

and

```text
d_h
 ~ C^2/Gamma(2 alpha) * h^(2 alpha-1) rho^h.             (2.2)
```

Thus the polynomial/span prefactor and the one-sided endpoint exponent are the same datum in two languages.

Special cases:

```text
alpha=1/2 : d_h ~ const * rho^h          (simple-pole order),
alpha=1   : d_h ~ const * h rho^h        (double pole / Erlang-2),
alpha=3/2 : d_h ~ const * h^2 rho^h      (third-order singularity).
```

The effective near-double-pole pattern in #800 therefore points naturally toward `alpha` near one, subject to direct operator confirmation.

## 3. Conditional branch split

Let

```text
U_h = A_-/(A_-+A_+)
```

conditional on total excess `A_-+A_+=h`.

Using (1.1), for `k=uh` away from the endpoints,

```text
P(A_-=k | A_-+A_+=h)
 proportional to
 k^(alpha-1)(h-k)^(alpha-1).
```

The Riemann-sum limit is

```text
boxed:
U_h => Beta(alpha,alpha),                                 (3.1)
```

with density

```text
u^(alpha-1)(1-u)^(alpha-1)/B(alpha,alpha).
```

Therefore the morphology and the singularity order cross-check one another:

```text
span prefactor h^(2alpha-1)
<=> generating singularity order 2alpha
<=> branch split Beta(alpha,alpha).                       (3.2)
```

## 4. Why alpha=1 is plausible for a free branch endpoint

A fixed-root to **fixed-point** two-dimensional Ornstein--Zernike connection typically has a longitudinal `h^-1/2` prefactor, suggestive of `alpha=1/2` if that fixed endpoint were the branch observable.

A complete-component branch endpoint is not pinned to one transverse site. In the simultaneous regime where the transverse Gaussian spread is `O(sqrt(h))` and is much smaller than the circumference, summing over the `O(sqrt(h))` typical endpoint locations cancels the point-to-point `h^-1/2` local-CLT factor. This heuristically restores

```text
b_h ~ const * rho^h,
```

that is `alpha=1`.

This is a mechanism argument, not a proof for square SITE. Near criticality, endpoint arm insertions can also change amplitudes in the correlation length; those must be separated from the power of h.

## 5. A three-way falsification test

The same `alpha` can be estimated three ways without fitting an arbitrary morphology model:

### Spectrum

Fit the late span response to

```text
h^(2alpha-1) rho^h
```

or the corresponding singularity order.

### Hazard

From (2.2),

```text
gamma_eff(h)
 = gamma - (2alpha-1)/h + O(h^-2).                       (5.1)
```

### Morphology

At fixed large excess-span bins, fit/test

```text
U = A_-/(A_-+A_+) ~ Beta(alpha,alpha).                    (5.2)
```

Agreement of the same alpha across (5.1)--(5.2) would be a substantially stronger mechanism test than a single span exponent.

In particular:

- `alpha=1`: Uniform U and double-pole/Erlang span;
- `alpha=1/2`: arcsine U, endpoint concentration, and simple-pole total span.

This makes the planned #762 U statistic directly informative about endpoint sewing/prefactor structure relevant to #740/#758.

## 6. Claim boundary

- Exact algebra under ansatz (1.1): Sections 2--3 and hazard expansion (5.1).
- Mechanism hypothesis: free transverse endpoint summation yields `alpha=1` for the long SITE branch.
- Not claimed: the one-sided branch factorization, the value of alpha for actual SITE complete components, or equality with a fixed-end two-point OZ amplitude.
