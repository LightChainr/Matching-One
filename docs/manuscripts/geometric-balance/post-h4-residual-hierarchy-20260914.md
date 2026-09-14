# What remains after the leading H4 matching-root correction is projected out?

Date: 2026-09-14

Status: synthesis / conjecture ledger for the **next** finite-size mechanism.  The leading H4/spin-four root displacement is already much more strongly constrained than the residual discussed here.  This note deliberately starts only after the same-circle H4 projector has been applied.

## 1. Empirical object

At one fixed physical circumference `ell`, let two orientations have

```text
h_i=H4(theta_i)=cos(4 theta_i),
p_i=p_ch(ell,theta_i).
```

Define the H4-annihilating scalar projection

```text
p_perp4
 = (h1 p2-h2 p1)/(h1-h2).
```

The finite residual relative to the infinite threshold is

```text
r(ell;theta1,theta2)=p_perp4-pc.
```

Existing deterministic controls give, after comparison to the diagnostic high-precision threshold only **after** forming the projector,

```text
ell=5, axis/(3,4):        r≈-9.81e-6
ell=sqrt(65), (1,8)/(4,7): r≈+5.66e-8
ell=sqrt(85), (2,9)/(6,7): r≈+4.33e-10
ell=10, axis/(3,4):       r≈-7.57e-8.
```

The N85 value is extraordinarily small but should not be interpreted by itself as an asymptotic error law.

## 2. The H4 projector removes the whole H4 tower

Write the fixed-ell D4 harmonic decomposition

```text
p_ch(ell,theta)
 = pc + S0(ell)
      + H4(theta) S4(ell)
      + H8(theta) S8(ell)
      + H12(theta) S12(ell)
      + ... .
```

The two-angle projection removes **all of `S4(ell)`**, independent of how many powers of `ell^-1` occur inside that coefficient.  Therefore an axial correction ladder

```text
ell^-4, ell^-6, ...
```

can disappear wholesale if those terms are successive corrections carrying the same H4 irrep.

The residual is

```text
r = S0 + C8 S8 + C12 S12 + ...,
```

where for the two-angle H4 projector

```text
C8=-(1+2 h1 h2).
```

This shifts the research question from “what is the next power of the axial root?” to “what angular irrep survives after the dominant one is removed?”

## 3. The simplest coordinate-artifact explanation is far too small

The TANGENT_SPIN4 picture says the leading H4 correction acts almost as a translation of the thermal scaling coordinate.  A nonlinear thermal field can then generate an `H4^2=(1+H8)/2` term even if no independent H8 operator exists.

The existing full-curve analysis calibrates the quadratic logit thermal metric as

```text
t = delta h + c2 (delta h)^2+...,
c2≈-0.0311.
```

The leading root shift has

```text
delta p≈-A_p H4 ell^-4,
A_p≈0.296--0.300.
```

At the threshold

```text
dp/dh=pq≈0.2414,
```

so the leading logit shift amplitude is `A_h≈A_p/(pq)≈1.2`.  Inverting the quadratic thermal field and then the logistic map gives only

```text
delta p_coord^(2)
 = O(1e-2) H4^2 ell^-8.
```

Thus the induced scalar/H8 pieces are roughly

```text
O(1e-10) at ell~8--10,
```

whereas N65 and ell=10 H4-projected residuals are `O(1e-8--1e-7)`.

Conclusion:

> The observed post-H4 residual cannot be explained by the already calibrated ordinary quadratic thermal-coordinate nonlinearity of the leading tangent spin-four shift.

More elaborate normal spin-four corrections or logarithmic mixing remain possible; only the cheapest analytic-coordinate explanation is excluded.

## 4. The residual already shows an H8 geometry pattern

For the three projectors with `ell` in the relatively narrow range `8.06--10`, the H8 weights are

```text
N65:   C8=-0.1484319457,  r=+5.66e-8
N85:   C8=+0.2224938303,  r=+4.33e-10
ell10: C8=+0.6864,        r=-7.57e-8.
```

If one ignores the modest ell variation only as an exploratory local diagnostic and writes

```text
r≈S0+C8 S8,
```

a least-squares line gives approximately

```text
S0≈+3.4e-8,
S8≈-1.59e-7,
```

with pointwise discrepancies only around `1e-9`.

This should **not** be promoted to a fitted scaling law: the three circumferences differ.  Its value is qualitative.  The sign progression and the near-zero N85 residual are naturally explained by a negative H8 contribution crossing a smaller scalar/other residual.

In this picture N85 is a geometry cancellation point, not miraculous proof of a `4e-10` asymptotic estimator error.

## 5. Candidate NEXT_H8_IDENTITY_x8

A structurally natural next angular field is an identity-family spin-eight lattice anisotropy with

```text
x_8=8,
spin=8,
angular harmonic H8=cos(8theta).
```

If the two microscopic safe regularizations have a nonzero difference coupling to this field, its magnetic-gap correction scales as

```text
Delta I_8 ~ H8 ell^(1-x8)=H8 ell^-7.
```

Dividing by the thermal root susceptibility `ell^-1/4` gives

```text
boxed:
delta p_8 ~ H8 ell^-27/4,
27/4=6.75.
```

This is distinct from a thermal-family spin-eight descendant, which would have `x=x_t+8=37/4` and root exponent `8`.

Why this candidate is interesting:

1. the raw axial `ell^-4,ell^-6,...` corrections may remain in the H4 tower and are removed by the H4 projector;
2. the same-angle axis/(3,4) projected residual drops by about a factor `129` between ell=5 and ell=10, corresponding to a naive two-point effective power near `7`; this is **not** an exponent measurement, but it is numerically compatible with `27/4`;
3. restoring the exploratory H8 piece around ell~9 with an `ell^-27/4` law gives an O(1) amplitude (`roughly 0.4--0.6`), whereas an `ell^-8` law requires a substantially larger amplitude;
4. N65 and ell10 residual signs agree with a negative H8 coefficient because their projected `C8` values have opposite signs.

This is presently a conjecture, not a field identification.  In particular, the fact that lower identity-family corrections are common/sector-even does not force every higher lattice coupling to have zero NN/matching difference.

## 6. Competing residual mechanisms

### 6.1 Angular scalar / spin-zero channel

A true angle-independent dual-odd correction survives every angular projector.  It could come from a scalar irrelevant/logarithmic sector or another microscopic coupling.  Same-circle tomography measures it as the H0 coefficient.

Its existence would **not** restore scalar x=21/4 as the leading explanation of the raw root shift; it would be a subleading channel exposed only after H4 removal.

### 6.2 Thermal-family spin eight

This gives H8 with root exponent 8.  N1105 can identify the H8 irrep at one ell but cannot distinguish exponent 6.75 from 8 by itself.

### 6.3 H12 / higher D4 harmonics

A two-angle H4 projector is generally sensitive to all higher harmonics.  The N1105 four-angle closure contrast has unit-normalized H12 gain about `0.705`, providing a clean detector if H0/H4/H8 are insufficient.

### 6.4 Logarithmic energy--hull mixing

At Q=1, logarithmic collisions can decorate an angular power with `log ell` and change apparent effective amplitudes.  This requires generic-Q/module evidence; it should not be inferred from three residual values.

### 6.5 Normal component of the leading spin-four field

Root/slope-normalized same-ell curves differ only at `O(1e-5--5e-5)` on the tested interval, so any leading normal H4 shape response is small at current sizes.  A higher-order H4 normal piece is nonetheless possible; importantly, the root H4 projector removes it from the scalar root estimator if it carries the same H4 angular character.

## 7. Two targeted discriminants

### N1105: direct same-circle H0/H4/H8/H12 tomography

Four primitive orientations at `N=1105` allow either

```text
H0/H4/H8 fit + one exact closure residual,
```

or exact four-column `H0/H4/H8/H12` decomposition.  This is the cleanest way to decide whether the present residual is actually H8 or scalar without cross-size assumptions.

### N377: H4/H8 double-notch

The same circle

```text
377=4^2+19^2=11^2+16^2
```

has two primitive orientations with

```text
C8=+0.0036146395.
```

Thus the ordinary two-angle H4 projector retains only about `0.36%` of H8 while remaining sensitive to later harmonics (`C12≈-0.1378`, `C16≈-0.9811`).  If computationally feasible, it is a direct scalar/H12-vs-H8 stress test.

A state-cap pilot shows its width-one safe automata exceed 200k states, so it is not a default cheap run.

## 8. Decision logic

### N1105 finds dominant P8, small P0/P12

Promote the next-channel picture to

```text
H4 leading regularization mismatch
 -> H8 subleading lattice anisotropy.
```

Then use existing different-ell projected data only as a weak exponent cross-check between identity-spin8 (`27/4`) and thermal-spin8 (`8`).

### P0 scalar dominates

A genuine subleading angle-independent channel survives.  Reopen its field/module classification, but keep it distinct from the already resolved leading H4 mechanism.

### P12 or closure residual large

The two-harmonic residual model is inadequate.  Keep a D4 tower; do not force an exponent fit.

### N377 remains large while N1105 says P8 dominates

Then the assumption that N377 nearly nulls the relevant H8 component is wrong, indicating either higher harmonics or a geometry-dependent non-Fourier nuisance.

## 9. Claim boundary

- H4-projected residuals, H8 geometry coefficients and thermal-metric coefficient are existing deterministic controls / algebra.
- The exclusion of the simple quadratic-coordinate explanation is an order-of-magnitude consequence of those measured coefficients.
- `NEXT_H8_IDENTITY_x8` is a new falsifiable conjecture.
- The exploratory local `S0+C8 S8` decomposition is not a same-ell fit and must not be treated as evidence of an asymptotic amplitude.
- N1105/N377 are mechanism-discriminating designs; only N1105 currently has an opened bounded compute issue (#807), with an explicit Phase-0 stop rule.
