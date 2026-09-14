# Post-H4 exponent competition: linear scalar `7` versus quadratic H4 composite `29/4`

Date: 2026-09-14

Status: scaling-mechanism derivation / model challenge.  It adds a missing post-H4 competitor to the historical `~L^-7` discussion.  The exponent algebra is standard RG bookkeeping once the leading H4 correction exponent is granted; the nonzero quadratic coefficient is a hypothesis to be tested.

## 1. Leading H4 correction fixes its RG exponent

The observed leading square-lattice charge-root correction behaves as

```text
p_root-pc ~ L^-4 H4(theta),
H4(theta)=cos(4theta).
```

Let the corresponding irrelevant coupling have correction exponent `omega_4>0`.  The thermal scaling exponent is

```text
y_t=3/4.
```

For a dimensionless balance/root equation, a first-order irrelevant correction shifts the thermal coordinate by

```text
t_root ~ u_4 L^[-(y_t+omega_4)].
```

Therefore the observed root exponent four implies

```text
omega_4 = 4-y_t = 13/4.
```

Equivalently the associated total scaling dimension is

```text
x_4=2+omega_4=21/4.
```

This is compatible with the current H4/spin-four shell interpretation but the argument below uses only the measured exponent/angular sector.

## 2. Second order in the same H4 coupling predicts root exponent `29/4`

At second order, the RG expansion contains terms quadratic in the same irrelevant coupling:

```text
u_4^2 L^(-2 omega_4).
```

Balancing against the thermal scaling variable gives

```text
t_root^(2)
 ~ u_4^2 L^[-(y_t+2 omega_4)].
```

Substituting

```text
y_t=3/4,
omega_4=13/4
```

gives

```text
boxed:
y_t+2 omega_4
 = 3/4+26/4
 = 29/4
 = 7.25.
```

Thus an entirely natural post-H4 correction is

```text
p_root^(2)-pc ~ L^-29/4,
```

provided the connected quadratic coefficient is nonzero.

This exponent was absent from the older #47 fixed-model challenge, which tested `11/2,6,7,8,10` plus a free power/log models.

## 3. Why the historical `~7.06` does not distinguish `7` from `7.25`

Mertens--Ziff observed an accelerated-root effective power near seven at small available sizes and explicitly warned that finite-size differences strongly distort apparent exponents.

Two current mechanisms are therefore close enough that small-size effective exponents cannot decide them safely:

```text
linear scalar candidate:
  root exponent = 7,
  omega_s = 7-y_t = 25/4,
  x_s = 2+omega_s = 33/4;

quadratic H4 composite:
  root exponent = 29/4 = 7.25.
```

A fitted `7.06` is compatible with substantial crossover between them or with one of them plus lower-order radial corrections/logs.

Hence future radial challenges should include **29/4 as a preregistered fixed exponent**, not leave it hidden inside the free-power model.

## 4. Angular structure of the quadratic H4 term: H0 and H8

Represent the real H4 lattice coupling using complex spin components

```text
u_+ ~ a exp(+i4 theta),
u_- ~ a exp(-i4 theta).
```

At second order there are two distinct angular tensor products:

```text
u_+ nu_- : spin 0 / H0,
nu_+^2 + nu_-^2 : spin +/-8 / H8.
```

Therefore the generic quadratic contribution has the structure

```text
L^-29/4 [
    B_0^(2)
  + B_8^(2) cos(8theta)
]
```

at root level (up to logs/contact mixing and other sources).

The scalar and H8 coefficients are **not required to be equal**; they are controlled by different connected/OPE/contact channels.  But the presence of H8 as a natural companion is a useful discriminator.

This is exactly why an H4-null projector does not eliminate the quadratic mechanism: it kills the linear H4 term, while H0/H8 survive.

## 5. N377 has unusually good leverage on this mechanism

The #808 pair has

```text
C4=0,
C8~=0.00361464,
C12~=-0.137808,
C16~=-0.98105.
```

Thus it almost removes the H8 companion of the quadratic H4 mechanism while preserving H0.

Consequently, if the quadratic mechanism dominates and `B_0^(2)` is nonzero, the N377 H4-null output can look very nearly scalar even though **no new linear scalar primary is present**.

This is a direct warning against interpreting a successful N377 `L^-7-ish` residual as V14 without a source-order test.

## 6. Strong discriminator: tune the microscopic H4 coupling

Suppose a local microscopic parameter `lambda` changes the leading H4 coupling through zero:

```text
u_4(lambda*)=0.
```

Then near the improved point,

```text
linear H4 amplitude ~ (lambda-lambda*),
quadratic H0/H8 composite ~ (lambda-lambda*)^2.
```

A genuinely independent linear scalar coupling generally need not vanish at `lambda*`.

Therefore an improved-action family is the cleanest discriminator:

```text
post-H4 scalar residual versus measured leading H4 amplitude.
```

If the residual scales quadratically with the leading H4 coefficient and collapses near its zero, the composite mechanism is strongly supported.

If the H4 amplitude crosses/tunes small while the scalar `L^-7` residual persists with nonzero intercept, an independent scalar channel is required.

## 7. Interface to PR #148 self-matching checkerboard family

PR #148 constructs the exact local family

```text
p_even=1/2+t+lambda,
p_odd =1/2+t-lambda,
```

with exact pair exchange `(t,lambda)->(-t,-lambda)`.  Its planned nontrivial improved-action search explicitly targets the **exchange-even H4 amplitude** `A_T4^+(N,lambda)` on same-norm orientation pairs.

That family is therefore conceptually ideal for the discriminator above **if** a nonzero H4-amplitude zero is found and transports across size.

Important boundary: PR #148's exact minimum quotient only proves that the exchange-odd response has no nonzero legal zero; it does not yet establish a nonzero improved point for the even H4 amplitude.  The N130/N170 protocol is a proposed search, not a completed root.

Even without a zero, several `lambda` points can test whether the post-H4 scalar coefficient correlates approximately with `[A_T4^+(lambda)]^2`.

## 8. A source-sign alternative if no improved point is available

For a physical H4 source `g`, compute normalized roots/free energies at `+g` and `-g`:

```text
R_odd(g)  =[R(g)-R(-g)]/2,
R_even(g) =[R(g)+R(-g)]/2-R(0).
```

Then

```text
R_odd/g -> linear H4,
R_even/g^2 -> quadratic H4 composite
```

at small `g`, after including thermal retuning, normalizer and connected contact terms.

This is the bounded finite-width test proposed in `post-h4-linear-vs-quadratic-scalar-mechanisms-20260914.md`.

## 9. Model challenge for future radial data

After angular H4 removal, compare at minimum:

```text
M1:  root residual ~ L^-7,
     independent linear scalar x=33/4 candidate;

M2:  root residual ~ L^-29/4,
     quadratic leading-H4 composite;

M3:  L^-7 (a+b log L) or nearby logarithmic collision model;

M4:  mixture of M1 and M2;

M5:  explicit higher-harmonic leakage model using the measured projector gains.
```

The primary score should be held-out prediction / cross-geometry consistency, not which intercept is closest to a preferred `pc`.

Because `7` and `7.25` are very close, no two-size extrapolation can discriminate them honestly without a remainder model.

## 10. Relation to the older `2 omega=3` speculation

Issue #47 also mentioned a speculative mechanism using a correction-to-scaling length exponent `omega=3/2`, whose second order gives a relative correction `q=3` and hence root exponent seven.

That remains a logically distinct composite route.  The current quadratic H4 mechanism instead follows directly from the **observed leading root correction** and predicts relative correction

```text
q=omega_4=13/4,
```

hence root exponent `4+13/4=29/4`.

The two composite mechanisms should not be conflated.

## 11. Claim boundary

Exact/scaling algebra once the leading H4 exponent is accepted:

```text
omega_4=13/4,
second-order root exponent=29/4,
spin4 tensor square contains H0 and H8.
```

Hypotheses:

- the corresponding connected quadratic coefficient is nonzero for square-site Matching One;
- its finite-size amplitude is large enough to explain the observed post-H4 residual;
- a particular improved-action family can tune the microscopic H4 coupling independently of scalar channels.

The main recommendation is immediate: **add `29/4` to every post-H4 model challenge before promoting exponent seven to a new linear scalar field.**