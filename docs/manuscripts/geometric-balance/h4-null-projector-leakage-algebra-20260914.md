# Exact leakage algebra for two-angle H4-null projectors

Date: 2026-09-14

Status: exact trigonometric/projector algebra.  It sharpens the interpretation of #808 and prevents an H4-null combination from being mislabeled as a pure scalar channel before higher harmonics are controlled.

## 1. Setup

For a square-lattice orientation `theta`, define

```text
H_(4n)(theta)=cos(4n theta).
```

Let

```text
h_i=H4(theta_i).
```

The normalized two-angle H4-null projector is

```text
L[f]
 = (h1 f(theta2)-h2 f(theta1))/(h1-h2).
```

It satisfies

```text
L[1]=1,
L[H4]=0.
```

Hence it preserves an angular scalar and annihilates every term proportional to exactly the same H4 function, including any radial dressing of that H4 coefficient.

## 2. Higher-harmonic leakage is an exact polynomial in `(h1,h2)`

Using Chebyshev identities

```text
H8  = 2 H4^2 - 1,
H12 = 4 H4^3 - 3 H4,
H16 = 8 H4^4 - 8 H4^2 + 1,
```

one obtains

```text
boxed:
C8 :=L[H8]  = -(1+2 h1 h2),

C12:=L[H12] = -4 h1 h2 (h1+h2),

C16:=L[H16]
 = 1+8p-8p s^2+8p^2,
```

where

```text
p=h1 h2,
s=h1+h2.
```

Thus the leakage budget can be computed exactly from H4 values alone; no numerical angular fit is needed.

## 3. Ideal simultaneous H4/H8/H12 notch

To annihilate H8 in addition to H4 requires

```text
h1 h2=-1/2.
```

Under this condition

```text
C12=2(h1+h2).
```

Therefore simultaneous H4/H8/H12 cancellation requires

```text
h1 h2=-1/2,
h1+h2=0,
```

or equivalently

```text
boxed:
h1=+1/sqrt(2),
h2=-1/sqrt(2).
```

This is the ideal two-angle triple-notch geometry.

It cannot also null every higher harmonic: at the ideal point

```text
C16=-1.
```

So a two-angle projector can suppress the first several angular contaminants but can never be called an all-harmonic scalar projector.

## 4. N377 is close to the ideal triple notch

For #808 directions

```text
u1=(4,19),
u2=(11,16),
```

the committed H4 values are approximately

```text
h1=+0.6748868985,
h2=-0.7435428378.
```

The formulas above give

```text
C8  ~= +0.00361464,
C12 ~= -0.137808,
C16 ~= -0.98105.
```

Thus N377 is much better described as

```text
H4 exact notch
+ H8 near-notch
+ partial H12 suppression
+ essentially unsuppressed H16.
```

This is exactly the right gate if the expected scalar residual is larger than any plausible H12/H16 contribution at that circumference.  That last clause must be demonstrated or bounded; it is not supplied by the trigonometric projector alone.

## 5. Consequence for naming #808 output

The primary combination should be called

```text
p_perp4
```

or

```text
H4-null root combination.
```

Calling it `p_H0` is justified only after showing that

```text
|C8 P8 + C12 P12 + C16 P16 + ...|
```

is below the claimed scalar signal/error budget.

A match to `ell^-7` across sizes is strong evidence for an angular-scalar contribution only if the allowed higher-harmonic radial powers cannot mimic that scaling at the tested sizes.

This is particularly important because #61 has downgraded continuum matching-parity assignments; “scalar” here should first mean angular H0, not an OPE parity label.

## 6. What the projector does prove if the residual is robust

Suppose the numerical root calculation is certified and the H8/H12/... leakage budget is controlled.  Then a nonzero residual proves an angular contribution orthogonal to H4 at that `ell`.

This directly falsifies the strongest version of the old hypothesis

```text
all visible post-leading corrections are just H4 with radial dressing.
```

because `L[H4 f(ell)]=0` for every radial function `f`.

This is why #808 is conceptually more valuable than another large-width H4 amplitude measurement.

## 7. Radial field identification comes later

If repeated/projected data support

```text
P0(ell)-pc ~ C ell^-7,
```

then a scalar field with total dimension `x` satisfying

```text
x-x_t=7
```

has

```text
x=33/4,
```

which is compatible with the current `V_<1,4>` candidate.  But the sequence of logical steps must remain

```text
H4-null residual
-> angular H0 isolation
-> radial exponent
-> continuum field / logarithmic block identification.
```

Do not reverse the arrows.

## 8. Claim boundary

Exact:

- projector formulas for C8/C12/C16;
- ideal triple-notch condition;
- H4 radial dressing is always annihilated by an exact H4 projector.

Conditional/programmatic:

- N377 higher-harmonic leakage is small enough relative to the target scalar signal;
- a measured `ell^-7` H0 coefficient maps to a particular simple `V_<1,4>` rather than a logarithmic/mixed block.

The useful correction is terminological and structural: `H4-null` is exact; `H0` is a conclusion that requires a leakage budget.