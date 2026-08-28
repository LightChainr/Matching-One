# Conditional Kac-parity ladder: a unifying correction hierarchy

## What is exact and what is conditional

On the critical c=0 Potts branch,

```text
h_(r,s)=((2r-3s)^2-1)/24.
```

For the diagonal Kac fields `V_<1,s>` this gives exactly

```text
V_<1,1>: x=0
V_<1,2>: x=5/4
V_<1,3>: x=4
V_<1,4>: x=33/4.
```

The alternating **matching parity** assignment

```text
eta_s = (-1)^(s-1)
```

is not an exact consequence of the finite matching identity.  It is conditional on
the stronger hypothesis isolated in Issue #61: the continuum matching exchange acts
as an OPE/interchiral automorphism, takes the thermal field `V_<1,2>` to minus itself,
and does not rotate a degenerate/logarithmic block in a way that invalidates the
naive parity recursion.

This note keeps those two levels of claim separate.

## Why the ladder is interesting

Under the conditional parity map, the first four diagonal fields read

```text
identity     V_<1,1>: x=0      matching even
thermal      V_<1,2>: x=5/4    matching odd
scalar       V_<1,3>: x=4      matching even
scalar       V_<1,4>: x=33/4   matching odd
```

The dimensions line up strikingly with three independent finite-size observations.

### 1. Relative q=2 corrections

A dimensionless torus correction from an `x=4` scalar scales as

```text
L^(2-x)=L^-2=N^-1.
```

This is exactly the relative correction structure now resolved in center-slope
quantities (P49 / PR #83).  It is also the first preregistered correction that
rescues the fresh 100M `P4[S']` replication in PR #73.

The empirical two-sector slope model

```text
Mbar'(N,cbar4)=Binf*N^(3/8) * [1 + (a+b*cbar4)/N + ...]
```

therefore has a natural CFT interpretation to test:

- scalar `a/N`: candidate `V_<1,3>` H0 correction;
- angular `b*cbar4/N`: candidate identity-family spin-4 `x=4` correction.

This is an interpretation, not yet an identification.

There is an additional harmonic selection clue.  The leading matching-odd field is
an H4 candidate.  To first order,

```text
H4 * H0 -> H4,
H4 * H4 -> H0 + H8.
```

Therefore a relative q=2 correction that is observed **inside the H4 projector**
`P4[S']` is naturally compatible with an even scalar correction such as `V_<1,3>`.
A single insertion of the even identity-family H4 field does not return H4.

### 2. Historical post-annihilator exponent near 7

The proposed leading matching-odd H4 field has

```text
x_T4 = 21/4,
M(pc) ~ L^(2-x_T4)=L^-13/4,
M'(pc) ~ L^3/4,
root bias ~ L^-4.
```

For `V_<1,4>`,

```text
x=33/4,
M(pc) ~ L^-25/4.
```

Relative to `L^-13/4`, that is exactly an additional

```text
q=3.
```

Dividing by `M'~L^3/4` gives a root contribution

```text
L^(-25/4-3/4)=L^-7.
```

Thus `V_<1,4>` is a concrete standard-spectrum H0 mechanism for the historical
Mertens-Ziff post-annihilator exponent near 7, provided its matching parity is odd.

### 3. The important ordering problem

The existence of the q=3 scalar does **not** imply that q=3 is asymptotically the
first correction after leading-term annihilation.

If the even `x=4` scalar coupling is nonzero, the nonlinear mixed term

```text
T4(H4, odd) * V13(H0, even)
```

is matching odd, has H4 angular character, and is relative q=2.  It would produce

```text
w_ann = 4 + 2 = 6,
```

which precedes the V14 q=3 / w=7 scalar.

Fresh PR #73 demonstrates that a q=2 correction is genuinely present in the
`P4[S']` channel.  Therefore a modern observation of `w≈7` in Issue #47 would need
an additional empirical fact: the q=2 **central-D/post-annihilator coefficient** is
zero, cancels, or is anomalously small in that observable.  It cannot be assumed
away merely because a q=3 field exists.

## Structural interpretation table

```text
w_ann ~ 6:
  first investigate H4 q=2 mixing, especially T4 * even-scalar(x=4)

w_ann ~ 7:
  first concrete standard-spectrum candidate is odd H0 V_<1,4>, x=33/4

w_ann ~ 8:
  nonlinear T4 * I4^2 gives q=4 and H4/H12 support

w_ann ~ 10:
  next ordinary thermal-family spin-4 quasiprimary (relative q=6)
```

A stable 7 would therefore be more interesting, not less: it would mean the lower
allowed q=2 H4 correction is suppressed in the annihilator observable while an odd
scalar survives.

## Direct experimental discriminants

The radial exponent alone is not the only discriminator.

- q=2 `T4*V13` correction has H4 character.
- q=3 `V14` is H0 scalar.
- q=4 `T4*I4^2` necessarily carries H4/H12 support.

Thus a future post-leading analysis should retain angular information wherever
possible rather than reducing everything to a single axis sequence.

The cleanest long-term separation is a same-N harmonic decomposition combined with
a leading-term annihilator.  The N=1105 exact four-orientation design already in the
repository can separate H0/H4/H8/H12 within its declared truncation.  A cheaper
first gate remains Issue #47: determine whether the radial post-annihilator exponent
is closer to 6, 7, 8, or 10 before paying for the full angular decomposition.

## Reproducibility

`scripts/kac_parity_ladder.py` checks the exact Kac dimensions and exponent arithmetic.
It deliberately labels the matching-parity assignment as a hypothesis rather than
an exact theorem.
