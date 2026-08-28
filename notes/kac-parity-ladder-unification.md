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

is not an exact consequence of the finite matching identity. It is conditional on
the stronger hypothesis isolated in Issue #61: the continuum matching exchange acts
as an OPE/interchiral automorphism, takes the thermal field `V_<1,2>` to minus itself,
and does not rotate a degenerate/logarithmic block in a way that invalidates the
naive parity recursion.

Under that conditional parity map the first fields are

```text
identity     V_<1,1>: x=0      matching even
thermal      V_<1,2>: x=5/4    matching odd
scalar       V_<1,3>: x=4      matching even
scalar       V_<1,4>: x=33/4   matching odd
```

## Relative q=2 corrections

An `x=4` scalar gives a dimensionless torus correction

```text
L^(2-x)=L^-2=N^-1.
```

This is exactly the relative correction structure resolved in center-slope quantities
(P49 / PR #83), and it is the first preregistered correction that rescues the fresh
100M `P4[S']` replication in PR #73.

The empirical slope model

```text
Mbar'(N,cbar4)=Binf*N^(3/8) * [1 + (a+b*cbar4)/N + ...]
```

therefore has a natural decomposition to test:

- scalar `a/N`: candidate `V_<1,3>` H0 correction;
- angular `b*cbar4/N`: candidate identity-family spin-4 `x=4` correction.

This is an interpretation, not yet an identification.

There is an additional harmonic selection clue. To first order,

```text
H4 * H0 -> H4,
H4 * H4 -> H0 + H8.
```

Therefore a relative q=2 correction observed inside the H4 projector `P4[S']` is
naturally compatible with an even scalar correction such as `V_<1,3>`. A single
insertion of the even identity-family H4 field does not return H4.

## Historical post-annihilator exponent near 7

The proposed leading matching-odd thermal spin-4 field has

```text
x_T4 = 21/4,
M(pc) ~ L^-13/4,
M'(pc) ~ L^3/4,
root bias ~ L^-4.
```

For `V_<1,4>`,

```text
x=33/4,
M(pc) ~ L^-25/4 = L^-13/4 * L^-3.
```

Thus it is relative `q=3`, and dividing by `M'~L^3/4` gives a root contribution

```text
L^-7.
```

So `V_<1,4>` is a concrete standard-spectrum H0 mechanism for the historical
Mertens-Ziff post-annihilator exponent near 7, provided its matching parity is odd.

## The important ordering problem

Existence of q=3 does not imply q=3 is the first post-annihilator correction.
If the even `x=4` scalar coupling is nonzero, the mixed term

```text
T4(H4, odd) * V13(H0, even)
```

is matching odd, has H4 character, and is relative q=2. It would produce

```text
w_ann = 4 + 2 = 6,
```

which precedes the V14 q=3 / w=7 scalar.

Fresh PR #73 demonstrates that a q=2 correction is genuinely present in the `P4[S']`
channel. Therefore a modern `w≈7` result in Issue #47 would need an additional empirical
fact: the q=2 central-D/post-annihilator coefficient is zero, cancels, or is anomalously
small in that observable. It cannot be assumed away merely because a q=3 field exists.

## Structural interpretation

```text
w_ann ~ 6:
  first investigate H4 q=2 mixing, especially T4 * even-scalar(x=4)

w_ann ~ 7:
  first concrete standard-spectrum candidate is odd H0 V_<1,4>, x=33/4

w_ann ~ 8:
  nonlinear T4 * I4^2 gives q=4 and H4/H12 support

w_ann ~ 10:
  next ordinary thermal-family spin-4 quasiprimary, relative q=6
```

The angular sectors supply an independent discriminator:

- q=2 `T4*V13`: H4;
- q=3 `V14`: H0;
- q=4 `T4*I4^2`: H4/H12.

A same-N harmonic decomposition combined with leading-term annihilation is therefore
more informative than a radial exponent alone. The N=1105 four-orientation design on
`main` can separate H0/H4/H8/H12 within its declared truncation. The cheaper first gate
remains Issue #47: determine whether the modern post-annihilator exponent is closer to
6, 7, 8, or 10 before paying for the full angular decomposition.

`scripts/kac_parity_ladder.py` checks the exact Kac dimensions and exponent arithmetic and
labels the matching-parity assignment explicitly as a hypothesis.
