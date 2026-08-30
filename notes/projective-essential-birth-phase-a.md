# Projective essential-birth spectroscopy: exact Phase A

Status: exact finite-quotient certificate for Issue 334.  This note identifies
what the projective line adds, and what it does not add, before any production
stream is changed.

## Exact state

For a site-addition path let `K1` and `K2` be the first steps at which the
ambient integral homology rank is at least one and two.  When the first jump is
`0 -> 1`, the unrestricted saturation theorem gives the state

```text
H(K1) = Z ell_1,
```

with `ell_1` a primitive vector up to sign.  Thus the exact archive coordinate
is

```text
(K1,K2,ell_1,site_1,site_2),  iota=1.
```

There is no free Smith/saturation variable.  A direct `0 -> 2` jump is a real
separate case: it has no canonical projective line and is recorded as the
typed atom `DIRECT_RANK2`.

The implementation uses a dynamic program on the Boolean subset lattice.  It
counts all `N!` orderings exactly while avoiding explicit enumeration of the
3,628,800 paths of the C4 control.

## Controls and result

| quotient | N | exact paths | direct `0->2` | line support | conclusion |
|---|---:|---:|---:|---|---|
| axis L=2 | 4 | 24 | 8 | `(1,0),(0,1)` | line is not fixed by timing |
| Gaussian `(2,1)` | 5 | 120 | 0 | `(1,0),(0,1)` | line is not fixed by timing |
| C4 self-matching `(3,1)` | 10 | 3,628,800 | 518,400 | `(1,0),(0,1)` | line is not fixed by timing |

The `(1,1)` C4 quotient has only two sites and only direct rank-two births, so
`(3,1)` is the smallest C4 control with a nonempty projective mark.

All exact gates pass:

1. every line-bearing state has `iota=1`, and all cycle generators reduce to
   the same canonical primitive vector under `ell ~ -ell`;
2. under a period-basis change `P -> P U`, the mark transforms as
   `ell -> primitive(U^{-1} ell)` for rotation, shear and reflection controls;
3. the axis quotient realizes all eight D4 elements: rotations preserve
   `chi4`, reflections complex-conjugate it;
4. complement/Alexander reversal maps
   `(K1,K2,ell,site1,site2)` exactly to
   `(N+1-K2,N+1-K1,ell,site2,site1)`;
5. the full `K1` histogram is recovered by summing over projective lines plus
   `DIRECT_RANK2`.  Summing over `ell` alone recovers exactly the line-bearing
   part, and cannot recover direct rank-two births.

## Does `ell` add a direction independent of `K1,K2`?

Yes in the information-theoretic, but not yet in the spin-4, sense.

In every line-bearing timing cell of all three controls, both projective lines
occur.  Therefore `ell` is not a function of `(K1,K2)`.  More strongly, after
conditioning on the event that a line is born, the exact joint counts factor:

```text
P(ell,K1,K2 | line birth)
= P(ell | line birth) P(K1,K2 | line birth).
```

This factorization is enforced by the tiny quotients' quarter-turn symmetry,
not proposed as an asymptotic law.  The two supported lines are related by a
quarter turn and hence have the same `chi4`.  Thus Phase A finds a genuine
projective-direction coordinate, but **no independent spin-4 value on these
minimal controls**.  A larger quotient with more than one D4 orbit of
primitive lines is the first informative directional-bias test.

## Exact crosswalk to Issue 156

Let

```text
A4(p)=E[1{tau1 <= p < tau2} chi4(ell1)].
```

Configuration by configuration, the event in brackets is simply “the current
configuration has ambient rank one”, and `ell1` is its unique primitive
homology line.  Hence `A4` is exactly the fixed-`p` primitive-sector character
of Issue 156.  Saving `A4` alone creates no new observable.

What is new is the oriented boundary of this state on the subset lattice.  The
oracle verifies, coefficient by coefficient in the degree-`N-1` Bernstein
basis,

```text
dA4/dp = j4_birth1(p) - j4_exit2(p).
```

At the sink `j4_exit2`, the record retains the line born at `K1`.  Internal
rank-one transitions contribute zero because their projective line is stable.
Direct `0 -> 2` transitions contribute neither source nor sink to `A4`.

The economical production consequence is sharp: reuse the Issue 156 plateau
character, but preserve both birth records if the scientific aim is to decide
whether an H4 effect enters through rank-one nucleation or rank-one exit.

## Claim boundary and next discriminator

The certificate is finite and exact.  It does not assign a continuum field or
radial exponent.  It does identify the smallest nonredundant next datum:

```text
(K1,K2,ell1,site1,site2,DIRECT_RANK2 flag), with iota omitted/fixed to one.
```

Before a large stream, use one modest quotient whose primitive-line support
contains at least two inequivalent D4 orbits.  A nonzero conditional difference
of `j4_birth1` and `j4_exit2` across those orbits would be genuinely new; a
plateau-only `A4` measurement would merely repeat Issue 156.

Reproduce with:

```bash
python3 scripts/projective_essential_birth_oracle.py \
  --json results/projective-essential-birth/latest.json \
  --markdown results/projective-essential-birth/latest.md
python3 -m unittest discover -s tests \
  -p 'test_projective_essential_birth_oracle.py'
```
