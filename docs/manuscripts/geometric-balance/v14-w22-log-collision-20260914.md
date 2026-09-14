# The post-H4 scalar at x=33/4 sits on a generic-Q logarithmic collision

Date: 2026-09-14

Status: exact Kac-dimension collision + literature-backed representation-theory warning + research synthesis.  The collision does **not** by itself identify the measured H4-projected root residual with a specific logarithmic partner, and it does not prove the matching parity of the `V_<1,4>` branch.

## 1. Why this matters now

The deterministic same-angle H4-null safe-root projector gives, for the same axis/(3,4) angular pair,

```text
ell=5:  (p_H0-pc_ref) ell^7 = -0.7660485,
ell=10: (p_H0-pc_ref) ell^7 = -0.7572352.
```

The historical post-H4 mechanism note had already identified the standard Potts scalar

```text
V_<1,4>:
h=hbar=33/8,
x=33/4,
spin=0,
```

as the natural `q=3` channel producing a root correction `ell^-7`, conditional on matching-odd parity and a nonzero lattice coupling.

The new deterministic result therefore makes the operator identity question concrete.  However, at `Q=1` the conformal weight `33/8` is not isolated.

## 2. Critical-Potts Kac parametrization

Use the standard generic critical-Potts Coulomb-gas parameter

```text
Q = 4 cos^2(pi beta^2),
```

with percolation at

```text
beta^2 = 2/3.
```

A convenient Kac-weight convention is

```text
h_{r,s}(beta)
 = (c-1)/24
   + (1/4) (r beta - s/beta)^2.
```

The repository's critical-branch notation `V_<1,s>` is the branch-swapped convention relative to papers that call the energy Kac field `Phi_{s,1}`.  Numerical weights are convention-independent.

For the candidate scalar,

```text
h_{1,4}
 = (c-1)/24
   + (1/4)(beta-4/beta)^2.
```

Now consider the generic-Q non-diagonal module `W(2,2)`.  The known Potts logarithmic representation has top/bottom diagonal fields with weight

```text
h_{2,-2}
 = (c-1)/24
   + (1/4)(2 beta+2/beta)^2.
```

At `beta^2=2/3`, both equal

```text
boxed:
h_{1,4}=h_{2,-2}=33/8.
```

Thus the scalar `x=33/4` Kac branch hits the diagonal eigenvalue of the `W(2,2)` logarithmic diamond precisely at percolation.

## 3. The collision is isolated at the physical Potts point

Put

```text
t = beta^2.
```

The common `(c-1)/24` piece cancels, and the exact chiral weight difference is

```text
Delta h(t)
 := h_{1,4}-h_{2,-2}
 = (1/4)(-3 t -16 +12/t).
```

The zero condition is

```text
3 t^2 +16 t -12 =0,
```

with roots

```text
t=2/3,
t=-6.
```

Hence the only physical positive collision is the percolation point `t=2/3`.

The scaling-dimension splitting is

```text
Delta x(t)=2 Delta h(t).
```

At percolation,

```text
d_t Delta x |_(2/3) = -15.
```

Meanwhile

```text
dQ/dt
 = -4 pi sin(2 pi t),
```

so at `t=2/3`,

```text
boxed:
d_Q Delta x |_(Q=1)
 = -15/(2 pi sqrt(3)).
```

Numerically this is about `-1.3783` per unit `Q`.

This is a sharp generic-Q splitting velocity for any future branch-resolved test.

## 4. Generic-Q representation input

Grans-Samuelsson--Liu--He--Jacobsen--Saleur, JHEP 10 (2020) 109, establish for generic Potts `Q` that non-diagonal modules with weights

```text
(h_{r,s},h_{r,-s}),
(h_{r,-s},h_{r,s})
```

belong to indecomposable diamond representations whose top and bottom fields have

```text
(h_{r,-s},h_{r,-s})
```

and form a rank-two Jordan cell of `L0` and `Lbar0`.

Therefore `W(2,2)` already carries a logarithmic pair at the diagonal weight `h_{2,-2}` for generic Q.  The exact equality in Section 2 says that the simple scalar Kac branch reaches that same diagonal eigenvalue at `Q=1`.

The detailed special-Q indecomposable extension is **not** fixed by this dimension equality alone.  It may enlarge/reorganize the generic-Q diamond, or the measured lattice observable may project mostly onto one bottom component.

## 5. Relation to modern c=0 energy-operator results

In the convention used by recent percolation LCFT work, the repository field `V_<1,4>` is the branch-swapped version of the **third energy Kac operator** often denoted `Phi_{4,1}`.

Recent work on logarithmic Kac operators at `c=0` argues that higher energy operators can become zero-norm bottom fields of higher-rank logarithmic multiplets.  In particular, the 2025 analysis by Yifei He discusses a potential rank-four Jordan structure whose bottom field is the third energy operator `Phi_{4,1} ~ epsilon''`, while emphasizing that higher-rank structures require observables beyond the simplest spin correlator to expose them.

This is strikingly aligned with the present situation:

```text
raw leading matching-odd observable:
    dominated by the ordinary thermal-Q4 tangent;

post-H4 scalar projector:
    first place where the x=33/4 energy-type bottom field becomes visible;

logarithmic partners:
    may appear only in a further normal / source-resolved residual.
```

The existence of a high-rank multiplet is therefore a plausible module-level explanation, not a reason to replace the observed pure-power bottom-field signal by an arbitrary log fit.

## 6. Why the clean ell^-7 signal does not rule out logarithmic structure

A logarithmic multiplet does not force every microscopic observable to show a large `log ell` coefficient.  A lattice source can project predominantly onto the bottom Kac field, with partner coefficients small or removed by a quotient/normalization.

This already happened one level earlier in the current research programme:

- the leading spin-four `x=21/4` channel has a clear ordinary thermal-Q4 tangent contribution;
- generic-Q/percolation representation theory still allows logarithmic energy--hull structure;
- historical pure-power H4 data show that any large log admixture is not required in that observable.

The post-H4 scalar may follow the same pattern recursively:

```text
bottom V_<1,4> pure-power piece dominates the H0 root residual,
while logarithmic partners survive only in a smaller normal/source-resolved component.
```

This is now a falsifiable working hypothesis.

## 7. Matching parity is empirically separated from continuum naming

The H4-null quantity is built from the primal/matching **charge-root difference**, so any nonzero H0 residual measured there is already an empirical matching-odd scalar block of the lattice observable.

Thus one need not prove the interchiral parity of `V_<1,4>` before establishing that the **measured block** is odd.

The remaining naming problem is:

> Does the empirical odd, spin-zero, `x=33/4` block correspond primarily to the `V_<1,4>` bottom Kac branch, to a special-Q logarithmic combination involving the `W(2,2)` diagonal Jordan pair, or to another degenerate scalar sector?

Generic-Q branch splitting is the cleanest way to answer this.

## 8. A targeted generic-Q prediction

Let

```text
A_Kac(Q),
A_W(Q)
```

be the source amplitudes of the scalar Kac branch and the `W(2,2)` diagonal logarithmic branch in a declared generic-Q lift of the same observable.

Their dimensions split linearly as

```text
x_Kac(Q)-x_W(Q)
 = -15/(2 pi sqrt(3)) (Q-1)
   + O((Q-1)^2).
```

Therefore:

### regular-bottom scenario

If `A_Kac(Q)` is regular and dominates,

```text
post-H4 scalar ~ ell^-7
```

with no compulsory leading logarithm.

### collision/log scenario

If the physical `Q=1` observable requires pole-cancelled mixing of the two branches, schematically

```text
A_Kac(Q) ~ +a/(Q-1),
A_W(Q)   ~ -a/(Q-1),
```

then the collision produces an `ell^-7 log ell` contribution whose coefficient is proportional to

```text
a * 15/(2 pi sqrt(3)).
```

The exact splitting velocity above fixes the conversion between a generic-Q pole residue and the logarithmic size coefficient.

This is much more informative than a free `A+B log ell` fit.

## 9. N1105 and generic-Q now have complementary jobs

### N1105 same-circle safe-root tomography

At fixed `Q=1` it should determine whether the post-H4 residual is primarily

```text
H0 scalar,
locked H0+H8 nonlinear mixing,
or higher D4 harmonic.
```

If H0 scalar wins, the x=33/4 branch becomes the leading candidate.

### generic-Q branch test

Only after H0 is established does it become worth spending resources to split

```text
V_<1,4>
vs
W(2,2) diagonal/log branch.
```

This ordering avoids using generic-Q machinery to solve an operator identity before the lattice has even established the angular scalar channel.

## 10. Current working hierarchy

For the post-H4 scalar residual:

```text
1. empirical odd H0 block with root exponent near 7: strong deterministic hint;
2. V_<1,4> / third-energy Kac bottom field: leading continuum candidate;
3. W(2,2)-related logarithmic collision at the same x=33/4: serious module-level ambiguity;
4. T4 x I4 H0+H8 ell^-6: exact allowed competitor, numerically small in the current two-scale same-projector decomposition;
5. higher H8/H12/scalar channels.
```

The important conceptual change is that `V_<1,4>` and logarithmic structure are **not mutually exclusive**.  The Kac field can be the bottom component of the very logarithmic structure that makes the full `Q=1` operator algebra non-semisimple.

## 11. Claim boundary

- The equality `h_{1,4}=h_{2,-2}=33/8` at `beta^2=2/3`, the isolation of the positive collision, and the splitting velocity are exact algebra.
- The generic-Q `W(2,2)` diamond/Jordan structure is a literature result.
- The precise special-Q extension when the scalar Kac branch collides with that diagonal Jordan eigenvalue is not derived here.
- The potential higher-rank third-energy logarithmic structure is a literature-motivated c=0 hypothesis/programme, not a proven identification of the Matching-One residual.
- The deterministic ell^-7 evidence uses only two same-projector scales and remains a strong diagnostic rather than an asymptotic theorem.
