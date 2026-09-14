# Cylinder Ward identity behind the thermal-tangent spin-four correction

Date: 2026-09-14

Status: **exact ordinary-CFT/Virasoro identity conditional on the ordinary `c=0,h=5/8` thermal Kac quotient**.  This note explains an important part of the finite safe-transfer `TANGENT_SPIN4` observation at the critical point.  It does not prove that the full `Q=1` logarithmic module reduces to the ordinary quotient, and it does not prove exact tangency at finite near-critical scaling variable `X`.

## 1. Question

The same-model oblique safe transfer now shows two finite facts very sharply:

1. the critical NN/matching safe-sector mismatch transforms as spin four and produces the orientation-dependent `ell^-4` charge-root shift;
2. after each orientation is recentered at its own root and normalized by its own thermal slope, the charge curves at equal physical circumference agree to `O(1e-5--1e-4)` on the tested local interval.

This suggests that the leading spin-four perturbation acts primarily by shifting the thermal scaling coordinate.

The ordinary thermal-family candidate has exactly the structure needed to test this analytically.  Let

```text
phi_t:  h=hbar=5/8,
Q4=40 L_-2^2 - 60 L_-3 L_-1 - 9 L_-4.
```

The question is whether, on a cylinder energy level, insertion of `Q4 phi_t` is proportional to insertion of `phi_t` itself.

## 2. Cylinder setup

Use a dimensionless cylinder coordinate `w` with imaginary period `2 pi`, and the plane map

```text
z=e^w.
```

Let `|m>` be any primary cylinder state of chiral weight `h_m=m`.  Insert a chiral primary `phi_h` at `w=0`.

On the plane,

```text
<m| phi_h(z) |m> = C z^-h.
```

The cylinder primary is

```text
phi_h^cyl(w)=z^h phi_h(z),
```

so

```text
<m| phi_h^cyl(w) |m> = C,
```

independent of `w`.  Thus every outer local `L_-1` insertion has zero one-point matrix element on the cylinder.

## 3. Exact one-stress-tensor Ward function

The plane Ward identity with external primaries at `0` and `infinity` gives

```text
<T(z) phi_h(1)>_m / <phi_h(1)>_m
 = m/z^2 + h/[z (z-1)^2].
```

The cylinder stress tensor is

```text
T_cyl(w)=z^2 T(z)-c/24.
```

At `c=0`, therefore,

```text
<T_cyl(w) phi_h^cyl(0)>_m / <phi_h^cyl(0)>_m
 = m + h e^w/(e^w-1)^2
 = m + h/[4 sinh^2(w/2)].
```

The exact local expansion is

```text
1/[4 sinh^2(w/2)]
 = 1/w^2 - 1/12 + w^2/240 - w^4/6048 + ... .
```

Comparing with the local Virasoro OPE

```text
T_cyl(w) phi(0)
 = h w^-2 phi
   + w^-1 L_-1 phi
   + L_-2 phi
   + w L_-3 phi
   + w^2 L_-4 phi
   + ...
```

gives the exact matrix-element ratios

```text
<L_-1 phi>_m / <phi>_m = 0,
<L_-2 phi>_m / <phi>_m = m-h/12,
<L_-3 phi>_m / <phi>_m = 0,
<L_-4 phi>_m / <phi>_m = h/240.
```

The crucial observation is that the `L_-4` ratio is **independent of the external cylinder state `m`**.

## 4. Reduce Q4 using the thermal null vector

Assume the ordinary thermal Kac quotient

```text
(L_-2 - 2/3 L_-1^2)|h> = 0,
h=5/8.
```

Cylinder translation invariance gives

```text
<L_-1 Psi>_m=0
```

for every local descendant `Psi` in this one-point matrix element.

Exactly as in the torus Ward reduction,

```text
<L_-3 L_-1 phi>_m = -2 <L_-4 phi>_m,
<L_-2^2 phi>_m    =  4/3 <L_-4 phi>_m.
```

Hence

```text
<Q4 phi>_m
 = (493/3) <L_-4 phi>_m.
```

Using Section 3,

```text
<Q4 phi>_m / <phi>_m
 = (493/3)(h/240).
```

For `h=5/8`,

```text
boxed:
<Q4 phi_t>_m / <phi_t>_m = 493/1152.
```

This coefficient is independent of the external primary state.

For the real bulk spin-four combination

```text
O4 = Q4 phi_t x phibar_t + phi_t x Qbar4 phibar_t,
```

both chiralities contribute the same number, so in this normalization

```text
<O4>_m / <phi_t phibar_t>_m = 493/576
```

before restoring the physical circumference factor and the microscopic spin-four coupling.

For a physical circumference `ell`, the level-four descendant contributes the expected extra factor proportional to `ell^-4`; a rotated lattice coupling supplies the real `cos(4 theta)` factor.

## 5. Why this is exactly a thermal-tangent statement at criticality

Let the critical Hamiltonian be perturbed thermally by

```text
H -> H + t int phi_t phibar_t.
```

First-order conformal perturbation theory gives the derivative of an excitation energy with respect to `t` from the matrix element of the thermal primary (with the common bulk/vacuum analytic piece subtracted in the excitation gap).

Now perturb instead by the ordinary thermal spin-four descendant

```text
u4 cos(4theta) int O4.
```

Section 4 says that, on every translation-invariant cylinder primary level, its first-order matrix element is a universal constant times the thermal-primary matrix element.  Therefore the singular excitation-energy correction at `t=0` obeys

```text
delta_u4 E_m(0)
 proportional to
cos(4theta) ell^-4 * partial_t E_m(0).
```

Equivalently, to first order at the critical point,

```text
E_m(t,u4)
 = E_m(t + c_mic u4 cos(4theta) ell^-4, 0)
   + higher-order / other-field terms,
```

where the nonuniversal microscopic coefficient is carried by `c_mic`, while the ordinary-CFT descendant/primary ratio is fixed by the Ward identity.

Thus:

> **For the ordinary thermal Kac Q4 candidate, critical tangency is not an accidental numerical property.  It is enforced by the cylinder Ward identity.**

This is precisely the mechanism observed in the safe-transfer root: a large orientation-dependent root translation with an almost orientation-independent leading thermal denominator.

## 6. What the identity does and does not prove about the full scaling curve

The Ward identity is a statement at the critical CFT point, i.e. at `X=0` in the near-critical scaling variable.

It proves the leading local relation

```text
G4(0) proportional to F'(0)
```

for the ordinary Q4 branch.

It does **not** by itself prove

```text
G4(X)=const * F'(X)
```

for finite `X`.  The equal-circumference safe-transfer observation that root/slope-normalized curves nearly coincide on `|y|<=0.5` is therefore genuinely additional information about the massive / near-critical continuation.

A useful decomposition remains

```text
G4(X)
 = [G4(0)/F'(0)] F'(X) + G4_perp(X),
G4_perp(0)=0.
```

The cylinder Ward identity forces the critical intercept into the tangent piece for the ordinary Kac descendant.  The numerical task measures how small `G4_perp` remains away from zero.

## 7. Logarithmic-module interpretation becomes sharper

At `Q=1`, the thermal/energy field can belong to an energy--two-hull logarithmic multiplet.  The ordinary Kac calculation above should then be interpreted as the bottom-field contribution.

A logarithmic partner or another spin-four module can add a matrix element not constrained to the same `493/1152` ratio.  Therefore a conceptually clean decomposition is

```text
spin4 response
 = ordinary-Kac tangent piece
 + logarithmic/other normal piece.
```

This changes the module question substantially.  The existence of a large root shift is **not** evidence for a large logarithmic component; the ordinary Q4 branch already predicts a thermal-tangent shift exactly.

The best places to look for Jordan/module information are instead:

1. the residual root/slope-normalized orientation dependence `G4_perp`;
2. a generic-Q branch splitting where energy and two-hull fields separate;
3. logarithmic size dependence of the H4 amplitude after ordinary power corrections are controlled;
4. observables whose moving-root/thermal-tangent projection has already been removed, such as the original-U contract.

## 8. Relation to the existing torus Ward identity

The earlier torus calculation on PR #151 found, in the same ordinary Kac quotient,

```text
<Q4 phi>/<phi> = (493/96) g2(tau).
```

The present cylinder calculation is the infinite-cylinder/local-energy analogue.  Both results have the same algebraic origin:

```text
Q4 -> (493/3) L_-4
```

under one-point translation invariance plus the level-two null relation.

The torus shape is then supplied by the torus Ward value of `<L_-4 phi>`, while the cylinder energy-level value is supplied by the local stress-tensor expansion above.

## 9. Immediate research consequence

The current evidence hierarchy should be reorganized as follows.

### Already strongly supported

```text
spin character: 4;
radial dimension: x=21/4;
critical response direction: thermal tangent, if ordinary Kac Q4 is the branch;
```

with the first two supported independently by historical prospective Gaussian tests and current deterministic safe-transfer angular controls.

### Still open

```text
amount of logarithmic-partner admixture;
finite-X normal response G4_perp;
microscopic NN-vs-matching difference coupling u4^-;
post-H4 residual irreps (H8/scalar/H12/... ).
```

The operator-identification problem should therefore no longer be phrased as “why can a spin-four field move the root?”  The ordinary thermal Q4 Ward identity already supplies that mechanism.  The harder question is:

> **what part of the observed spin-four response cannot be absorbed into this exact ordinary-Kac thermal tangent?**

That remainder is the appropriate target for LCFT/Jordan identification.

## 10. Claim boundary

- The cylinder Ward function and `493/1152` ratio are exact given the ordinary Virasoro/Kac quotient and standard state-field/cylinder map.
- Applying the ratio to the square-site safe-transfer correction still assumes the observed H4/x=21/4 lattice direction couples to this ordinary thermal descendant.
- The exact identity holds at the critical CFT point; finite-`X` tangency remains a scaling conjecture supported by deterministic transfer data.
- The argument does not remove the energy--two-hull logarithmic multiplet from the `Q=1` theory.  It identifies a precise bottom-field contribution and thereby relocates the logarithmic question to the normal residual rather than the existence of the root shift.
