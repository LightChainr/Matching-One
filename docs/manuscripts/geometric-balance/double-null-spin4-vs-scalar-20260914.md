# Double-null tomography for the `x=21/4` shell: orientation zero and hexagonal-modulus zero

Date: 2026-09-14

Status: scaling-design note combining the current oblique-cylinder result with the existing hexagonal/Pell modular programme.  The group-theory zeros are exact for a continuum spin-four one-point function; transfer to square-site finite-size amplitudes is a universality/scaling hypothesis with declared arithmetic leakage.

## 1. Why one null is not enough

The leading square matching-root correction currently lives in a degenerate total-dimension shell

```text
x=21/4,
```

which can contain at least

```text
spin 0 : scalar eight-arm-like contribution,
spin 4 : thermal-family anisotropic contribution.
```

An exponent-four root shift cannot separate them.

The oblique cylinder now supplies one powerful spin-four projector through `cos(4theta)`.  There is a second, logically independent projector already implicit in the hexagonal-modulus programme: the automorphism group of the equianharmonic torus.

Using both gives a **double-null test** of whether an unsuppressed scalar `x=21/4` floor survives.

## 2. Orientation null on a long cylinder

For a square-lattice spin-four correction,

```text
A4(theta) proportional to cos(4theta).
```

At

```text
theta = pi/8,
```

the leading spin-four contribution vanishes.

Primitive Pell directions satisfying

```text
a^2-2ab-b^2 = +/-1
```

approach that direction with

```text
cos(4theta)=O(|u|^-2).
```

Therefore a spin-four root term that is normally `ell^-4` is demoted to

```text
ell^-6
```

along the Pell sequence.

A scalar `x=21/4` contribution has no `cos4theta` zero and remains `ell^-4` unless a separate sector rule kills it.

## 3. Modular null at the hexagonal torus

Let

```text
tau_hex = exp(i pi/3)
```

(up to the chosen equivalent fundamental-domain convention).

The continuum torus has a 60-degree automorphism at this elliptic fixed point.  A one-point amplitude of spin `s` transforms by

```text
exp(i s pi/3).
```

For `s=4`,

```text
exp(i4pi/3) != 1.
```

Therefore the scalar torus one-point coefficient of a genuine spin-four field must vanish exactly at `tau_hex`:

```text
boxed:
F4(tau_hex)=0.
```

A spin-zero `x=21/4` field is not killed by this automorphism.

Thus the finite-size charge/root correction at the `x=21/4` shell should have the local structure

```text
V_21/4(tau,theta_lat)
 = b0 F0(tau)
 + b4 Re[e^{i4theta_lat} F4(tau)].
```

At `tau=tau_hex`, the second term vanishes while the first generically does not.

This is a different zero from the cylinder orientation node: here the microscopic lattice orientation can be fixed; the continuum torus modulus supplies the selection rule.

## 4. Arithmetic approximants predict an extra two powers

Existing Pell/Eisenstein-style integer-period approximants approach the elliptic modulus with a shape error of order

```text
delta tau = O(N^-1)
```

for area `N~L^2`, i.e.

```text
delta tau = O(L^-2).
```

If the spin-four one-point coefficient has a generic first-order zero in the local modular coordinate,

```text
F4(tau_L) = O(delta tau),
```

then its ordinary `L^-4` root correction becomes

```text
L^-4 * L^-2 = L^-6.
```

This mirrors the orientation Pell node, but the origin of the extra `L^-2` is completely different.

A scalar `x=21/4` component would remain `L^-4` after subtracting the known continuum shape baseline.

The order of the modular zero should be checked in the correct elliptic local coordinate; a higher-order zero only strengthens the suppression.

## 5. The double-null decision table

Suppose the ordinary square/axis sequence has an `L^-4` root correction.

### Both null sequences become `L^-6` or smaller

Strong evidence that the unsuppressed leading amplitude is spin four and that any scalar `x=21/4` component is small/zero in the rank-source channel.

### Orientation node suppresses but hex modulus does not

Then the long-cylinder `cos4theta` mechanism is real, but a scalar/map contribution can survive on finite-aspect tori.  The leading shell is genuinely multidimensional.

### Hex modulus suppresses but orientation node does not

Then the observed cylinder angular law has likely mixed geometry/metric effects or a different modular representation; re-audit the physical direction normalization.

### Neither suppresses

A scalar `x=21/4` floor or another angular-even mechanism is present at comparable scale; the pure spin-four interpretation is incomplete.

## 6. Why the tests are statistically and theoretically complementary

The two nulls have different nuisance directions.

### Orientation Pell

Keeps the semi-infinite cylinder logic and changes the homology direction relative to the microscopic square lattice.  Main nuisance: row-memory / finite-circumference corrections and the arithmetic approach to `theta=pi/8`.

### Hexagonal modulus

Keeps a finite-aspect torus and changes the modular shape.  Main nuisance: continuum shape subtraction, integer-period approximation, and mixing with other torus solution sectors.

A common scalar floor should survive both.  A genuine spin-four one-point coefficient is constrained by both independent symmetries.

## 7. Relation to #156 homology-character tomography

At the exact hexagonal modulus, the **scalar** spin-four one-point amplitude vanishes, but spin information can still be present in nontrivial homology/map characters.  This is precisely why #156's C3/projective-homology tomography is useful as a positive control rather than merely another scalar null.

However, the matching root itself is an aggregate rank-charge observable.  The positive character channel should not be substituted for the original root observable.  Its role is to verify that the spin-four sector exists even where the scalar projection kills it.

This yields the desirable pattern

```text
scalar rank-charge H4 response: null at tau_hex,
nontrivial homology character:  allowed positive spin signal.
```

subject to the exact lattice-to-character dictionary already being maintained in #156.

## 8. Relation to modular tomography (#585)

Any map-resolved torus basis proposed for the h-odd `x=21/4` correction must respect the elliptic selection rule:

```text
spin-four basis coefficient vanishes at tau_hex,
spin-zero basis coefficient need not.
```

Together with the large-aspect cylinder boundary, this gives two nonlocal constraints on the same torus solution:

```text
rho -> infinity : reproduce cos4theta cylinder amplitude,
tau -> tau_hex  : vanish in the scalar spin-four channel.
```

A modular basis that satisfies only one is not an adequate physical candidate.

## 9. Minimal production rule

Do not start a generic modulus ladder.

Use already-planned/available general-period infrastructure and ask only whether one carefully chosen sequence approaching `tau_hex` can distinguish

```text
L^-4 scalar floor
vs
L^-6-or-smaller spin-four leakage.
```

The continuum critical rank/homology baseline must be evaluated at the **actual** finite modulus before forming the lattice residual, as already required by #156.

Likewise, do not add a second orientation Pell width until the exact-equal-circumference `(4,3)` projector first constrains the scalar coordinate `B0`.

## 10. Claim boundary

Exact/group-theoretic:

- `cos4theta` orientation zero;
- 60-degree elliptic stabilizer kills a scalar one-point spin-four amplitude at `tau_hex`;
- spin-zero is not symmetry-killed by either spin rule.

Conditional/scaling:

- `L^-4 -> L^-6` suppression along the declared arithmetic approximants;
- application to square-site matching-root finite-size corrections;
- dominance/absence of the scalar eight-arm coordinate.

The double-null strategy is useful precisely because the two zeros arise from independent geometries while targeting the same spin representation.