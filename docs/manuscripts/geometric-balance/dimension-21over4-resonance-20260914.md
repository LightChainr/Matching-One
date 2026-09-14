# The `x=21/4` resonance: scalar eight-arm versus spin-four thermal descendant

Date: 2026-09-14

Status: conceptual reorganization of #768 after the oblique safe-transfer results.  The equality of scaling dimensions below uses the standard percolation arm-exponent input; the lattice angular data are deterministic finite-width controls already committed on #771.  No LCFT indecomposable structure is asserted merely from degeneracy.

## 1. The two leading stories are exactly degenerate in total dimension

The current square-lattice sector-odd candidates have often been described as competing mechanisms:

1. a scalar alternating eight-arm / eight-leg field;
2. a spin `+/-4` descendant in the thermal family.

But their total scaling dimensions are exactly the same.

For critical percolation, the alternating `j`-arm exponent is

```text
x_j = (j^2-1)/12.
```

Therefore

```text
x_8 = 63/12 = 21/4.
```

The thermal field has

```text
x_t = 5/4,
```

and a level-four chiral descendant has

```text
x_t+4 = 21/4.
```

Thus the pseudo-critical exponent four cannot distinguish these mechanisms even in principle.

This is stronger than saying that two unrelated fits happen to be numerically close: the continuum dimensions are algebraically identical.

## 2. They differ in spin, not in cylinder energy

A scalar eight-arm field has

```text
(h,hbar)=(21/8,21/8),
spin=0.
```

The thermal spin-four descendant has, schematically,

```text
(h,hbar)=(37/8,5/8) or (5/8,37/8),
spin=+/-4,
```

with reflection selecting the real combination.

Both have

```text
h+hbar=21/4.
```

Hence an axis-aligned cylinder energy sees the same `ell^-17/4` per-length power from both.  The only clean discriminator is their transformation under rotation / map-sector data.

On a square lattice, spin four is itself invariant under a `pi/2` microscopic rotation:

```text
exp(i*4*pi/2)=1.
```

Therefore the square point group also does **not** distinguish scalar spin zero from continuum spin four.  Both belong to the trivial `C4` lattice representation.

The cylinder orientation relative to the microscopic lattice supplies the missing analyzer.

## 3. The correct leading angular decomposition

At the `x=21/4` scale, reflection and square symmetry allow

```text
E_odd^phys(pc;ell,theta)
 = ell^(-17/4)
   [ B0 + B4 cos(4 theta) ]
 + smaller terms,
```

where, schematically,

```text
B0 : scalar eight-arm / any spin-0 contribution at this dimension,
B4 : spin-four thermal-family contribution.
```

More generally there can be multiple map/operator contributions inside each angular channel, but this is the minimal two-coordinate decomposition relevant to the current debate.

The important consequence is:

> `eight-arm` and `spin-four` should no longer be treated as mutually exclusive exponent hypotheses.  They are two angular components of the same leading total-dimension shell unless an additional parity/map rule removes one of them.

## 4. Existing oblique data already resolve most of the leading shell

The current same-model oblique safe-transfer results give

```text
E_odd^phys * ell^(17/4) / cos(4theta) ~= 1.0
```

across axis, diagonal, `(2,1)`, `(3,1)`, `(3,2)` controls, with the expected sign changes.

A near-equal-circumference Pell pair is already available:

```text
axis w=7:       ell=7,
diagonal n=5:   ell=5 sqrt(2),
ell_diag^2-ell_axis^2=1.
```

The scaled critical mismatch gives

```text
F_axis = +1.0358025515,
F_diag = -1.0206765597.
```

Hence

```text
F_spin4-like = (F_axis-F_diag)/2 = 1.0282395556,
F_even       = (F_axis+F_diag)/2 = 0.0075629959.
```

The angular-even leakage is only about `0.74%` of the angular-odd projector.  Because the circumferences are not exactly equal, ordinary radial finite-size drift contributes to this `0.74%`; it is therefore an upper-scale diagnostic, not a clean estimate of `B0/B4`.

The root projector independently gives essentially the same `~0.7%` leakage.

Thus a scalar `x=21/4` component with amplitude comparable to the observed spin-four component is already strongly disfavored.

## 5. Exact equal circumference is now a scalar-eight-arm measurement

The proposed `(4,3), n=2` versus axis `w=10` pair has exactly equal physical circumference `ell=10`.

Let

```text
c4 = cos(4 theta_(4,3)) = -527/625.
```

At equal `ell`, define

```text
E_axis = ell^(17/4) E^phys_axis,
E_43   = ell^(17/4) E^phys_43.
```

If only scalar and spin-four components are relevant at this order,

```text
E_axis = B0 + B4,
E_43   = B0 + c4 B4.
```

Therefore one can solve directly

```text
B4 = (E_axis-E_43)/(1-c4),
B0 = (E_43-c4 E_axis)/(1-c4).
```

This is a cleaner interpretation of the equal-circumference computation than merely calling it another spin-four check: it is an actual two-component projector on the degenerate `x=21/4` shell.

The same algebra can be applied to the root response after dividing by the common leading thermal slope, with the usual caveat that subleading thermal anisotropy can contaminate the ratio.

## 6. Pell near-node directions become especially decisive

At the exact spin-four node

```text
theta=pi/8,
cos(4theta)=0.
```

A **scalar** `x=21/4` component survives there at order

```text
ell^-4 in the root shift.
```

A pure spin-four component does not.

For Pell approximants with

```text
cos(4theta)=O(ell^-2),
```

the spin-four `ell^-4` root contribution is geometrically demoted to `ell^-6`.

Therefore:

- if `B0 != 0`, a scalar eight-arm term eventually dominates the Pell sequence again as `ell^-4`;
- if `B0 = 0`, the leading node response is `ell^-6` or smaller and comes from spin-four geometric leakage / another angular sector.

This is a much sharper use of the Pell experiment than simply “seeing the next correction power”.

The existing `(5,2), n=2` point already lies close to the node and tracks the small `cos4theta` prediction rather than showing an obvious unsuppressed scalar floor.  That is suggestive but not asymptotic evidence.

## 7. Why exact degeneracy does not automatically imply a Jordan/logarithmic pair

The equality

```text
x_scalar = x_spin4 = 21/4
```

is a degeneracy of **total scaling dimension**.

Under the full emergent rotation group the two operators carry different spin, so they remain distinct representation sectors.  Equality of `x` alone does not produce a Jordan cell, logarithm, or operator mixing.

At the lattice `C4` level both appear in the same discrete irrep, so a microscopic observable may couple to both.  This explains why exponent and point-group selection were insufficient.  But any logarithmic/indecomposable mixing requires an additional LCFT/module statement and must be demonstrated separately.

This prevents a repeat of the earlier Jordan over-interpretation elsewhere in the repository.

## 8. Reframing the six-arm kill test

The abundant theta/T3 simultaneous-birth geometry is matching-even pathwise and therefore does not by itself compete in the matching-odd one-point channel.

The true ordering of the theory questions is now:

1. **Lower-dimension exclusion:** does any six-arm / lower operator carry matching-odd map parity and a nonzero one-point sector difference?
2. **Leading-shell decomposition:** if lower channels cancel, what are `B0` and `B4` inside the degenerate `x=21/4` shell?
3. **Operator identity:** if `B4` dominates, is it specifically the level-four thermal quasiprimary after the actual percolation module/null quotient, or another spin-four field of the same dimension?

Only step 3 is an operator-naming question.  Steps 1--2 are observable/matrix-element questions and can be attacked directly.

## 9. A new pivotal interpretation of the spin-four channel

There is a plausible mechanism that makes the thermal-descendant interpretation more natural than a literal eight-arm event.

The exact local transition decomposition has

```text
alpha = rank 0->1 insertion mass,
beta  = rank 0->2 insertion mass,
gamma = rank 1->2 insertion mass.
```

Pathwise matching exchanges `alpha <-> gamma` and leaves `beta` even.  Hence

```text
alpha-gamma
```

is a local matching-odd observable, whereas the direct jump-two / theta-spine six-arm channel sits primarily in the matching-even `beta` sector.

Both `alpha` and `gamma` are naturally controlled at leading order by the same four-arm/thermal pivotal mechanism.  Their common isotropic amplitude can cancel in the difference, leaving the first *anisotropic correction to the thermal pivotal amplitude*.  A spin-four level-four thermal descendant has exactly the required dimension

```text
x_4arm + 4 = 5/4 + 4 = 21/4.
```

This gives a concrete alternative to the slogan “the root is controlled by an eight-arm event”:

> the `21/4` exponent may arise because the **difference between two four-arm thermal pivotal amplitudes** first appears in the spin-four anisotropic correction.

The equality with the eight-arm exponent is then a resonance of dimensions, not evidence that eight geometrically alternating arms are the microscopic event being counted.

This conjecture can be tested by the signed `alpha-gamma` output of #769 rather than by the absolute jump-two mass `beta`.

## 10. Minimal next tests

1. **Equal-circumference projector.**  Use axis `w=10` and `(4,3),n=2` to extract `(B0,B4)` directly at the same `ell`.

2. **Pivotal parity.**  In #769, report topology-resolved contributions separately to `alpha`, `beta`, `gamma`, and especially `alpha-gamma`.  A six-arm class that is large in `beta` but cancels from `alpha-gamma` is not a kill of the spin-four hypothesis.

3. **Pell node only after step 1.**  If `B0` is bounded tightly near zero, a second `(5,2)` width becomes an efficient probe of the next angular sector.  If `B0` is nonzero, the node instead becomes the best way to measure the scalar eight-arm amplitude.

4. **Module calculation.**  Only after the observable shell is decomposed should one spend theory effort deciding whether the spin-four component is the thermal level-four quasiprimary, a logarithmic partner, or another map-resolved spin-four field.

## 11. Claim boundary

- `x_8=21/4` and `x_t+4=21/4` use standard percolation/CFT scaling inputs.
- spin zero versus spin four is exact representation bookkeeping.
- the oblique finite-transfer values and near-equal-length projector are existing deterministic controls.
- dominance of `B4`, vanishing of `B0`, the pivotal-amplitude mechanism, and the precise thermal-module identification remain conjectural.

The main purpose of this note is to replace a false binary choice by a measurable two-coordinate leading shell.