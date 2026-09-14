# Two-stage angular/radial improvement of the matching charge root

Date: 2026-09-14

Status: deterministic finite-transfer construction plus a mechanism-conditioned extrapolation.  The angular projector is algebraic.  The second-stage `ell^-7` Richardson step uses the historically pre-stated `V_<1,4>` scalar mechanism and is therefore a sharp internal consistency test, not a rigorous threshold enclosure.

## 1. Stage one: remove the observed H4 irrep exactly at fixed ell

For two orientations with the same physical circumference `ell`, define

```text
h_i=H4(theta_i)=cos(4theta_i),
p_i=p_ch(ell,theta_i).
```

The H4-annihilating scalar combination is

```text
p_H0(ell)
 = (h1 p2-h2 p1)/(h1-h2).
```

It removes every finite-size contribution proportional to `H4(theta)` at that `ell`, without using `pc` or an amplitude.

Use the same angular pair at two scales:

```text
axis:       H4=1,
(3,4):      H4=-527/625.
```

Deterministic safe roots give

```text
ell=5:
p_H0(5)=0.5927362453692130,

ell=10:
p_H0(10)=0.5927459750664773.
```

## 2. Stage two: use the historically named scalar q=3 mechanism

The historical post-H4 operator analysis, written before these deterministic values were available, identified the scalar

```text
V_<1,4>,
x=33/4,
spin=0,
```

as the first named linear scalar candidate after the H4 leading root correction.  Since the thermal field has `x_t=5/4`, its pseudo-critical root exponent is

```text
x-x_t = 7.
```

Thus the pre-stated mechanism predicts

```text
p_H0(ell)
 = pc + C7 ell^-7 + higher.
```

For `ell=5` and `10=2*5`, eliminate the `ell^-7` term exactly:

```text
pc_hat^(H4,V14)
 = [2^7 p_H0(10)-p_H0(5)]/(2^7-1)
 = 0.5927460516782668.
```

The inferred scalar amplitude is

```text
C7_hat
 = [p_H0(5)-pc_hat] 5^7
 = -0.7661178948268482.
```

Neither number uses an external threshold.

## 3. After-the-fact comparison only

The independent high-precision diagnostic reference used elsewhere in the repository is

```text
pc_ref=0.59274605079210.
```

Only after constructing the reference-free estimator,

```text
pc_hat^(H4,V14)-pc_ref
 = +8.86e-10.
```

This level of agreement is a strong mechanism-level positive control.  It is not a certified error bar: with only two radial scales, higher H0 corrections and angular leakage can cancel.

## 4. An unrelated same-circle cross-check

The independent N85 angular pair

```text
(2,9),
(6,7)
```

has the same physical circumference `sqrt(85)` but a different H4/H8 geometry.  Its H4-projected scalar estimate is

```text
pc_hat_N85=0.5927460512227809.
```

It differs from the two-stage estimator above by

```text
-4.55e-10.
```

This should not be treated as a second sub-nanounit threshold determination.  The N85 geometry can have accidental cancellation among H8/scalar residuals.  It is useful as an independent arithmetic consistency check that did not use the axis/(3,4) lineage or the external reference.

## 5. Why this is more informative than a raw-width fit

The construction separates two mechanisms before extrapolating:

```text
raw pseudo-critical root
 -> exact angular H4 projection
 -> named scalar ell^-7 elimination.
```

A fit of raw axial roots would mix

```text
linear H4 tower,
scalar odd sector,
H8/higher harmonics,
nonlinear field products,
thermal-coordinate corrections.
```

The two-stage estimator instead removes the observed dominant irrep algebraically and only then invokes one historically predicted scalar exponent.

This is a finite-size version of improvement / projection rather than a larger-data regression.

## 6. The remaining falsification burden

The striking numerical agreement does not by itself identify `V_<1,4>` as a simple isolated field.  Two residual issues remain.

### Angular content

A two-angle H4 projector still contains

```text
H0,
H8,
H12,...
```

with pair-dependent geometry coefficients.  N1105 same-circle multi-angle tomography (#807) is the correct held-out test of whether the post-H4 residual is genuinely H0-dominated.

### Module content

At `Q=1`, the `x=33/4` Kac branch collides with the diagonal Jordan eigenvalue of the generic-Q `W(2,2)` module.  Thus an H0 `ell^-7` block can be dominated by the `V_<1,4>` bottom Kac component while still belonging to a larger logarithmic multiplet.

The new `ell^-7` success therefore supports the **bottom-field scaling law** more directly than it proves a simple-module identity.

## 7. Decision logic

### N1105 finds dominant H0

Then the two-stage estimator has a genuine same-circle angular justification.  The next field-theory question is `V_<1,4>` bottom vs logarithmic `W(2,2)`-related combination, not whether a scalar block exists.

### N1105 finds locked H0+H8 or large H12

Then the excellent `ell^-7` two-scale behavior is partly accidental/mixed.  Keep the numerical estimator as an improvement trick, but downgrade the direct V14 interpretation.

### A third same-angular scale breaks the ell^-7 law

Then the present two-scale agreement is preasymptotic.  Do not preserve the scalar label by adding more fit terms without a new mechanism.

## 8. Claim boundary

- The H4 projection is exact angular algebra.
- The two safe-root inputs are deterministic finite-state transfer outputs.
- The exponent seven was a historical mechanism prediction for the same scalar projector, not selected by fitting these two deterministic values.
- The resulting `pc_hat` is not a rigorous threshold enclosure.
- The close N85 agreement is an independent consistency check but can contain accidental higher-harmonic cancellation.
