# Post-compass status of the matching-odd spin-four operator question

Date: 2026-09-14

Status: evidence synthesis / research prioritization.  This note does not turn numerical support into a theorem and does not identify the full `Q=1` logarithmic module.  Its purpose is to stop treating logically different uncertainties as if they were one flat candidate list.

## 1. Four separate questions

The matching-root correction problem should now be decomposed into four questions.

```text
A. What angular irrep carries the leading observed correction?
B. What radial/scaling dimension carries that irrep?
C. How does that field project onto the thermal/root response?
D. Which c=0 Kac/logarithmic module realizes the field?
```

The current evidence has very different strength at these four levels.

## 2. A. Angular irrep: spin four is strongly supported

### Historical prospective Gaussian evidence

The frozen norm-five experiment in #57 compared the same-radial angular characters before child data were revealed.  The canonical 500M-per-child score was

```text
H4:   chi2 = 0.4163 / 2
H12:  chi2 = 35.1931 / 2
H8:   chi2 = 16.0120 / 2
zero: chi2 = 1.7764 / 2.
```

Thus H4 resolved the designed H4/H12 alias prospectively and strongly beat H12/H8.  The zero model was not rejected by that child block alone, so this block establishes angular character more strongly than nonzero amplitude.

### Current deterministic safe-transfer evidence

The same microscopic square-site NN/complementary-matching model, with only the integer cylinder direction changed, obeys

```text
critical safe-sector mismatch ~ cos(4 theta),
charge-root shift            ~ -cos(4 theta),
```

across axis, diagonal and several oblique directions.  Equal-physical-circumference H4 projectors remove most of the raw root bias without using a threshold reference.

This second evidence block supplies the strong nonzero amplitude missing from the #57 angular-only discrimination.

### Current verdict

```text
LEADING ANGULAR IRREP = H4 / spin 4: STRONG.
```

A scalar/spin-zero correction may survive as a smaller post-H4 residual, but it is strongly disfavoured as the dominant source of the observed `ell^-4` root displacement.

## 3. B. Radial dimension: x=21/4 has independent prospective support

The `x=21/4` law is not being inferred only from the root exponent.

### Gaussian doubling

Fresh prospective Gaussian doubling used the parameter-free prediction

```text
DeltaM(2N)/DeltaM(N)
 = -2^(-13/8).
```

The minus sign tests spin four and the magnitude tests the `N^-13/8`, equivalently `L^-13/4`, law.  Two fresh lineages gave a joint covariance-aware

```text
chi2 = 0.03445 / 2.
```

No exponent or amplitude was fitted to the child data.

### N=185/265 prospective radial test

The predeclared matching-odd H4 predictions gave

```text
x=21/4: chi2 = 3.04598 / 2,
x=17/4: chi2 = 30.24613 / 2,
zero:    chi2 = 29.40938 / 2.
```

The formal `x=17/4` non-diagonal competitor was therefore rejected on genuinely new geometries.  A proposed `x=14/3` V13 competitor was removed before scoring because its Kac-branch/parity construction was invalid for the declared channel.

### Evidence ledger

Across the independent primary matching-odd blocks currently registered in the canonical evidence ledger,

```text
H4_x21_over_4:
chi2 = 3.31346 over 5 dimensions.
```

### Current verdict

```text
LEADING H4 RADIAL DIMENSION x=21/4: STRONG NUMERICAL SUPPORT.
```

This still does not say which same-dimension `Q=1` module contributes.

## 4. C. Critical response direction: ordinary thermal Q4 is exactly tangent

The new cylinder Ward calculation gives a structural result unavailable from exponent fitting.

For the ordinary `c=0,h=5/8` thermal Kac quotient,

```text
Q4=40L_-2^2-60L_-3L_-1-9L_-4.
```

On a translation-invariant cylinder primary level,

```text
<L_-4 phi_t>/<phi_t> = h_t/240,
```

independently of the external primary state.  Null reduction plus translation invariance gives

```text
Q4 -> (493/3)L_-4,
```

hence

```text
<Q4 phi_t>/<phi_t> = 493/1152
```

for the chiral descendant in the circumference-`2pi` normalization.

Therefore, at the critical point, an ordinary thermal-Q4 insertion changes every cylinder excitation energy in a direction proportional to its thermal-primary matrix element.  In other words:

```text
ordinary thermal Q4 is exactly thermal-tangent at X=0.
```

The equal-circumference transfer data add the nontrivial massive statement that this tangency remains extremely accurate after root-centering and slope normalization over a finite local near-critical interval.

### Current verdict

```text
CRITICAL TANGENCY, CONDITIONAL ON ORDINARY THERMAL Q4 IDENTITY: EXACT.
FINITE-X TANGENCY: STRONG FINITE DETERMINISTIC EVIDENCE.
```

This explains the root shift without invoking a logarithmic partner.

## 5. D. Logarithmic/module identity: open, but the large-log alternative is not equally supported

At `Q=1`, the energy field can collide with the two-hull field and form a logarithmic multiplet.  Their level-four spin-four descendants share the same `x=21/4` at the collision point, so angular character and radial dimension alone cannot distinguish them.

However, the data do not currently require a large logarithmic admixture in the matching-odd H4 amplitude.

- The historical operator note records that adding a free logarithmic parameter worsened held-out prediction relative to the pure `13/8` law on the then-available range.
- The fresh Gaussian doubling test passed the exact pure-power ratio in two independent lineages with `chi2=0.03445/2`.
- The canonical multi-block `H4_x21_over_4` pure-power score remains good.

This is not a theoretical exclusion of a Jordan module.  A logarithmic partner can have a small coefficient in this particular microscopic source/observable even if it is required elsewhere in the `c=0` theory.

### Working hypothesis

The most economical current hypothesis is

```text
observed leading matching-odd spin4 response
 = ordinary thermal-Q4 bottom-field component
 + smaller logarithmic/other normal component.
```

The burden of proof should now be on a large logarithmic component, not on the existence of the ordinary tangent piece.

## 6. Where to look for the logarithmic partner now

Because the ordinary Q4 branch already explains the leading root displacement, the logarithmic question should be targeted at observables where the exact tangent piece is removed or can be separated.

High-information targets are:

1. **root/slope-normalized angular residual**

```text
G4_perp(X)
 = G4(X)-[G4(0)/F'(0)]F'(X);
```

2. **generic-Q branch splitting** between energy and two-hull families before the `Q=1` collision;
3. **logarithmic size dependence** of the H4 amplitude after ordinary H4 power corrections are projected/controlled;
4. **moving-root-normalized original-U**, where the leading thermal tangent is explicitly removed by contract;
5. a map/projector-resolved observable known to couple differently to the two generic-Q branches.

A free `A+B log L` fit to the raw root is now a low-information test because the dominant ordinary piece is already known and power corrections remain.

## 7. Relation to the post-H4 residual hierarchy

The leading H4 root tower can be removed algebraically using same-circle angular projectors.  The remaining residual then probes genuinely different information:

```text
H8 / spin8,
angle-independent scalar/log channel,
H12 or higher D4 harmonic,
normal component of the spin4 field.
```

The simple quadratic thermal-coordinate nonlinearity of the leading tangent H4 shift is too small to explain the observed post-H4 residual at `ell~8--10`.

Thus there are now two orthogonal residual programmes:

```text
spin4-normal / logarithmic residual,
post-H4 higher-irrep residual.
```

They should not be conflated.

## 8. Updated candidate ordering

For the **observed leading root correction**, the current ranking is:

```text
1. ordinary thermal-family Q4 bottom component, spin4, x=21/4;
2. same spin/dimension logarithmic admixture as a correction to (1);
3. other spin4 x=21/4 map/defect realization;
4. leading scalar x=21/4 explanation -- strongly disfavoured by same-ell angular projection.
```

For the **post-H4 residual**, the ranking is separate and currently less settled:

```text
H8 identity-family spin8 candidate,
small scalar/log channel,
H12/higher harmonic,
other normal spin4 contribution.
```

## 9. Research consequence

The main root-mechanism question has changed from

```text
which field can produce an L^-4 shift?
```

to

```text
what survives after subtracting the ordinary thermal-Q4 tangent that already explains the leading shift?
```

This is a materially smaller and more falsifiable problem.  It also aligns with the research-compass principle: do not keep solving the already-explained leading effect; use the residual to identify information the leading effective description necessarily discards.

## 10. Claim boundary

- Historical chi-square values are evidence-ledger / preregistered numerical results, not new calculations in this note.
- The cylinder Ward tangency is exact only conditional on the ordinary thermal Kac quotient and standard conformal perturbation interpretation.
- The ranking of logarithmic admixture is a research judgement based on current evidence, not an LCFT theorem.
- A small logarithmic coefficient in this observable does not imply absence of the energy--two-hull logarithmic multiplet from percolation.
