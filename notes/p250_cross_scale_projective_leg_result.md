# P250 N101 cross-scale projective-leg result

## Answer

The larger quotient resolves the propagator but rejects every frozen
one-component shape.  This is a stronger result than another d1--d3 exponent
fit: the projective-leg charged row is real and portable across scale, yet it
is **not** behaving like one pure thermal cylinder primary or one finite
correlation length.

The independent N505 80k stream passes the full `d=1..5` transfer gate with a
weakest real-row resolution of 7.247.  The three predictions were all fixed
before target collection; in particular, the parameter-free thermal
`2x=5/4` sine law was already committed at `a740bdd`, before either N101 run.

| frozen target | all N101 d2--d5 | held-out d5 |
|---|---:|---:|
| fixed thermal sine, `2x=5/4` | 90.406/16, p=2.11e-12 | 30.828/4, p=3.32e-6 |
| source-fitted sine, `2x=1.41454` | 96.813/16, p=1.37e-13 | 6.673/4, p=0.154 |
| source-fitted exponential, `m=0.71210` | 77.526/16, p=4.63e-10 | 40.908/4, p=2.81e-8 |

Thus the fixed lattice-correlation-length alternative is decisively wrong,
but the simple conformal-cylinder image is also wrong for this observable.
The source-fitted sine happens to reach the antipodal `d=5` row, yet fails the
joint intermediate-distance shape much more strongly; that isolated holdout
agreement cannot rescue it.

## New mechanism map

The raw local effective powers show systematic curvature rather than noise.
Across the four channels the `d1->d2` values are about 0.90--1.04, while the
`d2->d3` values are 1.30--1.42 and the better-resolved `d3->d4` values are
1.18--1.67.  No constant power, cylinder exponent, or lattice mass can absorb
that motion with one amplitude removed.

The natural update is a **mixed charged transfer spectrum**: the projective-leg
indicator overlaps at least two states whose relative weight changes with
distance, or it samples a genuinely two-dimensional torus function not
captured by the cylinder sine kernel.  This explains why the N65 three-point
window looked approximately power-like while the cross-scale held-out row
breaks every scalar continuation.  It also says what not to do next: adding
replicas at the same N101 will not turn this into one primary.

The complex target row additionally rejects zero phase jointly
(`39.559/16`, p=`9.03e-4`).  The individual phases are small, so this is best
recorded as a resolved deck-sensitive transfer component, not as a cubic/OPE
phase or field identity.

## Next discriminant

The next useful object is a two-state, amplitude-free matrix-pencil/Prony
test across the existing N65 and N101 pair rows, frozen with one singular-value
rank gate and a held-out distance.  A rank-two closure would turn the present
curvature into two transfer scales; failure would point toward the full torus
two-point kernel.  No new cubic observable is needed.

## Reproduction capsule

- branch: `analysis/p250-cross-scale-projective-leg-propagator-20260830`
- protocol / initial freeze / resolution refreeze:
  `56c228d` / `a740bdd` / `dc3d036`
- host/path: `Huawei-CodeBuddy-XPk2PZ`,
  `/workspace/p250-projective-leg-cross-scale-n505-80k`
- seed/counters: `25050510120261031`, `[0,80000)`
- raw batch SHA-256:
  `8b0e06f3fdc577c362e6f2404db60933d0cf489ee532c053288a03c44dc7fe5c`
- runtime: 56 seconds; 160 batches; complete complex `T/A` covariance retained
