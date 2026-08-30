# P250 projective-leg minimal state dimension

## Result: rank two is minimal

The apparent exponent drift is not a failed propagator.  It is the signature
of a second charged transfer state.

The covariance-aware block-Prony analysis uses N505 `d=1..4` for the shared
complex recurrence, keeps `d=5` held out, and applies the resulting roots to
the independent N325 geometry.  Four hand-charge channels share the roots but
have independent complex amplitudes.

| recurrence | N505 fit | N505 d5 | N325 fixed roots | N325 `1/L` roots |
|---|---:|---:|---:|---:|
| rank 1 | 76.725/22, p=5.56e-8 | 16.505/8, p=0.0357 | p=0.00461 | p=0.0218 |
| rank 2 | 18.830/12, p=0.0927 | 8.894/8, p=0.351 | p=0.298 | p=0.296 |
| rank 3 | 1.118/2, p=0.572 | 8.340/8, p=0.401 | not identifiable | not identifiable |

Rank one cannot be rescued by its marginal d5 or conformal-transport p-value:
it fails the N505 fit itself by seven orders of magnitude.  Rank two closes
the fit, the untouched d5 row, and both declared source-geometry transports.
Rank three buys almost no held-out improvement and cannot be checked on the
short N325 sequence.  Therefore the minimal identifiable dimension is two.

## Complex spectrum

The two N505 roots are

\[
\lambda_1=0.57018-0.03649i,
\qquad
\lambda_2=-0.00451-0.40174i.
\]

The first has magnitude `0.5713 +/- 0.0278` and phase
`-0.0639 +/- 0.0449`: it is the ordinary nearly real decay.  The second has a
noisier magnitude `0.4018 +/- 0.2646` but phase
`-1.5820 +/- 0.4937`, strikingly close to `-pi/2`.  Its channel amplitudes are
roughly 6--16% of the leading amplitude.  A weak quarter-turn oscillatory mode
is exactly the kind of component that produces the observed distance
curvature and the small but resolved complex pair phase without invoking a
cubic/OPE phase.

The second block-Hankel singular value is only 4.57% of the first, explaining
why short scalar fits looked almost rank one.  Its effect is nevertheless
statistically necessary because it is phase-orthogonal rather than merely a
small correction to the leading positive decay.

## Periodic-image alternative

The frozen nearest `3x3` two-dimensional image kernel fits an exponent
`alpha=1.5948`, but fails both independent checks: N505 d5 gives
`60.649/8`, p=`3.48e-10`, and N325 d3 gives `31.155/8`, p=`1.32e-4`.
Thus the curvature is not explained by the simplest periodic-image sum.

What remains open is the geometry transport of the two roots.  Current N325
precision and only one rank-two recurrence equation give essentially equal
scores for unchanged lattice roots and `1/L`-scaled complex logs.  A future
larger quotient should target that transport directly, not re-estimate the
state count.

## Scientific card

- **Mechanism changed:** scalar exponent drift becomes interference between a
  leading real mode and a weaker approximately quarter-turn charged mode.
- **Not proved:** exactly two continuum primaries, a thermal field identity,
  or a universal OPE phase.
- **Observer / sector / source / geometry:** complex projective-leg pair row;
  Z5 charges 1/2, both hands; existing independent N325 and N505 batches.
- **Dependency:** no new random block; this is a covariance-aware reuse of the
  two archived streams.
- **Next upweighting observation:** a third geometry predicts the two root
  transports and their phase without changing rank or fitting new states.
