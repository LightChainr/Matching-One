# Issue 253: continuous annulus spectrum and matching-sector split

## Outcome

The existing N325/N425/N365 archive now gives a genuine continuous-transfer
test rather than another rendering of the dyadic recurrence.  All four radii
enter the profile, so `R=7` probes the noninteger coordinate
`n=log2(R/2)`.  Rank one is disfavoured in the matching-odd channel over the
declared exponent window (`chi2=18.290/8`, `p=0.0192`, optimum at `x=-8`),
whereas J2, fixed-gap R2, generic R2 and C2 all close.  The plus channel does
not require rank two (`R1 chi2=8.954/8`, `p=0.346`).

The most informative new mechanism test is the nested sector-sharing
comparison.  For J2, one common radial eigenvalue gives `chi2=10.822/11`; two
parity-specific eigenvalues give `7.504/10`.  The improvement is
`Delta chi2=3.318` for one additional spectral parameter (`p=0.0685`).  The
fixed dimension-gap-one adversary gives almost the same result (`Delta
chi2=3.216`, `p=0.0729`).  Thus there is weak, unresolved evidence that
matching-even and matching-odd insertions see different radial generators.
It is not evidence for path memory, and the likelihood reference is
descriptive whenever an optimum hits a profile boundary.

## Spectrum/Jordan map

The 24-observation profile assigns independent amplitudes to each
`(matching parity, geometry)` readout and profiles only the transfer spectrum.
The continuous classes are:

- R1: one eigenvalue `lambda=2^-x`;
- J2: a repeated eigenvalue with the second basis vector
  `n lambda^n`;
- R2-gap1: two positive eigenvalues whose scaling dimensions differ by one;
- R2: two unrestricted positive eigenvalues;
- C2: a conjugate pair with modulus `2^-x` and phase `theta` per log2 step.

The full-data matching-odd C2 optimum is interior (`x=-1.549`,
`theta=1.653`, `chi2=1.502/4`, `p=0.826`).  Its point chi-square improves on
matching-odd J2 by only `1.854` while spending one additional nonlinear
parameter.  In contrast, the plus C2 optimum runs to the declared small-phase
boundary, as does the shared-sector C2 fit.  The scientifically defensible
statement is therefore: a finite rotating pair is an economical *point
realization for the odd sector*, but the archive does not identify it over a
Jordan or ordinary two-real spectrum.

Generic R2 also exposes the identifiability boundary rather than resolving
it: plus runs to the common-shift boundary with a large gap; minus runs to the
minimum-gap boundary.  J2 and the exact dimension-gap-one adversary differ by
only `0.0347` chi-square units in the joint shared fit.

## Independent-geometry holdout

The held-out spectrum-transfer score is stricter than the earlier recurrence
score:

1. fit the continuous J2 or gap-one spectrum to all four radii of the old
   N325/N425 block;
2. expose N365 only at R2/R4 to solve its two amplitudes;
3. predict N365 R7/R8 in both matching sectors;
4. obtain a full four-by-four residual covariance by delete-one jackknife in
   each independent dependency group and add the two covariance components.

For a common parity spectrum, J2 gives `chi2=2.944/4`, `p=0.567`; gap-one
gives `2.939/4`, `p=0.568`.  The data transfer, but do not separate the two
mechanisms.  The parity-specific old fit also passes (`p about 0.90`) only
because the plus-sector exponent hits `x=-8`, making its extrapolation and
uncertainty enormous (J2 plus residuals `-3.72,-13.44` with standard errors
`7.65,28.26`).  That is an instability diagnostic, not stronger support for
sector separation.  The common-spectrum holdout is the useful stable score.

## Dependency and archive boundary

`G_old` contains N325/N425, all radii and both parities; its full 16-by-16
delete-one covariance is retained.  `G_n365` contains N365, all radii and both
parities; its full 8-by-8 covariance is retained.  Only the covariance between
these disjoint counter domains is set to zero.

P43, P49, P50 and P57 full-curve/rank archives have no annulus-radius or
marked-insertion label.  The C4 marked-pivotal pilot has a different geometry,
size and radius grid.  They remain useful operator-context metadata but are
not numerical evidence rows in this score.

## Scientific layers

**Exact analysis fact.** The named GLS profiles, degrees of freedom, full
within-block covariance, independent-block sum, and N365 holdout construction
are explicit in the machine-readable result.

**Mechanism inference.** The odd readout needs more than a single mode within
the profile window; a parity-specific radial generator has weak information
gain over a shared generator.

**Exploratory conjecture.** The odd sector may be a rotating two-state plane
while the even sector lies near its zero-phase/Jordan degeneration.  A new
acquisition should target the sector split only after freezing a radius range
that breaks the current exponent/amplitude cancellation; this result does not
justify selecting C2 today.

## Reproduction

```bash
python3 scripts/analyze_p253_annulus_spectrum_sectors.py \
  --output results/p253-annulus-spectrum/latest.json
python3 -m unittest discover -s tests -p 'test_p253_annulus_spectrum_sectors.py'
```
