# Issue #275: C3 amplitude identifiability

## Status and evidence boundary

This note is a zero-new-sample synthesis of one exact conditional design and
one completed production block.  It does not rerun either source, fit a
continuum field, or combine observables that have different physical
normalizers.

Both source results are `branch_only` relative to the current integration
branch `research/navigation-priority-refresh-20260829` at
`f6a867592e2eb4ca9bb14d787b64d5b1b051681c`:

- exact phase design:
  `origin/review/p275-c3-amplitude-contract-20260901` at
  `3f3d670d6904919ecf41860ae47ab478d0c7f606`, file
  `experiments/p275-c3-phase-contract-20260901/phase_design.py`;
- rho-child production:
  `experiment/p267-rho-c3-etop-20260830` at
  `2402a3330b421595d3573337a5723ff3dbdcb7e9`, files
  `results/server-20260830/P267-rho-C3-Etop-N112-2M/batches.csv`,
  `run.json`, `score.json`, and producer/scorer
  `scripts/p267_rho_child_etop_mc.py` and
  `scripts/score_p267_rho_child_etop.py`.

The exact design is conditional algebra.  Its own boundary is: no existing-C3
data reanalysis, no choice of a new geometry, and no assertion that an
unnormalized trace or normalized response obeys a pure-spin transport law.
The production result is a finite-`N=112` square-bond topology response, not a
continuum field identity, exponent measurement, state-count decision, or
square-site matching transfer.

## Production block and dependency

The production schemas are
`matching-one/p267-rho-child-etop-c3-run/v1` and
`matching-one/p267-rho-child-etop-c3-score/v1`.  The entire result is one
dependency group and therefore one evidence unit:

```text
seed                 2672751123001
replica interval     [27500000000, 27502000000)
samples              2,000,000
aligned batches      100
children             2omega, omega_over_2, omega_plus_1_over_2
common field         the same counter-derived 224-bit bond mask in child edge order
```

Each CSV batch retains, for every child, rank-0/rank-1/rank-2 counts and the
real and imaginary primitive-H4 sums.  `run.json` retains the mean and full
`9 x 9` covariance in the order

```text
2omega_Etop, 2omega_H4_re, 2omega_H4_im,
omega_over_2_Etop, omega_over_2_H4_re, omega_over_2_H4_im,
omega_plus_1_over_2_Etop,
omega_plus_1_over_2_H4_re, omega_plus_1_over_2_H4_im.
```

The primary nontrivial C3 coordinate in `score.json` is

```text
Etop_r1 = (-0.0017446307362800637, -0.000047198384506250046)

Cov(Etop_r1) =
[[ 2.142175049102127e-8,  -1.8622340595521373e-11],
 [ -1.8622340595521373e-11, 2.100517066498321e-8 ]]

chi2 = 142.199239404773 on 2 df
p    = 1.3238156099792117e-31.
```

Thus the finite-lattice Alexander-even rank-redistribution C3 response is
resolved.  Resolution of this one complex coordinate is not, by itself,
identification of its spin label.

The same-stream observer-ray determinant is

```text
D = H4_r0 * Etop_r1 - H4_r1 * Etop_r0
  = (-3.2971460760194775e-6, -4.797244050665041e-7)

Cov(D) =
[[ 5.259100420032095e-13, -5.574025668283198e-15],
 [ -5.574025668283198e-15, 4.987579124085081e-13 ]]

chi2 = 21.202316278511617 on 2 df
p    = 2.4887170225840352e-5.
```

This rejects a common child-character ray for the two retained observers.
It does not provide a second physical rotation of one observer, and it does
not turn primitive H4 and `E_top` into two readings with a shared amplitude
or normalizer.

## Exact one-angle nonidentifiability

A real C3 orbit has three real readings.  Its first discrete Fourier
coordinate is one complex number.  At physical angle `theta`, the two
conditional pure-spin descriptions may be written

```text
H4: z(theta) = A4 exp(i 4 theta),
H8: z(theta) = A8 exp(-i 8 theta),
```

with the conjugation convention absorbed into the unknown complex amplitude.
At one angle, each model maps its arbitrary complex amplitude onto all of
`C`.  Consequently both real design images are `R^2`, their intersection is
`R^2`, and no positive-definite covariance can separate them.  Applying this
to the measured `Etop_r1` gives a strict `NONIDENTIFIABLE` spin decision even
though the coordinate itself is very far from zero.

The observer-ray determinant does not alter this rank statement: adding a
second observer with its own unrestricted complex amplitude changes both
model images to the same product space rather than supplying a calibrated
rotation of the first observer.

## When a second angle contains information

Let the second reading be taken at `theta + delta`.

Under a shared-complex-amplitude transport contract, the two predictions are

```text
H4: z2 = exp(i 4 delta) z1,
H8: z2 = exp(-i 8 delta) z1.
```

The combined real design has joint rank four exactly when
`sin(6 delta) != 0`.  Full-rank covariance whitening changes the numerical
separation but not this algebraic rank.

Under independent but nonzero signed-real gains, the informative quantity is
relative phase.  The exact design supplies the two restrictions

```text
H4: Im[z2 conj(z1) exp(-i 4 delta)] = 0,
H8: Im[z2 conj(z1) exp( i 8 delta)] = 0.
```

They alias when `sin(12 delta) = 0`.  In particular, a `15 degree` rotation
(`delta = pi/12`) has an exact counterexample: H4 gain `+1` and H8 gain `-1`
produce the same two readings.  A `7.5 degree` rotation
(`delta = pi/24`) has maximal signed-real phase separation and distinguishes
the two predictions for nonzero real gains.

If an arbitrary complex gain is allowed independently at the second angle,
both models span all of `C^2` at every `delta`.  No choice of angle then makes
the spin label identifiable.  A prospective second reading is therefore
scientifically meaningful only after the cross-rotation gain law is frozen.

## Normalizer firewall

The following integrated archives are useful controls within their own
observable bases, but they cannot be appended to the rho-child C3 vector as
extra rows without an explicit physical source-to-readout map and the same
normalizer:

- norm-4 K1/K2 reuse:
  `analysis/norm4_two_activation_h4_manifest.yaml`,
  `results/norm4-two-activation-h4/latest.json`, and
  `scripts/analyze_two_activation_h4.py`, schema
  `matching-one.two-activation-h4.v1`.  Its coordinates are the pooled-root,
  exact-`Delta cos(4 theta)` normalized `angular_delta_F1`,
  `angular_delta_F2`, and derived `angular_delta_M`.  N65/N85/N130/N170 share
  dependency group `crn-2026104501-5100000000-7000000000`; N260 and N340 are
  the separate groups `crn-2026105401-8200000000-9200000000` and
  `crn-2026105402-8200000000-9200000000`.
- P43/P57 matching-odd reuse:
  `results/server-20260828/P43-heldout-fullcurve-500m/analysis/primary_score.json`,
  `results/server-20260829/P57-norm5-500m/primary_score.json`, and
  `results/evidence-ledger/issue212-matching-odd-synthesis.json`.  The latter
  explicitly treats `issue43_n185_n265_500m_histograms` and
  `issue57_norm5_production` as distinct raw-data groups with block-diagonal
  joint covariance.  Their raw `DeltaM` orientation contrasts are
  matching-odd rows, not the rho-child Alexander-even C3 observable.

Numerical equality of dimensions, spin vocabulary, or the existence of a
covariance matrix is not a normalizer identity.  In particular, neither the
norm-4 pooled-root normalization nor the P43/P57 raw matching contrast may be
used to restrict the complex C3 amplitude post hoc.

## Decision and unique missing input

The completed one-angle result supports the finite statement

```text
resolved nonzero Etop_r1 + distinct Etop/primitive-H4 observer rays
does not identify H4 versus H8 under arbitrary complex amplitudes.
```

The unique missing physical input for a positive H4-versus-H8 decision is:

> a phase-calibrated second physical rotation of the same real-C3 observable,
> at the recommended non-alias angle `delta = pi/24`, retaining the identical
> physical normalizer and conditioning, with the cross-rotation gain law
> frozen in advance as shared-complex or nonzero signed-real transport.

If that gain remains arbitrary complex, stop at `NONIDENTIFIABLE`; additional
archive concatenation or a different uncalibrated observer cannot repair the
rank degeneracy.
