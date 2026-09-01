# Issue #275: C3 amplitude identifiability

## Status and evidence boundary

This note began as a zero-new-sample synthesis of one exact conditional design
and one completed rho-child production block.  A subsequent prospective paired
physical-rotation gate has now supplied the previously missing second reading
for the **primitive real-C3 observer**.  The later result is recorded here as
`branch_only`; it is not silently transferred to `E_top`, original `q/E`, or
pooled-root `U`.

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
- frozen paired primitive-C3 contract:
  `origin/analysis/p275-gaussian-c3-phase-20260901` at
  `9eeb6700bb59832d21b68993dc9e041ea0ce2e76`, file
  `experiments/p275-gaussian-c3-phase-20260901/CONTRACT.json`;
- paired primitive-C3 result and interpretation:
  `0b9e89c9528ab283a8175adcd596d9b0ac5047c1`, files
  `results/p275-gaussian-c3-phase/RESULT.json`, `REPORT.md`, and current
  interpretation head `93e470669e15fcec3eaee24538fcbaa80c510d31`,
  `notes/p275-real-c3-harmonic-split-20260901.md`.

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

## Completed paired primitive-C3 phase gate

The required signed-real gain law was frozen before production.  Two equal-area
Gaussian quotients, `g1=8+i` and `g2=7+4i`, have `N=65` and physical rotation
`delta=atan2(5,12)=22.619864948 degrees`, close to the other maximal
signed-real separator at `22.5 degrees`.  The same counter-derived 130-bit
field was applied to both edge orders for 2,000,000 paired replicas in 100
batches; no top-up was allowed.

```text
mean z1 = +0.00208049662712 - 0.00262405113148 i
mean z2 = -0.00110625337288 + 0.00219494722788 i

H4: chi2=73.6412/1, p=9.36904e-18, signed gain=+0.76842
H8: chi2= 1.1122/1, p=0.291603,  signed gain=+0.73660
```

The frozen binary decision is `H8_SELECTED_H4_STOP`.  It excludes the pure
`+4 delta` signed-real transport law for this primitive observer and retains
the `-8 delta` alias.  A post-reveal signed-scalar H0 line also survives
(`chi2=1.32069/1`, `p=.250468`) because `-8 delta` is almost real for this
arithmetic pair.  Thus the gate identifies an H8/even branch against H4, not a
unique local H8 field.  It does not exclude global `A_top` H4, mixed H4/H8
transport, complex geometry gains, the rho-child `E_top` observer, or
original-U/Jordan mechanisms.

The same branch notes that a literal local weight-8 insertion should decay no
slower than `N^-3`.  Transporting that radial law into the global K1/K2
residual relative to the proposed weight-21/4 Q4 leading term fixes
`kappa=2^(-11/8)`.  The existing-data forward audit now excludes even the
four-independent-amplitude envelope of that transplant
(`chi2=13.360937/4`, nominal `p=.00964044`).  Thus the observed primitive H8
phase cannot be copied into the global residual as one unmixed canonical
weight-8 radial mode; an observer form factor or sector-mixing map is required.

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

The completed sequence supports the finite statements

```text
one angle with arbitrary complex amplitudes is nonidentifying;
under the frozen binary signed-real transport, primitive real-C3 selects H8 over H4;
post reveal, a signed-scalar H0 alias also survives;
the naive canonical weight-8 member's radial transplant into global K1/K2 is excluded.
```

The second physical rotation is no longer the missing input.  The unique next
physical input is:

> an explicit observer transport map from the finite primitive-C3 H8 form
> factor to rank-0/rank-2 restricted traces of the same original `q/E` source
> and pooled-root physical normalizer.

Without that map, the primitive H8 and global H4 results are a resolved sector
split, not two votes on one amplitude.  Additional sizes of the primitive
observer may test persistence but cannot by themselves identify original-U.
