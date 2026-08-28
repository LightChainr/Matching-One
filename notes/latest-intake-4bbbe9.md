# Intake after angular-root confirmation and cross-size covariance audit

Source server checkpoint: `4bbbe90180a9960dfae613b47548897ae1defa8b`.

This note updates the post-doubling decision layer with two results that arrived
while PR #56 was being prepared. It does not change any frozen prediction.

## 1. The angular-normalized root-amplitude test passed

Issue #45 preregistered

\[
A_p(N)=-\frac{N^2\Delta p_N^*}{\Delta\cos4\theta}
=\frac{A_M}{B}
\]

with frozen target

```text
0.4510066 +/- 0.02013.
```

A clean-source, 100-million-per-size threshold-rank run gives

```text
N=65: 0.4203381 +/- 0.0215723
N=85: 0.3949486 +/- 0.0307810
frozen joint chi-square: 2.42667 / 2
zero-effect chi-square: 461.28 / 2
```

The direct root amplitude and independently reconstructed `A_M/B` agree at
both sizes:

```text
N=65: 0.4203381 versus 0.4203053
N=85: 0.3949486 versus 0.3949832.
```

Thus the previously noted drift of bare `N^2 DeltaRoot` was caused by omitting
the changing angular leverage. It was not evidence against the residual-to-root
mechanism. The correctly normalized frozen primary test is complete.

This adds a second prospective mechanism closure after Gaussian doubling:

1. fixed-p matching residual follows the angular/radial law;
2. the predicted angular-normalized root amplitude is observed independently.

Secondary sizes and the parameter-free full-curve ratio `-1/4` remain open.

## 2. Cross-size covariance is real but not decisive

PR #46 reconstructs the existing implicitly coupled P33 batches with aligned
batch covariance and delete-one root pseudo-values. The largest observed
absolute cross-size correlations are about `0.22` and have mixed signs.

The held-out constant-amplitude scores are

```text
A_M: full covariance 5.5300 / 2; diagonal 5.2733 / 2
A_p: full covariance 5.5516 / 2; diagonal 5.2929 / 2
```

The low-stat root-doubling score is

```text
full covariance 3.4625 / 2; diagonal 3.4417 / 2.
```

Therefore accidental cross-size coupling did not create or remove the old P33
radial tension. The low-stat data were simply noisy. The audit is scientifically
useful, but the generic tool still requires:

- one common aligned-batch weight assertion or explicit weights;
- covariance eigenstructure and stable factorization/pseudoinverse;
- finite-batch/bootstrap calibration;
- an exact synthetic full-versus-diagonal regression.

## 3. Frozen RNG policy

Use shared random fields inside one same-N orientation pair. Across distinct
`N`, domain-separate by default. Deliberate parent/child coupling is allowed
only when a prespecified pilot demonstrates lower variance for the exact
lineage residual and the run retains aligned batches and full covariance.

This avoids accidental dependence in ordinary radial fits while preserving
common-random-number gains where they are actually measured.

## 4. Engine readiness

Server commit `4bbbe9` adds the prospective primitive designs

```text
N=185: (13,4)/(11,8)
N=265: (16,3)/(12,11)
N=290: (13,11)/(17,1)
```

and tests the frozen order. The first two are the prospective held-out H4
sizes; the third is the preregistered third `1+i` doubling descendant.

The norm-5 H4-versus-H12 children `N=325,425` are deliberately not yet added to
the production engine. Their run remains gated on a billion-scale power table,
paired-control evaluation, and the clean RNG policy.

## 5. Updated immediate order

1. harden PR #46 and finish the clean/RNG-policy contract in #39;
2. run the full-curve doubling triptych #49;
3. run prospective `N=185,265` held-out tests #43;
4. run the third `N=145 -> 290` doubling lineage #50;
5. complete paired exact motif controls #40;
6. power and execute the norm-5 H4-versus-H12 test #57;
7. measure the linked derivative-parity spectrum #48;
8. perform self-dual/self-matching parity controls #42/#44;
9. pursue LCFT operator identification only after these empirical gates.
