# P267 Gaussian x annulus context rectangle

Status: frozen production acquisition and heldout score complete.

## Outcome

This is the first numerical rectangle in which both contexts use the same
fixed-`p` root-toggle source and the same conditional landing-H4 readout.  The
shared-generator model is not rejected, and there is no evidence for the
context-enriched alternative.

The best shared candidate is nominally `lambda=1`, with
`chi2=3.9339275579` at effective rank 6 (chi-square survival reference
`0.6856175289`).  The best unrestricted context pair is
`lambda_Gaussian=1/2`, `lambda_annulus=1`, with `chi2=3.7177074013` at rank 6.
The frozen improvement is only

`Delta = min(shared diagonal) - min(all pairs) = 0.2162201565`.

The 100,000-draw Gaussian parametric bootstrap gives worst-case
`p=0.4002359976` over the three shared-null candidates.  This does not identify
`lambda=1`: the three candidates have weak separation at this acquisition
level.

## All frozen fixed scores

Rows are Gaussian `lambda`; columns are annulus `lambda`.

| Gaussian / annulus | 0 | 1/2 | 1 |
|---|---:|---:|---:|
| 0 | 5.1881566363 | 4.1925259426 | 3.7626197409 |
| 1/2 | 5.1432442967 | 4.1476136030 | **3.7177074013** |
| 1 | 5.3594644532 | 4.3638337596 | 3.9339275579 |

The shared diagonal scores `(lambda=0,1/2,1)` are respectively
`5.1881566363`, `4.1476136030`, and `3.9339275579`, each at effective rank 6.
Their chi-square survival references are `0.5199171964`, `0.6567079889`, and
`0.6856175289`.

Bootstrap plus-one p-values under the individual shared nulls are:

| shared null | p |
|---|---:|
| 0 | 0.3752862471 |
| 1/2 | 0.4002359976 |
| 1 | 0.3770662293 |

The preregistered composite-null report is their maximum, `0.4002359976`.

## Acquisition and gates

- Source implementation: commit `ee9190c9996031cb1c088b98668088c29cb2e86d`.
- Preregistration authorization: commit `f8cdb6b`.
- Huawei host: `Huawei-CodeBuddy-XPk2PZ`; 8 OpenMP threads.
- Production: 16 designs, 200,000 samples per design, 200 batches, radius 2.
- Frozen RNG: seed `26725360829`, counters
  `[26725300000,26725500000)`.
- ARM64 binary SHA256:
  `f273763dea4736db894f0074a125c52debe78a2eb1c6aa4ecef53481f096fdbb`.
- Wall time: `11.8819 s`; user time `94.86 s`; maximum RSS `5824 KiB`.
- The smallest per-design pivotal denominator is 5,822 events (N680 first);
  the largest is 25,820, so no fitted cell is sustained by a vanishing event
  count.
- Equal-N orientation pairs passed the per-batch common-field digest check.
- All 16 R2 landings are injective and retain both axis and diagonal C4
  direction orbits.  Their scalar/spin-4 response therefore remains rank two,
  as required by exact gate `83e98fc`.
- The tiny oracle exhausts all 1,024 primitive N10 fields and all 256
  nonprimitive Smith-(2,4) N8 fields; pivotal and H4 flags are invariant under
  the tested relabellings and period-basis change.
- Synthetic scoring recovers both the shared `(1/2,1/2)` truth and enriched
  `(1/2,1)` truth.

## Scientific interpretation

Exact: the 24-point base vector uses a complete 16x16 Gaussian delete-one
covariance and the existing complete 8x8 P253 N425 covariance.  The two blocks
are independent by disjoint seed/counter domains.  Every one of the nine
candidate scores is the GLS norm of the same six frozen residuals.

Mechanism inference: one candidate effective transfer generator remains
adequate across Gaussian cover doubling and annulus radius transfer for these
two same-semantics rows.  Allowing one generator per context gains too little
to resolve against the composite shared null.

Exploratory boundary: neither the nominal `lambda=1` minimum nor the
`(1/2,1)` unrestricted minimum is a lambda identification.  The candidate
lattice still has weak discrimination.  Absence of context gain also does not
prove that path/state memory is absent outside this local pivotal readout.

