# P37: parameter-free Gaussian doubling test

The frozen prediction is

```text
Delta M(2N) / Delta M(N) = -2^(-13/8) = -0.3242098886627524.
```

This scores the protocol committed as `ceb7a57` before the run. Multiplication
by `1+i` doubles `N` and rotates the orientation by `pi/4`, so the spin-4
harmonic changes sign. No amplitude or exponent was fitted to these four
points.

## Fresh-seed result

| lineage | Delta M(N) | lineage Delta M(2N) | ratio | ratio SE | fixed residual | residual SE | z |
|---|---:|---:|---:|---:|---:|---:|---:|
| 65 -> 130 | +1.08608e-3 | -3.40830e-4 | -0.31382 | 0.0908 | +1.12879e-5 | 9.89e-5 | +0.114 |
| 85 -> 170 | +9.17570e-4 | -3.12850e-4 | -0.34095 | 0.1118 | -1.53647e-5 | 1.02e-4 | -0.150 |

The covariance-aware joint score of the two fixed-prediction residuals is
`chi2 = 0.03445` for 2 degrees of freedom. Both sign flips occur and both
magnitudes land within `0.15` standard errors of the parameter-free prediction.

The stored display order at `N=130` and `N=170` is the reverse of the doubled
lineage order. The analyzer explicitly flips those two differences; ignoring
this genealogy would erase the predicted spin-4 sign.

## Independence and provenance

- New seed: `2026100101`.
- Replica counters: `[2000000000, 2100000000)`.
- Sampling: 100,000,000 paired replicas per size, 100 equal batches.
- Sizes: `N=65,85,130,170`; no size was dropped.
- Fixed coordinate: `p_ref=0.592746050790`.
- Execution source commit: `80fbdd1e9a380a87a3c56dec7795ceebb0ada23e`.
- Source SHA-256: `9f0bb3c539d93fe6a968f80412bccf63edcca8be0f29919b5908fae8075c8489`.
- Hardware: Huawei ARM64, 16 CPUs; four concurrent 4-thread jobs.
- Per-size wall times: 212.2 s (`N=65`), 269.6 s (`N=85`), 408.7 s
  (`N=130`), and 519.9 s (`N=170`).

All eight raw batch/metadata checksums were verified after download. Empty
stderr files and the per-size stdout completion records are retained.

## Retrospective comparison

The earlier P31 seed was known before this test was written and is therefore
only a retrospective pilot. It gave residual z-scores `-0.688` and `+0.929`
and joint `chi2=1.468/2`. The fresh run is the actual frozen score; the pilot
is retained in `RETROSPECTIVE_PILOT.md` rather than pooled into it.

## Interpretation

This is strong independent numerical support for the joint spin-4 sign and
`N^-13/8` magnitude relation at the two exact doubling lineages. It is not by
itself a proof of the `x=21/4` operator assignment: the test is at one fixed
thermal coordinate and two lineages, and it does not exclude finely aligned
subleading or logarithmic sectors. The next operator-level check remains the
same ratio for the corrected full-curve thermal-even matching-odd projector.
