# Issue #43 DeltaS channel-sign audit

## Finding

The reported literal frozen `DeltaS` failure is retained unchanged, but it
does **not** show a physical sign reversal.  The frozen positive amplitude was
estimated from the P31 `either` wrapping channel, while the N=185/265
threshold-rank engine records only rank-2 `cross` wrapping.

For every complementary primal/matching configuration on the torus,

```text
R_G_cross  + R_hat_either = 1
R_G_either + R_hat_cross  = 1
```

Therefore `S_cross=1-S_either` for each orientation and

```text
DeltaS_cross = -DeltaS_either
```

for a common first-minus-second order.  This is a deterministic observable
translation, not a fitted sign choice.

## Code and artifact trace

- The frozen artifact explicitly names `channel: either` and stores positive
  `DeltaS` predictions.
- Both production metadata files explicitly name `rank-2 cross wrapping`.
- The primary scorer reconstructs `S=(R_G+R_hat)/2` from those cross-channel
  thresholds and keeps first-minus-second ordering.
- P31's five source sizes satisfy the cross/either sign identity down to the
  stored batch-integer precision.  Their inverse-variance amplitudes are
  `+0.010603216462677735 +/- 0.000936687018246324` for `either` and
  `-0.010603216462677737 +/- 0.000936687018246330` for `cross`.
- The P48 derivative analyzers also reconstruct the cross channel and use the
  same first-minus-second and `DeltaCos4=cos4(first)-cos4(second)` convention.

Thus the orientation convention, angular factor, threshold reconstruction,
and P48 definitions are mutually consistent.  The defect is the frozen
prediction-generator/protocol bridge: it transferred an `either`-channel even
amplitude to a `cross`-channel target without the exact minus sign.

## Deterministic post-reveal transport

Applying only the exact channel identity, with the frozen amplitude magnitude,
exponent, covariance, sizes, and target observations unchanged, gives

| N | observed cross DeltaS | transported cross prediction | residual z |
|---:|---:|---:|---:|
| 185 | -6.0815376e-5 | -6.7521637e-5 | +0.667 |
| 265 | -7.0249508e-5 | -6.8919447e-5 | -0.119 |

The two-size transported statistic is `0.57003 / 2 df`, compared with the
already reported literal-positive score `240.24721 / 2 df` and zero score
`112.53891 / 2 df`.

This transported statistic is a post-reveal deterministic audit, not a
replacement preregistered endpoint.  It uses zero target-fit parameters.  The
original prospective failure remains part of the immutable chronology; its
physical interpretation must be narrowed from “even-sector sign reversal” to
“frozen source/target wrapping-channel contract mismatch.”
