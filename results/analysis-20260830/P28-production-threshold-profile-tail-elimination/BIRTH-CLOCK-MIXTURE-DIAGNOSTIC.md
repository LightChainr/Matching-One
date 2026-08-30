# P28 K1/K2 birth-clock mixture diagnostic

This post-reveal mechanism diagnostic uses the same held-out sizes, sides, orientations, and
synchronous delete-one covariance as the frozen composite result.  It adds no simulation and no
new tail-law vote.

## Exact reconstruction gate

The histogram convention is:

- `K1 = K_minus`: first essential ambient-H1 birth;
- `K2 = K_plus`: second essential ambient-H1 birth.

For every archived curve, the equal mixture reconstructed from the two marginal Beta
order-statistic distributions reproduces the composite CDF and density.  The largest numerical
error among CDF, density, and the log-mixture identity is below `1.6e-15`.

## Mechanism result

Separately standardized, neither birth clock has a pure `4/3` tail:

| component | post-reveal diagnostic chi-square / df |
|---|---:|
| K1 | `14,918,477.15 / 48` |
| K2 | `15,965,901.67 / 48` |
| composite | `6,280,338.87 / 48` |

The clock/side split is strongly complementary:

| clock | left marginal chi-square | right marginal chi-square | dominant side |
|---|---:|---:|---|
| K1 | `1,809,197.88` | `11,950,834.22` | right |
| K2 | `12,517,817.38` | `2,730,295.64` | left |

These marginal quadratic forms retain their clock/side covariance but do not add to the joint
score.  They identify the direction of the internal curvature, not independent evidence counts.

## Does mixing create the curvature?

At each composite-standardized point, let
`w_j=rho_j/(rho_1+rho_2)`.  Then exactly

```text
log[(rho_1+rho_2)/2]
  = w_1 log rho_1 + w_2 log rho_2
    + [-w_1 log w_1-w_2 log w_2-log 2].
```

Applying the frozen `4/3` residual projector and composite covariance gives an exactly additive
but potentially signed GLS attribution:

- responsibility-weighted component shape: `+9,110,751.10`;
- separation-entropy term: `-2,830,412.23`;
- composite: `6,280,338.87`;
- additivity error: below `8e-7`, or roughly `1e-13` relative.

The separation term cancels about 31% of the component-shape contribution.  The mixture therefore
**hides** curvature; it does not generate the rejection.

## Size drift

| channel | N=265 beta-effective range | N=425 beta-effective range |
|---|---:|---:|
| K1 left | `1.618–1.657` | `1.568–1.599` |
| K1 right | `2.195–2.238` | `2.116–2.118` |
| K2 left | `1.917–1.939` | `1.880–1.914` |
| K2 right | `1.873–1.997` | `1.767–1.814` |

All 32 clock/orientation/side constant-beta diagnostics reject at `p<0.01`.  The two `N=290` K1
right-tail cells fall below the original count gate (`6.6–6.8k`, or `66–68` per batch), so those
cells remain descriptive.  Every K2 cell and the independent K1 blocks at `N=265,325,425` pass.

The mechanistic reading is a pair of complementary birth clocks with opposite dominant tail
curvatures and a mixing entropy that partially cancels them.  This is naturally relevant to the
canonical two-birth/topological-source decomposition, but it does not identify a continuum field.
