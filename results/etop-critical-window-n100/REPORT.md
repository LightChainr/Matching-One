# N100: locate the thermal redistribution before assigning its mechanism

This consumes PR484's archived 200 aligned batches. No new sampling, curve fitting, model vote or transfer-engine run is added.

## Primary descriptive window and odd dipole

The center window is p in [0.414918110, 0.770573992], or |z|<=1 with z=N^(3/8)(p-p_ref). The exponent is a coordinate convention, not a new estimate.

| region | integral z R_A dp | delete-one SE |
|---|---:|---:|
| full | -0.000293635291 | 4.69623e-05 |
| w1_lower | -0.0004198789744 | 4.95388e-05 |
| w1_core | -1.011331359e-05 | 2.24112e-05 |
| w1_upper | 0.000136356997 | 2.0759e-05 |

These are additive **signed** contributions, not shares of positive signal mass. They retain all cross-region covariance and the uncertainty in the same-stream clock ratio.

## First/second activation in the center window

| readout | integral R dp | delete-one SE |
|---|---:|---:|
| A | -0.0003460293249 | 4.62677e-05 |
| E | 1.619451641e-06 | 5.54905e-05 |
| F1 | -0.0001738243883 | 3.50209e-05 |
| F2 | -0.0001722049366 | 3.71953e-05 |

F1=(A-E)/2 and F2=(A+E)/2; the birth rows are a change of coordinates, not two independent replications.

## Width sensitivity of the odd dipole

| z half-width | lower | core | upper |
|---|---:|---:|---:|
| 0.5 | -0.00037742089 ± 5.35e-05 | 1.1898858e-06 ± 5.85e-06 | 8.2595712e-05 ± 2.5e-05 |
| 1 | -0.00041987897 ± 4.95e-05 | -1.0113314e-05 ± 2.24e-05 | 0.000136357 ± 2.08e-05 |
| 1.5 | -0.00029931901 ± 2.81e-05 | -5.8312511e-05 ± 3.41e-05 | 6.3996232e-05 ± 6.23e-06 |

## Interpretation and next discriminant

This is a localization of the existing finite-N response. A nonzero center-window response does not identify a critical field; a larger outer contribution does not establish a thin-geometry mechanism. The next scale comparison should preserve the three homothetic shapes and these named regions, then ask whether the distribution contracts in p or remains at an off-critical location. PR484's N400 design is an existing acquisition option, not a new approval gate.

No transport test is repeated: PR484 already contains the failed joint finite-Jacobian transport comparison. A-only quantile transport is not a discriminating model.

## Reproduction and dependence

Source: open PR484, `894b3d800c5aeaad3dd8b0f893b6f17d85d234c6`; all three shape pairs share the same 2,000,000 permutation counters, seed 20260831125401, offset 267100000000. Histograms and their SHA256 values are listed in latest.json.

Clock ratio -0.2779817479 ± 0.0193769. Each delete-one removes one aligned batch from all three shape pairs and refits this ratio. JSON stores all 200 vectors, the complete covariance, and window bounds; singular directions and overlapping windows are preserved.

The Bernstein integrals use incomplete-beta identities evaluated in float64, not quadrature and not a rigorous numeric certificate. This follows disclosure of the full-p curve; all window claims are retrospective.

```bash
git fetch origin analysis/etop-modulus-survivors-20260831
python3 scripts/analyze_etop_critical_window.py --output-dir /path/to/fresh-output
```
