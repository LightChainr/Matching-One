# Clock-compatible integral, nontrivial thermal redistribution

The new N100 stream contains more information than its fixed-p readout. The exact integral of the odd orientation contrast is -2C; that of the even contrast is -W. Therefore using C as the shape coordinate forces the remaining odd response to have zero area, but does not force its p-profile to vanish.

The fitted clock secant is -0.277981748 +/- 0.01937. This is a post-reveal decomposition of the existing correlated block, with no new sampling.

## Thermal moments of the clock-quotient residual

Here z=N^(3/8)(p-p_ref). Integrals are evaluated algebraically from threshold-rank histograms, not by numerical quadrature. The zero odd area is imposed by an exact identity, not a failed-to-reject test.

| field | power of z | signed integral | SE |
|---|---:|---:|---:|
| A_top | 0 | 0 | 0 |
| A_top | 1 | -0.000293635291 | 4.6967e-05 |
| A_top | 2 | 0.000848557117 | 0.00010592 |
| A_top | 3 | -0.00102589827 | 0.00010415 |
| A_top | 4 | 0.00296054821 | 0.00029904 |
| E_top | 0 | -0.000157451869 | 6.5885e-05 |
| E_top | 1 | 0.000420769264 | 7.1136e-05 |
| E_top | 2 | -0.000534819096 | 6.324e-05 |
| E_top | 3 | 0.00158761919 | 0.00017692 |
| E_top | 4 | -0.00203922427 | 0.0001916 |

## Empirical profile structure

- A_top: core zeroes [0.42427987773344655, 0.7719349733430698]; largest mean lobe 0.00149561 near p=0.32125; positive/negative areas 0.000347135/-0.000347135.
- E_top: core zeroes [0.42493013060374935, 0.5971521948240728, 0.7714296980272011]; largest mean lobe -0.00149502 near p=0.32125; positive/negative areas 0.000234364/-0.000391816.

Roots, quadrature lobe areas and peak locations are descriptive. All pointwise readouts share one dependency group; the JSON retains their full compact-grid covariance. A separately propagated fixed-E4 profile is included for comparison.

Clock calibration is data-adaptive; uncertainties use the correlated ratio influence. No new field count, exact E4 identity, independent p-grid evidence, or continuum exponent is inferred.
