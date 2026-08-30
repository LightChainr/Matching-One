# N680 same-lineage heldout reveal

The production was generated only after preregistration commit `ba4ca6f`.

Observed H4 amplitude: `-0.002167556 +/- 0.000557`. Exact pair: `-0.003456090 +/- 0.000888` (`z=-3.892` versus zero).

| frozen forecast | target | residual / measurement SE | residual / predictive SE | full vector q/2 |
|---|---:|---:|---:|---:|
| two_mode_recurrence | -0.001841297 | -0.586 | -0.295 | 0.280 |
| single_frozen_lambda0 | -0.001027154 | -2.048 | -2.024 | 4.987 |
| single_free_lambda | -0.003037819 | +1.563 | +0.747 | 0.652 |
| scale_neutral | -0.008911964 | +12.110 | +7.159 | 51.748 |

Projective scalar: `-0.000311650 +/- 0.000785` (`z=-0.397`).

Closest frozen point forecast: `two_mode_recurrence`. Source-fit uncertainty is retained separately in the predictive column; it was not used to alter the frozen model order.

Continuity closes to `1.05e-13`. No model, exponent, harmonic, or basis was refit.
