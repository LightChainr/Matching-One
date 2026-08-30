# Four-generation covariance-aware H4 recurrence

No new simulation is used. Exact H4 geometry signs are removed before fitting the one Gaussian lineage.

Fixed `lambda0=2^-13/8`; fit `lambda1=0.270681 +/- 0.128009` with 95% profile interval `[2.78862e-05,0.641662]`. The lower endpoint is boundary-close and only marginal evidence for positive lambda1.

Recurrence GOF is `0.077/1` (`p=0.781`). Correction/leading magnitude falls `0.870 -> 0.727 -> 0.607 -> 0.506 -> 0.423` through predicted N1360.

| model | q/df | GOF p | descriptive AIC | delta AIC | N1360 A_H | SE |
|---|---:|---:|---:|---:|---:|---:|
| fixed_lambda0_single | 15.843/3 | 0.00122 | 17.843 | 11.863 | -0.0003416 | 2.75e-05 |
| free_single_lambda | 1.979/2 | 0.372 | 5.979 | 0.000 | -0.0011778 | 0.000334 |
| scale_neutral | 68.940/3 | 7.2e-15 | 70.940 | 64.961 | -0.0045247 | 0.000449 |
| fixed_lambda0_plus_correction | 0.077/1 | 0.781 | 6.077 | 0.098 | -0.0007956 | 0.000299 |

The recurrence passes the genuinely heldout fourth-generation shape. Free-single also has acceptable GOF and is essentially AIC-tied; fixed single H4 and scale-neutral fail. AIC is descriptive, not the core claim.

N1360 forecasts and exact `(4,340)` child geometry are frozen in the machine-readable result. No N1360 production is authorized or running.
