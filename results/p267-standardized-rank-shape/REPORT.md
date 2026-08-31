# Standardized rank shape beyond center and width

| Frozen readout | N100 | N400 | N400-N100 | SE |
|---|---:|---:|---:|---:|
| rank_step.standardized_mu3 | 0.40772305 | 0.20547979 | -0.20224326 | 0.068732 |
| rank_step.standardized_mu4 | 1.6271517 | 1.5751699 | -0.051981739 | 0.061799 |
| rank_step.standardized_mu5 | 1.1289886 | 0.48496291 | -0.64402572 | 0.21635 |
| rank_step.standardized_mu6 | 3.2865657 | 3.2370862 | -0.049479532 | 0.25902 |
| canonical.standardized_mu3 | 0.3655846 | 0.1908877 | -0.1746969 | 0.065907 |
| canonical.standardized_mu4 | 1.702392 | 1.6293852 | -0.073006848 | 0.05933 |
| canonical.standardized_mu5 | 1.0797351 | 0.46398504 | -0.61575008 | 0.21681 |
| canonical.standardized_mu6 | 3.6382149 | 3.4840885 | -0.15412648 | 0.26165 |
| canonical.standardized_p.first_peak | -0.6488184 | -0.69334776 | -0.044529367 | 0.04023 |
| canonical.standardized_p.valley | 0.40228288 | 0.30509793 | -0.097184946 | 0.069876 |
| canonical.standardized_p.second_peak | 1.2602201 | 1.0644501 | -0.19576997 | 0.047476 |
| canonical.unit_area_scaled_height.first_peak | 0.59656234 | 0.57498562 | -0.021576723 | 0.033692 |
| canonical.unit_area_scaled_height.valley | 0.030522946 | 0.014480021 | -0.016042925 | 0.040506 |
| canonical.unit_area_scaled_height.second_peak | 0.43231157 | 0.49582094 | 0.063509365 | 0.033888 |
| canonical.valley_over_first_peak | 0.051164721 | 0.025183275 | -0.025981446 | 0.07158 |
| canonical.second_over_first_peak | 0.72467124 | 0.86231885 | 0.13764761 | 0.059885 |

## Fixed blocks

- rank_step_mu3_to_mu6_primary: chi2=73.739323/4, nominal p=3.6813e-15.
- canonical_mu3_to_mu6_smoothing_control: chi2=111.37854/4, nominal p=3.69792e-23.
- canonical_ordered_peak_coordinates_auxiliary: chi2=53.034018/6, nominal p=1.15493e-09.

Fixed orders and ordinal landmarks, not selected windows. N100/N400 inputs had previously informed the width conjecture; this is a mechanism decomposition, not independent evidence from new samples. Signed-profile moments, not probability-law cumulants. The six peak coordinates exclude redundant height ratios from their omnibus score. Canonical peaks retain finite-N smoothing and are not literal extrema of the noisy rank-step function.
