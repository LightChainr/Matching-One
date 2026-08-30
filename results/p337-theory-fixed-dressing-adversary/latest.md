# Theory-fixed identity-dressing adversary

No new simulation is used. The identity-dressing eigenvalues are frozen to `2^-13/8` and `2^-21/8`; only their amplitudes are fitted.

| model | q/df | GOF p | descriptive AIC | delta AIC | N1360 A_H | target SE |
|---|---:|---:|---:|---:|---:|---:|
| fixed_identity_dressing | 0.803/2 | 0.669 | 4.803 | 0.000 | -0.0005830 | 6.8e-05 |
| rank3_same_base_jordan | 0.084/1 | 0.772 | 6.084 | 1.281 | -0.0007896 | 0.00035 |
| free_lambda_recurrence | 0.077/1 | 0.781 | 6.077 | 1.274 | -0.0007956 | 0.000311 |
| free_single_lambda | 1.979/2 | 0.372 | 5.979 | 1.176 | -0.0011778 | 0.000298 |
| fixed_single_lambda0 | 15.843/3 | 0.00122 | 17.843 | 13.039 | -0.0003416 | 2.75e-05 |
| scale_neutral | 68.940/3 | 7.2e-15 | 70.940 | 66.137 | -0.0045247 | 0.000449 |

The fixed identity dressing passes GOF and is the descriptive AIC leader. Its dressing/leading magnitude halves each generation by construction and is `0.728 -> 0.364 -> 0.182 -> 0.091 -> 0.045`.

N1360 source-covariance ceiling against identity dressing:

- `rank3_same_base_jordan`: maximum `0.594` sigma; 3-sigma is impossible without reducing source uncertainty.
- `free_lambda_recurrence`: maximum `0.695` sigma; 3-sigma is impossible without reducing source uncertainty.
- `free_single_lambda`: maximum `2.308` sigma; 3-sigma is impossible without reducing source uncertainty.
- `fixed_single_lambda0`: maximum `3.878` sigma; 3-sigma requires measurement SE below `5.1e-05`.
- `scale_neutral`: maximum `9.681` sigma; 3-sigma requires measurement SE below `0.00125`.

Therefore N1360 is not yet a universal discriminator: recurrence and rank-3 Jordan are nearly forecast-identical, and free-single remains source-limited. It can efficiently reject neutral, while separating fixed-single would require unusually small measurement error. No N1360 production is started.
