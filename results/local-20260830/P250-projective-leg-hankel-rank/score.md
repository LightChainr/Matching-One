# P250 model-free multivariate Hankel rank score

| group | rank<=1 p | rank<=2 p | rank<=3 p | rank<=4 p | rank<=5 p | lower bound |
|---|---:|---:|---:|---:|---:|---:|
| plus_charge1 | 0 | 8.37412e-111 | 0.116602 | 0.692866 | 0.400156 | 3 |
| plus_charge2 | 0 | 3.01446e-125 | 1.22234e-07 | 0.00017114 | 0.3886 | 5 |
| minus_charge1 | 0 | 1.22642e-39 | 0.214544 | 0.114472 | 0.0528766 | 3 |
| minus_charge2 | 0 | 1.11401e-79 | 7.23594e-07 | 6.99768e-06 | 0.0125675 | 5 |
| plus_block | 2.14373e-280 | 9.27628e-191 | 1.66878e-36 | 2.8049e-11 | 0.0542623 | 5 |
| minus_block | 8.06243e-284 | 5.70084e-213 | 3.83241e-17 | 0.00235376 | 0.0655065 | 5 |
| shared_block | 1.08591e-183 | 2.32075e-151 | 4.04334e-137 | 1.79702e-13 | 2.8402e-05 | 6 |

Decision: `shared_multivariate_Hankel_rank_at_least_6`.

The primary probabilities use the finite-400-batch Hotelling correction; full residual covariances and asymptotic chi-square scores are in JSON.
