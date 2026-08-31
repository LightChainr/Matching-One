# P398 eta=0 analytic response and zero-frequency susceptibility

Analytic derivative equations, evaluated in float64 on the frozen 1430-state/186-sector finite model. No new eta points or MC.

| lag | Uprime minus-plus | Uprime plus-minus | pi-only plus-minus | generator-only plus-minus |
|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 |
| 0.05 | -0.0133802074491 | -0.0581332080775 | -0.071755099596 | 0.0136218915185 |
| 0.1 | -0.0237078199818 | -0.0922982186338 | -0.116485947027 | 0.0241877283932 |
| 0.25 | -0.0412196777872 | -0.114823884919 | -0.15813577064 | 0.0433118857215 |
| 0.5 | -0.0451228476943 | -0.0680760618381 | -0.119497063464 | 0.0514210016257 |
| 1 | -0.0274468079511 | -0.00245883912064 | -0.0385847649865 | 0.0361259258659 |
| 2 | -0.00539314196998 | 0.00540524681744 | -0.00264617549784 | 0.00805142231528 |
| 4 | -0.000125271454755 | 0.000186519130542 | -9.78787482755e-06 | 0.000196307005369 |

Zero-frequency integrated Uprime (rows/columns minus, plus):

```
[[-2.65384469e-16-4.32729412e-17j -5.30034146e-02+7.09271727e-18j]
 [-5.02337767e-02-1.39204669e-16j -1.02469776e-16+3.25397182e-17j]]
```

Stationary-only contribution:
```
[[-2.80574210e-16-3.89103265e-17j -1.07689783e-01+3.10034036e-18j]
 [-1.15754092e-01-1.05495911e-16j -6.63533251e-17+1.16603927e-17j]]
```
Generator-only contribution:
```
[[ 1.51897405e-17-4.36261463e-18j  5.46863681e-02+3.99237691e-18j]
 [ 6.55203155e-02-3.37087576e-17j -3.61164514e-17+2.08793255e-17j]]
```

| model | zero-freq minus-plus | zero-freq plus-minus | max Uprime error on old lags | crossing plus-minus |
|---|---:|---:|---:|---|
| instantaneous_two_source | -0.0242880748705 | -0.117750925699 | 0.0573883 | [] |
| triplet_incidence | -0.0535346880308 | -0.0483789504958 | 0.00115135 | [{'old_lag_bracket': [1.0, 2.0], 'root': 1.0241329352575979}] |
| triplet_incidence_plus_T4 | -0.053250889263 | -0.0499567443069 | 0.000155338 | [{'old_lag_bracket': [1.0, 2.0], 'root': 1.0444784394468107}] |

| slow modes requested per ray | retained counts | zero-freq minus-plus | zero-freq plus-minus | max Uprime error |
|---:|---|---:|---:|---:|
| 1 | [1, 1] | -0.185948709864 | -0.239748051141 | 0.776858 |
| 2 | [2, 2] | 0.0537319745586 | -0.215583664121 | 0.442343 |
| 4 | [4, 4] | -0.0244964003174 | -0.0525096900714 | 0.24041 |
| 8 | [8, 8] | -0.0807720264516 | -0.0837093431267 | 0.280849 |
| 16 | [16, 16] | -0.0455778194804 | -0.0526025674085 | 0.0436921 |
| 32 | [32, 32] | -0.0530440513354 | -0.0505998021747 | 0.00123812 |
| 93 | [93, 93] | -0.0530034145838 | -0.0502337767379 | 3.13151e-15 |

Full derivative crossings: {'minus_plus': [], 'plus_minus': [{'old_lag_bracket': [1.0, 2.0], 'root': 1.047989651970103}]}

Checks: {"stationary_residual": 6.938893903907228e-16, "stationary_tangent_residual": 1.0269562977782698e-15, "stationary_derivative_sum": -1.1102230246251565e-16, "stationary_K_odd_residual": 1.8041124150158794e-16, "old_baseline_max_difference": 0.0, "same_ray_first_derivative_max": 5.054867013339658e-16, "Uprime_zero_lag_max": 0.0, "spectral_vs_Frechet_max_difference": 3.1315131483054602e-15, "spectral_vs_resolvent_integral_max_difference": 1.7677489893439373e-15}

Actual analysis wall time: 4.282644 seconds; environment: {'python': '3.9.9', 'numpy': '1.26.4', 'scipy': '1.13.1', 'machine': 'aarch64', 'hostname': '83750ac4eae34bbbbc64b894b854dd8c', 'threads': {'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'}}

This is finite-process response, not a continuum-field or square-site identification. The old geometric models consume pi and pi-prime from the full model. Modal budgets are signed and need not converge monotonically.
