# Canonical regular-pair spatial Q activation

Decision: **contact_only_zero_rejected_at_L64_r16**.

| L | r | C | Monte Carlo SE | 99% interval |
|---:|---:|---:|---:|---:|
| 32 | 8 | 3.659179687e-05 | 1.3776e-06 | [3.300910094e-05, 4.017449281e-05] |
| 64 | 16 | 6.85546875e-06 | 6.3523e-07 | [5.203397276e-06, 8.507540224e-06] |

Primary two-sided p=1.1328120879027334e-21. Fixed ratio C64/C32: {'point_ratio': 0.18734987990392346, 'status': 'bounded', 'confidence': 0.99, 'covariance_between_sizes': 0, 'bounded_interval': [np.float64(0.14016164110332746), np.float64(0.2381648820438993)]}.

The kernel is the Q derivative of the actual connected two-insertion colour contraction. It is not a covariance of separately closed one-site marks. The full shared-component mean/covariance and input provenance are in score.json. No exponent is fitted and no continuum field is identified. This fixed-budget result receives no top-up or new completion coefficient.
