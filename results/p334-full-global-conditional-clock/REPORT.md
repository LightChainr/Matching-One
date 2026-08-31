# Complete A_top under the paired prefix-safe conditional policy

## N325

{'eligible_with_R1_replacement': 9055, 'both_R2_no_change': 1702, 'blocked_any_R0': 9207, 'blocked_solver_failure': 36}

- baseline p_ref A_H4 = 0.0001691064399 +/- 0.00564169
- baseline p_integral A_H4 = 0.000110893907 +/- 0.000775132
- safe p_ref A_H4 = -0.000119738456 +/- 0.00534348
- safe p_integral A_H4 = 0.0001086716523 +/- 0.000751799

| Observer | canonical between-R1-state fraction | integral between-R1-state fraction |
|---|---:|---:|
| old_gated_R1_F2 | 84.23901% +/- 0.2844pp | 99.84895% +/- 0.0037498pp |
| baseline_full_A | 0.02036156% +/- 0.023125pp | 0.01849675% +/- 0.020332pp |
| safe_full_A | 0.02038979% +/- 0.02376pp | 0.02229144% +/- 0.022583pp |

Safe/baseline variance ratios (canonical, integral): [0.9765503805049196, 0.9540506150998739]
Paired residual-noise fractions: [0.023171717210361577, 0.046318721586539116]

## N425

{'eligible_with_R1_replacement': 8903, 'both_R2_no_change': 1572, 'blocked_any_R0': 9413, 'blocked_solver_failure': 112}

- baseline p_ref A_H4 = 0.009260019261 +/- 0.00746148
- baseline p_integral A_H4 = 0.0007625364637 +/- 0.000879756
- safe p_ref A_H4 = 0.01010731787 +/- 0.00694886
- safe p_integral A_H4 = 0.0009430344634 +/- 0.000819725

| Observer | canonical between-R1-state fraction | integral between-R1-state fraction |
|---|---:|---:|
| old_gated_R1_F2 | 83.07235% +/- 0.30556pp | 99.8638% +/- 0.0037413pp |
| baseline_full_A | 0.007595141% +/- 0.012411pp | 0.005108005% +/- 0.010798pp |
| safe_full_A | 0.007138115% +/- 0.013274pp | 0.005454106% +/- 0.011425pp |

Safe/baseline variance ratios (canonical, integral): [0.977647145592864, 0.9553890395929978]
Paired residual-noise fractions: [0.023229174641186415, 0.04611210774425943]

Full A_top on the original paired counter blocks. Baseline, safe hybrid, gated strata and centering are correlated readouts of the same archive, never independent evidence or separately summed errors. Four-state fractions condition on the old R1 flags only. No high-dimensional covariance inverse or new omnibus test.
