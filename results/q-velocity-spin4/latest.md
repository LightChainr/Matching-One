# Q=1 spin-sector velocity oracle

Exact continuum fingerprints; no lattice target or fitted field normalization is used.

| field | family | legs | spin | x(1) | dx/dQ at 1 |
|---|---|---:|---:|---:|---:|
| loop_V_2_2 | generic_loop_primary | 4 | -4 | 17/4 | `-5*sqrt(3)/(16*pi)` |
| loop_V_2_minus2 | generic_loop_primary | 4 | 4 | 17/4 | `-5*sqrt(3)/(16*pi)` |
| loop_V_2_4_spin8_control | generic_loop_primary | 4 | -8 | 53/4 | `-41*sqrt(3)/(16*pi)` |
| loop_V_2_6_spin12_control | generic_loop_primary | 4 | -12 | 113/4 | `-101*sqrt(3)/(16*pi)` |
| thermal_Q4_epsilon | thermal_energy_descendant | 0 | 4 | 21/4 | `-9*sqrt(3)/(16*pi)` |

## Primary discriminator

- dimension gap `x_22-x_Q4 = -1`
- velocity gap `x'_22-x'_Q4 = sqrt(3)/(4*pi)`

The spin-8/spin-12 generic-loop rows are separately declared `V_(2,4)` and `V_(2,6)` controls. They are not assignments for the experiment-design H8/H12 angular aliases.

Potts multiplicities, explicit field-definition derivatives, and lattice overlaps remain unresolved inputs to any Q-score measurement.
