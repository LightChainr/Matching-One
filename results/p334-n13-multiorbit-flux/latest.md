# Exact N=13 multi-orbit projective birth/exit flux

Status: `exact_multiorbit_birth_exit_flux`.

Gaussian quotient 3+2i (N=13) has two inequivalent primitive-line orbits with exactly opposite chi4. Both have nonzero birth and exit flux, and their orbit-resolved birth-minus-exit imbalances reinforce the same total dA4/dp direction at p_ref.

## Exact frontier

- subset states: 8,192
- directed addition edges: 53,248
- orbit count: 2
- coefficientwise source/sink identity: `True`

| orbit | primitive lines | chi4 | birth edges | exit edges |
|---|---|---|---:|---:|
| axis_orbit | [[0, 1], [1, 0]] | (-119/169, 120/169) | 8554 | 4602 |
| diagonal_orbit | [[1, -1], [1, 1]] | (119/169, -120/169) | 858 | 650 |

## At p_ref = 0.59274605079

Total dA4/dp = `{'real': '-0.0885811311990293042', 'imag': '0.0893255104528026597'}`.

- axis_orbit: birth-exit incidence 0.0950721525657492276; signed reinforcing share 0.755739917417081006.
- diagonal_orbit: birth-exit incidence -0.0307279413219478515; signed reinforcing share 0.244260082582918994.

The diagonal orbit has the opposite chi4 and the opposite incidence imbalance, so its contribution reinforces rather than cancels the axis-orbit H4 derivative.

## Boundary

- This localizes the exact finite-volume H4 derivative into line orbits; it does not identify a continuum field.
- The p_ref shares are an exact evaluation coordinate, not a fitted asymptotic amplitude.
- No 13! ordering enumeration, Monte Carlo sample, or new Huawei production is used.
