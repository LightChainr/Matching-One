# Four-channel operator-mixing parity audit

The matching involution fixes an exact zero pattern, but it does not reduce the four
allowed Taylor coefficients to two amplitudes.

| channel | nonzero coefficient | L power | N power |
|---|---|---:|---:|
| P4_S | `I:f0` | `-2` | `-1` |
| P4_D | `T:f0` | `-13/4` | `-13/8` |
| P4_S_prime | `T:f1` | `-5/2` | `-5/4` |
| P4_D_prime | `I:f1` | `-5/4` | `-5/8` |

## Exact identifiability conclusion

With columns `[I:f0, I:f1, T:f0, T:f1]`, the structural matrix has rank `4`.
Therefore parity supplies zeros and channel assignments, not the ratios `f_I1/f_I0` or
`f_T1/f_T0`. A genuine two-amplitude joint prediction needs those two relations from an
independent dynamical calculation or frozen training data.

## Boundary

This exact audit fixes selection zeros and powers only. It does not fit amplitudes, identify fields, or supply the missing Taylor-ratio dynamics.
