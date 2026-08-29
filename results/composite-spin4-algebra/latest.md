# Exact composite spin-4 algebra

## Named cases

| case | parity | q | accelerated w | harmonics | H12/H4 |
|---|---:|---:|---:|---|---:|
| `T4` | `-1` | `0` | `4` | `H4:1` | `--` |
| `T4*S0` | `-1` | `2` | `6` | `H4:1` | `--` |
| `T4*I4` | `-1` | `2` | `6` | `H0:1/2, H8:1/2` | `--` |
| `T4*I4^2` | `-1` | `4` | `8` | `H4:3/4, H12:1/4` | `1/3` |
| `T4*I4*V4` | `-1` | `14/3` | `26/3` | `H4:3/4, H12:1/4` | `1/3` |
| `T4*V4^2` | `-1` | `16/3` | `28/3` | `H4:3/4, H12:1/4` | `1/3` |

`T4*I4` is matching odd but has only H0/H8, so one even spin-4 insertion cannot
correct the H4 channel. Every listed cubic spin-4 product has exact support
`(3/4) H4 + (1/4) H12`, giving the elementary ratio `H12/H4=1/3`.

## q=3 exclusion

For generator counts `(t,i,v,s)`, the relative exponent obeys

```text
12q = 39(t-1) + 24(i+s) + 32v.
```

Matching oddness requires positive odd `t`. At `q=3`, `t>=3` already contributes
at least 78 to a target of 36. For `t=1`, division by four requires
`6(i+s)+8v=9`, an even/odd contradiction. Thus `q=3` is absent at every degree
inside the declared generator semiring, not merely in the finite enumeration.

## Additional exact warning

The ordinary thermal-tower row `q=6, w=10` is not exponent-unique after optional
analytic composites are admitted. The degree-limited artifact lists the colliding
composite monomials explicitly; harmonic sidebands and independent amplitude controls
are required to distinguish them.

## Evidence boundary

The algebra does not prove that any generator exists with nonzero lattice coupling.
`V4` parity and `S0` existence remain conditional. Continuum response tensors or
mixing may change the literal H12/H4 amplitude ratio, although the elementary harmonic
support and rational exponent arithmetic are exact.

## Reproduction

```bash
python scripts/composite_spin4_algebra.py --format json
python scripts/composite_spin4_algebra.py --format markdown
python -m unittest tests.test_composite_spin4_algebra
```
