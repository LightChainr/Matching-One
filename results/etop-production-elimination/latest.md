# Production E_top model elimination

The existing ten-size threshold-rank production block rejects a pure
Alexander-odd state response. Both activation-resolved directional H4
response components are nonzero on the declared production block,
and neither one common `E=lambda A` line nor one uncorrected
`E=c N^(-13/8)` amplitude describes the complete declared archive set.

No Monte Carlo samples are generated here. The exact state transform is

```text
A_top = delta_F1 + delta_F2
E_top = delta_F2 - delta_F1
```

and every score uses the pinned full cross-size aligned-delete-one covariance.

| model | fitted parameter | chi-square / df | survival p | alpha=.05 |
|:--|:--|--:|--:|:--|
| `M0_PURE_ALEXANDER_ODD` | none | 445.618411 / 10 | 1.79698e-89 | excluded |
| `M1_SECOND_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO` | none | 182.904518 / 10 | 5.84193e-34 | excluded |
| `M2_FIRST_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO` | none | 1041.048616 / 10 | 2.67955e-217 | excluded |
| `M3_COMMON_PROJECTIVE_RANK_PLANE_LINE` | lambda=-0.432710041 | 28.593006 / 9 | 0.000758788 | excluded |
| `M4_SINGLE_FIXED_H4_POWER` | c=-0.311937114 | 37.482032 / 9 | 2.15838e-05 | excluded |

## Direction-normalized E_top production coordinates

| N | E_top | SE | dependency group |
|--:|--:|--:|:--|
| 65 | -2.822838027e-04 | 2.927e-05 | `crn-2026104501-5000000000-5100000000` |
| 85 | -2.173859034e-04 | 2.301e-05 | `crn-2026104501-5000000000-5100000000` |
| 130 | -1.336148057e-04 | 2.747e-05 | `crn-2026104501-5000000000-5100000000` |
| 145 | -1.169780619e-04 | 2.103e-05 | `crn-2026105003-7000000000-7100000000` |
| 170 | -1.065363607e-04 | 2.605e-05 | `crn-2026104501-5000000000-5100000000` |
| 185 | -1.032326313e-04 | 1.373e-05 | `crn-2026104301-7000000000-7500000000` |
| 265 | -8.156724854e-05 | 1.090e-05 | `crn-2026104301-7500000000-8000000000` |
| 290 | -4.192429588e-05 | 2.217e-05 | `crn-2026105004-7000000000-7100000000` |
| 325 | -4.542720661e-05 | 2.387e-05 | `crn-2026105701-10000000000-10500000000` |
| 425 | -4.663976706e-05 | 2.061e-05 | `crn-2026105701-10500000000-11000000000` |

## Certificate form and boundary

The pure-odd distance is `445.618411`. At the
stored strong reference alpha `1.0e-12`, the
chi-square critical value is `78.471647`
and the separation margin is `367.146765`.

The redundant matching-coordinate audit has maximum mean residual `5.421e-20` and covariance-row residual `2.895e-24`.

This artifact is a hash-bound exact linear transform followed by a
floating jackknife Mahalanobis confidence-set separation. It is not an
exact probability bound, an interval-reconstructed LDL certificate, or
an SOS certificate. M3 excludes only one common projective line across
the declared sizes/geometries; M4 excludes only one uncorrected fixed-power
amplitude. Neither result excludes every H4 or multi-field mechanism.
The F1/F2 rows are activation-resolved directional responses; their
nonzero values do not assert the mere existence of K1 or K2, which was
already part of the input construction.
