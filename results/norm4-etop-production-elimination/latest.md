# Production E_top model elimination audit

The existing 6-size threshold-rank production block is scored
against five declared model images below. The table records which images
are excluded on this block; non-exclusion is not model confirmation.

No Monte Carlo samples are generated here. The exact state transform is

```text
A_top = delta_F1 + delta_F2
E_top = delta_F2 - delta_F1
```

and every score uses the pinned full cross-size aligned-delete-one covariance.

| model | fitted parameter | chi-square / df | survival p | alpha=.05 |
|:--|:--|--:|--:|:--|
| `M0_PURE_ALEXANDER_ODD` | none | 5324.014637 / 6 | 2.84856153229e-1150 | excluded |
| `M1_SECOND_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO` | none | 3761.100630 / 6 | 3.43065482651e-811 | excluded |
| `M2_FIRST_ACTIVATION_DIRECTIONAL_RESPONSE_ZERO` | none | 21141.018681 / 6 | 1.07985761793e-4583 | excluded |
| `M3_COMMON_PROJECTIVE_RANK_PLANE_LINE` | lambda=-0.378145886 | 106.664717 / 5 | 2.07518807199e-21 | excluded |
| `M4_SINGLE_FIXED_H4_POWER` | c=-0.280066348 | 177.527945 / 5 | 1.80440349908e-36 | excluded |

## Direction-normalized E_top production coordinates

| N | E_top | SE | dependency group |
|--:|--:|--:|:--|
| 65 | -2.890888046e-04 | 5.378e-06 | `crn-2026104501-5100000000-7000000000` |
| 85 | -2.190824903e-04 | 5.128e-06 | `crn-2026104501-5100000000-7000000000` |
| 130 | -1.360201721e-04 | 6.646e-06 | `crn-2026104501-5100000000-7000000000` |
| 170 | -1.092840259e-04 | 5.357e-06 | `crn-2026104501-5100000000-7000000000` |
| 260 | -7.513283594e-05 | 1.017e-05 | `crn-2026105401-8200000000-9200000000` |
| 340 | -5.913167138e-05 | 8.349e-06 | `crn-2026105402-8200000000-9200000000` |

## Certificate form and boundary

The pure-odd distance is `5324.014637`. At the
stored strong reference alpha `1.0e-12`, the
chi-square critical value is `68.104748`
and the separation margin is `5255.909888`.

The redundant matching-coordinate audit has maximum mean residual `5.421e-20` and covariance-row residual `9.500e-25`.

This artifact is a hash-bound exact linear transform followed by a
floating jackknife Mahalanobis confidence-set separation. It is not an
exact probability bound, an interval-reconstructed LDL certificate, or
an SOS certificate. M3 excludes only one common projective line across
the declared sizes/geometries; M4 excludes only one uncorrected fixed-power
amplitude. Neither result excludes every H4 or multi-field mechanism.
The F1/F2 rows are activation-resolved directional responses; their
nonzero values do not assert the mere existence of K1 or K2, which was
already part of the input construction.
This completes only the canonical Phase-D E_top production scoring. It
does not implement the Phase-E `J_top` versus `J_bulk` test and is not
the proof-carrying outward-interval/SOS certificate proposed in #370.
