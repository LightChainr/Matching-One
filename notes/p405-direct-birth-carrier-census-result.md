# P405 exact priority-weighted direct-birth carrier census

## Decision

The theta and figure-eight direct-birth channels are now separated exactly by
predecessor size and converted from directed insertion edges to uniform random
permutation probabilities.

For carrier type `c`, the checked identity is

```text
D_c(N) = sum_k E_c(k) / (N * binom(N-1,k))
       = sum_v integral_0^1 P_p(E_c(v,A)) dp.
```

The resulting finite-volume masses are:

| quotient | N | theta probability | figure-eight probability | figure-eight share of D |
|---|---:|---:|---:|---:|
| `2+i` | 5 | 0 | 0 | -- |
| `2+2i` | 8 | `17/105` | 0 | 0 |
| `3` | 9 | `1/14` | `1/70` | `1/6` |
| `3+i` | 10 | `5/63` | 0 | 0 |
| `3+2i` | 13 | `304/3465` | 0 | 0 |
| `4` | 16 | `2644/45045` | `1/273` | `165/2809` |
| `4+i` | 17 | `10411/180180` | `61/90090` | `122/10533` |

For every geometry the two typed masses sum exactly to the previously
certified total direct-birth probability.

## Edge shares are not path probabilities

The raw figure-eight edge shares and the priority-weighted probability shares
are different:

| N | raw edge share | permutation-probability share |
|---:|---:|---:|
| 9 | `9/45 = 0.2` | `1/6 = 0.1667` |
| 16 | `336/4624 = 0.07266` | `165/2809 = 0.05874` |
| 17 | `119/8823 = 0.01349` | `122/10533 = 0.01158` |

The difference comes entirely from the exact predecessor-size Beta weights.
Counting directed subset edges without them does not estimate the probability
seen by a uniform site-order filtration.

## Geometry dominates this tiny-volume table

The figure-eight channel is present on the square quotients N9 and N16,
absent on the nearby oblique quotients N10 and N13, and much smaller on N17
than N16.  This is an exact finite-geometry selection effect, not evidence for
or against the candidate `N^(-2)` figure-eight asymptotic.

Consequently, the proposed relative law

```text
D_figure8 / D_theta ~ N^(-7/6)
```

cannot be tested by pooling these quotient shapes.  It needs a matched-shape
lineage with the carrier tag frozen before reveal.  The present census is the
exact tiny-control and weighting oracle for that acquisition.

## Scientific card

- Mechanism space changed: the nonnegative figure-eight correction amplitude
  is an exactly measurable path channel, not a free scalar nuisance term.
- Not proved: either arm exponent, a stable typed/untyped landing frequency,
  or a regular-variation amplitude.
- Observer/sector/source/geometry: direct `0->2` ambient-rank birth under a
  uniform random site order; theta and figure-eight carrier topology; seven
  Gaussian square-site quotients through N17.
- Dependency group: exhaustive subset enumeration and universal-cover carrier
  classification; no Monte Carlo.
- Next lift: record theta/figure-eight tags jointly on one matched-shape
  production lineage, preserving within-path covariance and per-k counts.
