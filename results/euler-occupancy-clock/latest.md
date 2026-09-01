# Euler occupancy-clock and source residual

This is a zero-new-sample reanalysis of the committed N325/N425 path
aggregates. `n_occ` is the pre-insertion occupation count, not K1/K2.

## Decision

The residual after the exact occupancy clock and the fixed-`n_occ` JS source projection has joint `chi2=390067/4`, `p=1.97418090738e-84697` (`log10 p=-84696.705`).

| N | residual complex | zero p | clock / raw | residual / raw | residual / fixed-`n_occ` |
|---:|---:|---:|---:|---:|---:|
| 325 | `-10.192415+10.25924i` | `5.19921607448e-42107` | `0.531395` | `0.466395` | `0.995272` |
| 425 | `-11.907344-7.8519268i` | `1.94687257449e-42596` | `0.561588` | `0.436546` | `0.996019` |

## Exact population decomposition

For each orientation and delete-one replicate:

```text
O_ext = mu_ext(n_occ) + fixed-n_occ residual,
mu_ext(k) = k - 2N (k)_2/(N)_2 + N (k)_4/(N)_4.
```

The frozen radius-2 local nuisance also has an exact root-absent
conditional mean on these locally injective period quotients. The report
then centers JD and JS within each `n_occ`, recomputes their same-next-site
Gram coefficient, and scores the far-Euler coupling to the remaining source.

The finite archive records the small sample closure residual between the
analytic clock and the empirical fixed-occupation covariance in the JSON.
The clock, within-occupation and source-projected rows are coordinates of
the same production block, not independent votes. A surviving residual
escapes only the declared `sigma(n_occ)+span(JS)` nuisance. It is not a
Q4 field identification, exponent, or proof of a second microscopic source.
