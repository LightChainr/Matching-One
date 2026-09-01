# P275 rho-child observer-model elimination

Decision: `ALL_THREE_DECLARED_PRODUCTION_PARAMETERIZATIONS_EXCLUDED`.

This is a post-reveal reanalysis of one frozen N112 square-bond dependency block. It generated no new Monte Carlo samples.

## Production reconstruction

- Commit: `2402a3330b421595d3573337a5723ff3dbdcb7e9`
- Dependency group: `p267-rho-C3-Etop-N112-fresh2m`
- Batches / samples: `100` / `2,000,000`
- Covariance condition number: `39.483433`
- The nine-dimensional mean and full covariance of the mean were reconstructed directly from the pinned `batches.csv`; the archived Etop-r1 and determinant points were replayed against the pinned `score.json`.

## Declared model scores

| parameterization | T2 / constraints | Hotelling F reference p | decision at 0.01 |
|---|---:|---:|---|
| `normalizer_only` | 215.105242 / 6 | 1.52127e-21 | excluded |
| `rank1_mass_only` | 334.368604 / 6 | 5.16481e-28 | excluded |
| `independent_real_rescalings` | 19.087779 / 3 | 0.000646405 | excluded |

The fixed probability identities rejected here are: common-denominator-only scaling, rank-1-total-mass-only reweighting with fixed internal winding composition, and the broader class in which E is free but primitive H4 may only undergo a real rescaling on each child.

## Interpretation boundary

All three scores reuse the same 100 aligned batches and form one evidence unit. The Hotelling values are finite-batch references under iid Gaussian batch means; the hypotheses were defined after reveal and their p-values are not prospective or additive.

These exclusions are parameterization-specific. `E_top` is an Alexander-even topology coordinate, not an identified energy operator; primitive `H4` is a direction-weighted observer, not a local spin-4 identification. The result does not decide square-site original U, Q4/Jordan identity, H4 versus H8, or an asymptotic exponent.
