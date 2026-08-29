# P275 Cartan contact background

This is a post-reveal mechanism analysis. It does not modify or replace the frozen P275 score.

## Exact regression constraint

Commit `8bef10b` establishes the connected rank-gate identity

`Cov_p(q,J_D) = E_p[J_S]/2 + (p-1/2-E_p[q])E_p[J_D]`.

Here it is used as a numerical regression constraint, not reproved as a new result. The current Gamma is a contact term, not a contact term plus an identifiable small remainder.

## Nine-geometry reproduction

| geometry | Re Gamma | Re thermal | Re relative | max component residual |
|---|---:|---:|---:|---:|
| N50/i | 0.3426131802 | 0.3420694833 | 0.000543696845 | 2.34e-17 |
| N50/2i | 0.4341650783 | 0.4347767266 | -0.0006116483383 | 5.55e-17 |
| N50/5i_over_2 | 0.4651537308 | 0.4654385205 | -0.0002847896773 | 5.55e-17 |
| N130/i | 0.3357896157 | 0.3356953227 | 9.429307881e-05 | 0 |
| N130/2i | 0.4313774845 | 0.4311814424 | 0.0001960420846 | 0 |
| N130/5i_over_2 | 0.4621522334 | 0.4623691352 | -0.0002169017845 | 8.33e-17 |
| N170/i | 0.3349478329 | 0.3348235845 | 0.0001242483416 | 1.11e-16 |
| N170/2i | 0.4308121757 | 0.4306990862 | 0.0001130894903 | 5.55e-17 |
| N170/5i_over_2 | 0.4622541533 | 0.4621266978 | 0.0001274555338 | 0 |

The relative-source term is statistically resolved but small: its largest complex-magnitude fraction is `0.1587%`. The order-one background is therefore specifically thermal-contact dominated.

## Post-reveal discovery GLS

| model | chi2 | dof | survival p |
|---|---:|---:|---:|
| constant_by_modulus | 9275.92 | 12 | 0 |
| constant_by_modulus_plus_Q4_shape_tail | 4184.73 | 10 | 0 |
| constant_by_modulus_plus_free_tail | 74.745 | 6 | 4.33141e-14 |

These fits describe the scale-zero contact shape and finite-size drift only. They are not a re-score of P275.
Even the free per-modulus N^-13/8 tail fails, so the same statistic does not support a constant-plus-small-field decomposition.

## Scientific layers

- **exact:** The observed Gamma is exactly the Cartan contact term at every finite N; contact subtraction leaves algebraic zero.
- **continuum_shape:** A nonzero scale-zero limit by modulus is compatible with a conditional/projective-line polarization. Its limiting modulus function is inferred, not fixed by the Ward identity.
- **N_minus_13_over_8:** No independent remainder is identifiable in this same-site statistic: the exact identity annihilates it before asymptotics.
- **exploratory_conjecture:** Any genuine H4 field propagation must be measured with a separated or typed source/readout rather than a global same-site q-J_D covariance.

## Scientific card

- Question: why does revealed Gamma approach a nonzero modulus-dependent constant?
- Exact: three-state Cartan increment algebra plus Bernoulli edge balance makes Gamma a finite-N contact Ward term.
- Discovery: nine-geometry GLS separates modulus constants from optional N^-13/8 tails using the full covariance.
- Boundary: the same Gamma contains no separately identifiable H4 remainder after exact contact subtraction.
- Next test: freeze a background-annihilated held-out size and a separated typed-source observable; never reinterpret the frozen selector.
