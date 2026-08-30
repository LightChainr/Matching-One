# P250 augmented old-plus-degree-five joint operator

This is a zero-new-sample combination of the pinned N505 radius-four 80k/400
batch stream and the independent radius-five 1.2M/400 batch stream.  It uses
one `40 x 6` candidate matrix and one 70-real-coordinate rank-at-most-five
Wald gate rather than combining the two published p-values.  The fixed-map
family predates the fresh reveal, but this augmented gate was designed after
the marginal R2/R3 conflict was observed.  It is a primary **retrospective**
zero-sample synthesis, not a prospective confirmation.

All five maps frozen before the radius-five reveal reject.  Their finite-batch
Hotelling p-values are:

| map | p |
|---|---:|
| identity + conjugation | `3.670e-14` |
| Alexander R0 + conjugation | `2.463e-12` |
| Alexander R1 + conjugation | `1.031e-6` |
| Alexander R2 + conjugation | `3.550e-8` |
| Alexander R3 + conjugation | `1.240e-8` |

All eleven secondary maps already frozen in the joint-annihilation manifest
also reject.  This is not merely a cross-hand bridge failure: the plus-hand
`20 x 6` full-Schur rank-at-most-five gate rejects (`p=6.007e-7`), and no
primary coordinate has both hand-specific gates survive.  The correct
result is therefore `declared_fixed_bridges_rejected_and_full_hand_rank5_support_incomplete`,
not a selected R2/R3 bridge or a supported general `5+5` fallback.

Every candidate exactly replays its old 38-real-coordinate point and
covariance block before the extension rows are appended.  The largest point
replay error is `6.51e-19`; the largest covariance replay error is `7.94e-22`.
Old and fresh delete-one covariance matrices are stored separately and added
once.  Equal batch numbers across the two independently seeded streams are
never paired.

`score.json` is the compact decision summary.  `covariance.npz` contains all
residual points and complete old/fresh covariance components; every total is
derived as `old + fresh` and is not duplicated.  Shared hand blocks are stored
once: one canonical plus block and eight minus-geometry blocks are referenced
by the sixteen candidate summaries.  The NPZ SHA256 is pinned in the summary.

The input audit pins both response JSONs, both exact gates and the runner
sources.  It verifies common N505 geometry, `p=0.59274605079`, hands
`plus/minus`, charges `1/2`, affine-fiber C4 gauge, translated-origin
convention, intentionally distinct translation salts, and the five explicit
radius-five aliases.

The earlier radius-five support result asked whether old SVD annihilator-line
directions were compatible with fresh rows.  This analysis asks the stronger
locally complete question by retaining the old residual magnitude and its
cross-covariance with the extension.  It does not support the current general
`5+5` fallback because the plus-hand full-Schur gate fails, but it also does
not identify rank six, a closed transfer algebra, a microscopic quotient, or
a continuum field.  The saved
finite-batch Hotelling calibration follows the existing project approximation
for summed independent jackknife covariance estimates; it is not an exact
two-sample pivot law.
