# P250 result: the quarter-turn root is not a scalar character

The archived N505 stream already permits the exact reconstruction `X=T+A`,
`Y=T-A`, so this tomography used no new simulation.  The scientific answer is
negative but sharp: the approximately `-pi/2` translation root cannot be
named as an independent one-dimensional spatial C4 or internal Z5 character.

The exact alphabet comparison alone is insufficient.  The second-root point
phase is `-1.5820` with jackknife SE `0.4937`; it is compatible with `-i`
(`p=0.982`) and, at this precision, also with the nearest Z5 phase
`-2 pi/5` (`p=0.510`).  Arithmetically, however, `-pi/2` is not one of the
one-step Z5 phases `2 pi j/5`.  An internal-only story would have to make one
separation step a compound/non-generator operation.

The direction data are decisive against the scalar reduction:

- the initial strict x/y character models all fail; even the best has
  `p=1.21e-63`;
- after transparently repairing the unjustified zero leading-A restriction,
  all second-mode character amplitudes look acceptable on `d1..4`
  (`p=0.094..0.517`), but every one fails held-out `d5`
  (`p=2.81e-19..1.21e-17`);
- most importantly, allowing *both* A amplitudes to be completely free at the
  frozen T roots still fails d5: `chi2=61.744/8`, `p=2.12e-10`.

Thus the failure is not “the wrong q was chosen.”  The axis-difference row is
not closed by the axis-average T-row two-state spectrum.  Picking C4, Z5, or
`q=1` by the smallest chi-square would attach a false alphabet to a projection
that is not one dimensional.

## Scientific card

- **Mechanism space changed:** scalar spatial/deck character for the second T
  root is removed on this observable and distance window.
- **Not proved:** the continuum field identity, a unique latent dimension, or
  whether the additional state is spatial, deck-internal, or mixed.
- **Observer/sector/source/geometry:** projective-leg charged pair; Z5 charges
  1/2, plus/minus hands; existing N505 Gaussian child stream; exact X/Y
  reconstruction from T/A.
- **Dependency group:** the same 80k N505 batches used by the cross-scale and
  rank-two analyses; this is not independent evidence for rank selection.
- **Next discriminator:** realize a common vector state carrying two transfer
  generators `T_x,T_y`, plus spatial rotation `R` and deck action `D`.  Mixed
  displacement rows should test `T_y=R T_x R^-1`, `T_xT_y=T_yT_x`, `R^4=I`,
  and `D^5=I`.  A bivariate block-Hankel realization is the minimal next
  object; another scalar phase vote is not.
