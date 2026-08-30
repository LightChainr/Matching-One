# Scientific card: the second child returns toward nominal H4 decay

- Mechanism changed: N340 `(18+4i,14+12i)` restores the exact positive H4
  pair direction, `Delta K_A=+0.00774472 +/- 0.00199130` (`3.889 sigma`).  The
  sign is right, although the preregistered 5-sigma resolution gate was not met
  because the observed N340 variance exceeded the N170 projection.
- Fixed scale discriminator: `A_H=-0.00485726 +/- 0.00124889`.  It is `-1.005`
  measurement SE from nominal `2^-13/8` decay, `+2.270` SE from the secondary
  N85-to-N170 effective continuation, and `+5.008` SE from scale-neutral.
- Source uncertainty: the effective-continuation target inherits the noisy N85
  source, reducing its predictive residual to `0.942` SE.  Therefore this block
  sharply rejects no decay and points back toward nominal H4, but does not by
  itself separate nominal H4 from a curved two-state continuation.
- Orthogonal control: the charged/projective scalar is
  `+0.00106775 +/- 0.00211093` (`0.506 sigma`), again consistent with zero.
- Mechanism-space update: the N170 excess amplitude does not persist as a
  scale-neutral state.  The same-lineage sequence bends toward nominal H4
  decay while retaining the exact alternating geometry sign.
- Does not prove: a continuum exponent, a unique correction law, or an H4/H8
  distinction.  No exponent was fit and no harmonic vote was reopened.
- Observer/sector/source/geometry: `K_A=d_eta log W_A`, second `1+i` Gaussian
  child, Smith `(2,170)`, 12M fresh samples/shape, 80 aligned batches and full
  covariance.
- Dependency: all targets were frozen in `5369c21`; the N340 seed/counter block
  is independent of N85/N170.  The effective target is explicitly secondary.
- Next discriminator: use the exact same lineage to fit a correction coordinate
  only after treating N85 source noise and N-dependent variance; do not add a
  generic fourth geometry merely to refit an exponent.
