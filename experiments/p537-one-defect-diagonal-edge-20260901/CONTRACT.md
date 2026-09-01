# Frozen full-root one-defect diagonal-edge contract

Status: `FROZEN_BEFORE_OUTPUT`

This consumes the existing exact N25 radius-one collar fibres.  It does not
generate a new population and does not re-estimate a root or Schur coefficient
inside a cell.

A physical site-flip row is a **diagonal one-defect edge** exactly when:

1. the thermal site has the alternating collar landing;
2. the digital rank changes (`rank0 != rank1`);
3. the source is present; and
4. the canonical source Bell state changes (`bell0 != bell1`).

No corner word, source-contact word, source orbit, or geometry is excluded.
The complete axis/tilted pooled-root `R`, componentwise `beta`, source means,
and Schur allocation are imported unchanged from the complete population.

The first lexicographic coarse class with a signed interval excluding zero is
then refined by exact `(bell0,bell1)`.  The first lexicographic Bell transition
whose complete pooled-root Schur weight excludes zero is the certificate.
Aggregation cannot create a nonzero total from zero physical edges, so this is
an existence certificate even though the old raw fibre does not retain a
single background mask.

Stop rule:

- nonzero certificate: `TWO_INDEPENDENT_DEFECT_GAIN_REJECTED`; retain the
  surviving leading four-arm signed functional;
- no certificate: only then study whether row and column defects localize in
  separable annuli.

Required outputs are the selected transition class, Bell before/after and
`g16` before/after, positive state masses `P0/P1`, signed Schur jets `S0/S1`,
their total, and the rank-transition by source-orbit matrix with `S*1` and
`1^T*S`.

