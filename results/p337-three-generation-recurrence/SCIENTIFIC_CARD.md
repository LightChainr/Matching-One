# Scientific card: a two-mode cancellation explains the overshoot/return pattern

- Mechanism changed: after exact H4 sign alignment, the same-lineage amplitudes
  are `(-0.016051,-0.011111,-0.004857)` at N85/N170/N340.  Fixing the leading
  transfer to `lambda0=2^-13/8` makes the exact two-mode identity return
  `lambda1=0.2124 +/- 0.2364`.
- Interpretation: the point solution has `c0=-0.06889`, `c1=+0.05284`.  The
  correction is opposite in sign and faster than the leading term, so it
  cancels 76.7%, 50.2%, and 32.9% of the leading magnitude over the three
  observed generations.  This produces the N170 overshoot and N340 return
  without a scale-neutral charged state.
- Uncertainty: `lambda1` lies in `(0,1)` at the point estimate but is only
  `0.898 sigma` above zero; its 95% delta interval is `[-0.251,0.676]`.  N340
  contributes 79.9% of its propagated variance.  This is a useful mechanism
  coordinate, not a resolved second eigenvalue.
- N680 discriminator: the two-mode prediction is
  `A_H=-0.0018413 +/- 0.0009551`, versus `-0.0010272` for a globally fitted
  frozen single mode, `-0.0030378` for one free single transfer, and
  `-0.0089120` for scale-neutral.  The exact H4 child flips the pair negative.
- Model accounting: the two-mode model has three parameters for three data and
  therefore zero residual degrees of freedom.  Its scientific content is the
  frozen N680 prediction, not its interpolation score.
- Does not prove: a unique correction field, a continuum exponent, or a
  cross-geometry scaling law.  All algebra is confined to one Gaussian
  `1+i` lineage after exact covector alignment.
- Dependency: N85/N170/N340 are independent random blocks; each amplitude uses
  the full paired-orientation covariance, and the inter-generation covariance
  is exactly block diagonal by construction.
- Natural upgrade: an N680 block would distinguish the four frozen predictions
  directly; no fourth unrelated geometry is needed.
