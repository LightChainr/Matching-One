# Scientific card: N170 flips the charged current and localizes the excess to scale curvature

- Mechanism changed: the exact N85-to-N170 angle flip is strongly resolved.
  `Delta_K_A(N170)=-0.0177169 +/- 0.0015556`, or `-11.389 sigma` against the
  scalar/geometry-blind zero contrast.
- Frozen H4 vector: the preregistered prediction was
  `(+0.0039355,-0.0058502)`; observed was
  `(+0.0064347,-0.0112822)`.  Full predictive score is `10.252/2` and the pair
  exceeds the frozen magnitude by `3.035` predictive SE.
- Localization: in the preregistered exact basis, H4 amplitude is
  `-0.0111115` versus frozen `-0.0061373` (`-3.035 SE`), while the orthogonal
  charged/projective scalar is `-0.0013810 +/- 0.0017648` (`-0.783 sigma`).
- Mechanism-space update: geometry controls the sign, but the H4-only radial
  amplitude is insufficient.  The residual is scale curvature along the same
  geometry covector, not a projective common scalar.
- Does not prove: a new exponent, asymptotic curvature form, or an H4/H8
  distinction.  This experiment deliberately performs no harmonic revote.
- Observer/sector/source/geometry: `K_A=d_eta log W_A`; exact angle-flip child
  N170 `(11+7i,13+i)`; 8M fresh samples per shape, 80 aligned batches and full
  covariance.
- Dependency: H4 vector and fit covariance frozen in `cf1bdf8` from the
  N65/N85-only model `186d72a`; N170 is a new independent seed/counter block.
- Next discriminator: model the same-lineage H4 amplitude curvature with the
  already available N85/N170 pair before choosing another geometry.  Keep the
  projective scalar as a zero control.
