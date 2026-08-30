# Scientific card: the old N65 archive already resolves both explicit F3 charges

- Mechanism changed: the explicit A/B1 and D/B2 sources frozen in `539b629`
  require no new simulation.  The archived `(tau1,ell,tau2)` cells identify
  their state, birth and exit response with one existing covariance block.
- New result: both charged activations are precisely measured.  The frozen
  same-N orientation triplets give `12.153/3 df` for A and `1.509/3 df` for D;
  A timing, but not D, carries the visible orientation modulation.  The joint
  six-vector is `15.530/6 df`.
- Exact controls: A-D cross response is zero statewise;
  `dW_C/dp=J_C,birth-J_C,exit` holds to `4.1e-15`; F3 phase is fixed by
  `O_C(omega)=(omega-omega^2)W_C/2`.  Internal T-shear transport agrees to
  `4.2e-17`.
- Does not prove: large-N survival, continuum field identity, or an independent
  shear result.  The two N65 shapes are not an identity/T-source pair.
- Observer/sector/source/geometry: F3 projective rank-one plateau and birth/exit
  currents; C4 charge-2 A/B1 and D/B2; explicit `omega^q`; Gaussian N65 `8+i`
  and `7+4i` on the same 20k counter block.
- Dependency: archive `1714141`, flat-twist placement `a7cb19a`, source freeze
  `539b629`.
- Next discriminator: reuse the same six-vector contract at a larger existing
  archive or in the planned common rank-birth stream.  Preserve the full
  covariance; the A signal is a joint timing response, not a `W_A` amplitude.
