# Gaussian x annulus context rectangle

The semantic crosswalk first rejected a numerical splice between PR277 global
threshold derivatives and P253 local pivotal amplitudes.  The missing Gaussian
cells were therefore acquired with the P253 source/readout itself.  Only the
quotient-coordinate backend changed: cyclic Gaussian labels were replaced by
the existing arbitrary integer-period HNF/Smith backend.

This adapter is deliberately narrow.  It records the same integer pivotal and
landing-H4 sufficient statistics, uses the same fixed-`p` root toggle, and has
no generic observable API.  Primitive cyclic equivalence and nonprimitive
direct-coset equivalence are exact finite checks.  The R2 constructor also
enforces the exact direction-character gate: axis and diagonal orbits must both
survive, since either orbit alone aliases scalar and spin 4.

The production result does not favor a richer context model.  The unrestricted
best pair improves the best shared score by only 0.2162, with composite-null
bootstrap p=0.4002.  Scientifically this is a positive transport result for the
shared-generator description at current resolution, not an identification of
the nominal best lambda and not a general no-memory theorem.

