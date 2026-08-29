# Q=1 four-leg spin-4 competitor: exact preflight

The generic-loop spectrum and the physical Potts state space are separate gates. This patch makes that distinction executable.

Using the published `V_(r,s)` weights at `beta^2=2/3`, exact rational arithmetic gives for `V_(2,2)`:

`(Delta, DeltaBar) = (1/8, 33/8)`, `x=17/4`, spin `-4`, and four legs.

Changing `s` to `-s` reverses chirality while preserving `x`. The thermal `Q4 epsilon` candidate has `x=21/4`, so the dimensions differ exactly by one. Consequently, in a continuum two-field mixture with a common observable normalization, the Q4-to-four-leg amplitude ratio changes by one inverse length under dilation, or `Q^-1/2` under a Gaussian area multiplier `Q`.

That last statement is deliberately not a score for PR #247's normalized shell coordinate. The issue itself requires the lattice-to-radial normalization first. The generated result therefore keeps four explicit unresolved gates: the Q=1 Potts multiplicity, the global matching matrix element, the local character-weighted pivotal overlap, and the normalized-shell transfer law.

This narrows the next analysis: `x=17/4` is a mandatory live adversary, but its mere presence in the generic loop spectrum is not evidence that either repository observable couples to it.
