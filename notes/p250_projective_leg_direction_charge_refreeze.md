# P250 direction tomography protocol repair

The first strict reveal rejected every frozen character.  That result exposed a
specific protocol defect rather than selecting a character: the strict model
silently set the leading root's axis-difference amplitude to zero.  The prior
rank-two analysis was fitted only to the axis-average T row and never supplied
that restriction.

Preserve the strict-null output and repair only this nuisance condition.  At
the same two frozen roots, fit both A-row amplitudes on `d1..4`.  Leave the
leading A amplitude free in each hand-charge channel.  Test the candidate q
only through the second-mode relation

\[
A_{2,c}=\frac{1-q}{1+q}T_{2,c}.
\]

The primary statistic is the joint four-complex-vector residual with
delete-one covariance.  Under each q, refit only the leading A amplitude on
`d1..4` and predict A at held-out `d5`.  A character closes only if both scores
have `p>=0.01`.

The candidate and exact phase alphabets remain unchanged: direction-invariant
`q=1`, spatial C4 `q=+/-i`, and Z5 `q=exp(2 pi i j/5)`.  Per-channel unrestricted
q estimates are descriptive diagnostics, not extra candidates.  This repair
is transparently post-strict-null; it removes an unjustified leading-mode
constraint but does not reopen rank, roots, or phase selection.
