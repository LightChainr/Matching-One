# P250 direction versus internal-charge tomography freeze

The N505 archive already contains both axis-average
`T=(X+Y)/2` and axis-difference `A=(X-Y)/2` at every `d=1..5`, with
delete-one batches.  Therefore `X=T+A` and `Y=T-A` are exact reconstructions;
no new simulation or rank selection is needed.

Keep the two roots selected at `c28fd7d`.  The leading, nearly real root is
frozen as direction invariant.  For the weaker second mode, write its y-axis
amplitude as `A_y=q A_x`.  Its observable axis-difference/axis-average
amplitude ratio is then

\[
\frac{A_x-A_y}{A_x+A_y}=\frac{1-q}{1+q}.
\]

The pre-reveal candidates are:

- direction-invariant/internal: `q=1`;
- spatial quarter-turn: `q=-i,+i`;
- an explicit one-step Z5 alphabet comparator:
  `q=exp(2 pi i j/5)`, `j=1,2,3,4`.

The exact phase alphabets are frozen separately.  A spatial C4 step has phase
`+/-pi/2`, while one Z5 deck-generator step has phase `2 pi j/5`.  In
particular, the prior second-root point phase near `-pi/2` is not in the Z5
one-step alphabet.  Root phase alone is not an identity: a compound or
non-generator transfer step can evade that alphabet.  Conversely, an
internal-only interpretation that keeps the observed non-Z5 root phase must
say explicitly that one separation step is not the deck generator.

Score the reconstructed complex axis-difference rows jointly across all four
hand-charge channels and `d=1..5`, with the complete delete-one covariance.
Refit the frozen rank-two recurrence in each deletion.  Report `d=5`
separately as a diagnostic, plus the distinct `Y=conjugate(X)` comparator.
Use `alpha=0.01`; never choose a direction character merely because its exact
phase is closest to the noisy root phase.

This test distinguishes a displacement-direction character from a
direction-invariant state.  It does not rename the state, reopen the rank
decision, or revisit cubic/OPE closure.
