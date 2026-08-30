# P250 minimal charged state-dimension freeze

The cross-scale scalar shapes failed, so the next question is state dimension,
not another exponent.  This retrospective reanalysis uses only the archived
N325 and N505 complex pair batches and their complete covariance.

For every hand-charge channel, freeze a shared complex recurrence

\[
T_c(d)=a_1T_c(d-1)+\cdots+a_KT_c(d-K).
\]

The coefficients, hence the complex transfer eigenvalues, are common to all
four channels; each channel retains its own complex amplitudes.  N505
`d=1..4` determine the recurrence and `d=5` is held out.  The score order is
rank one, rank two, then rank three only if necessary.

N325 is an independent geometry constraint, not extra target-fit data.  Rank
one is tested on its `d2,d3` recurrences and rank two on `d3`.  Two eigenvalue
transports are declared: unchanged lattice-unit eigenvalues, and complex-log
eigenvalues scaled as `1/L`.  Rank three cannot be identified from only three
source distances and must say so explicitly.

As a non-Prony comparator, freeze one minimal two-dimensional image kernel: the
axis-average of the nearest `3x3` period images of `|r|^-alpha`.  N505 d1--d4
fit its common alpha and channel amplitudes, d5 is held out; on N325 only d1--d2
fit amplitudes and d3 is held out.  This is a deliberately simple finite-image
model, not a claim that conformal symmetry fixes the full torus two-point
function.

The first rank whose N505 holdout and at least one N325 transport both have
`p>=0.01` is the minimal identifiable cross-geometry dimension.  All complex
phases and delete-one batch covariance are retained.  No simulation, cubic
row, OPE phase, or single-power rescue is allowed in this score.
