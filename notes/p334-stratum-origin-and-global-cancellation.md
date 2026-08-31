# A complete orientation contrast cancels the arbitrary stratum origin

The full-birth readout `9c495ab1` permits an exact separation between a
checkpoint-stratum contribution and the complete topology contrast. The
large R1-membership step seen in the previous integral is not a new noise
source that must survive summing the other strata.

## The observable origin enters each stratum, but not their sum

Let A_i be either the original or safe whole-pair hybrid integrated A_top
readout in orientation i. Write I_ir=1{rank_i(k0)=r}, r=0,1,2, and let
delta be the fixed first-minus-second cos(4theta) normalization. Define

\[
 G_r=\frac{I_{fr}A_f-I_{sr}A_s}{\delta},\qquad
 V_r=\frac{I_{fr}-I_{sr}}{\delta}.
\]

Subtract a common fixed origin alpha from both orientations *before*
stratifying. Then

\[
 \widetilde G_r=G_r-\alpha V_r,\qquad
 \sum_r V_r=0,\qquad
 \boxed{\sum_r\widetilde G_r=\sum_rG_r
        =\frac{A_f-A_s}{\delta}.}
\]

This is a per-counter identity, not just equality of expectations. The
complete means, every batch contrast, covariance and suffix-noise residual
are unchanged. Individual stratum means and variances generally do change:

\[
 \operatorname{Var}(\widetilde G_r)
 =\operatorname{Var}(G_r)-2\alpha\operatorname{Cov}(G_r,V_r)
  +\alpha^2\operatorname{Var}(V_r).
\]

There is no general assertion that every individual stratum variance falls.
The invariant object is the complete contrast. One must use the same alpha
in the two orientations; independently fitting their origins would change
the target. The original uncentered stratum observable is still well-defined,
but its noise fraction cannot be promoted to an origin-independent global
mechanism.

## A fixed physical origin exposes the centered birth clock

For the original integrated observable,

\[
 A_i=1-\frac{K_{1,i}+K_{2,i}}{N+1}.
\]

Use the already declared p_ref=0.59274605079 and the non-fitted common
origin alpha=1-2p_ref. The centered observable is exactly

\[
 A_i-\alpha=-\frac{(K_{1,i}-(N+1)p_{ref})
                         +(K_{2,i}-(N+1)p_{ref})}{N+1}.
\]

Thus the centered stratum readout allocates the birth-clock displacement,
without allocating the large common baseline alpha. The same equation holds
with the conditional second-birth mean when the safe global gate is used.
This does not improve the total estimate by itself: its value and variance
are literally unchanged. It identifies which apparent stratum fluctuations
were introduced by distributing a common constant across random membership.

The interpretation is complementary to the exact C/W relation. With
C=(K1+K2)/2 and W=K2-K1, integrated A_top is1-2C/(N+1), while the rank-one
plateau integral is W/(N+1). Completing both birth terms removes W from the
first expression. A lifetime-rich contribution inside one stratum cannot
alone establish a lifetime mechanism for the complete integral.

## The covariance cancellation is essential

Since the three strata are mutually exclusive in one orientation,
Cov(I_ir A_i,I_it A_i)=-E[I_ir A_i]E[I_it A_i] for r!=t. For an orientation
contrast, additional cross-orientation terms use the same paired permutation.
All must be retained. The variance of the sum is not the sum of marginal
stratum variances, and an apparently large contribution can cancel through
these off-diagonal terms.

The same bookkeeping applies to any exhaustive source partition: shifting
the common observable origin reallocates alpha times source probabilities,
while leaving the completed orientation contrast invariant. This does not
invalidate a physically specified marked observable. It explains why its
source assignment must name the observer origin and why related decompositions
are not extra independent evidence.

Scientific card: exact finite-population identity, common observer origin,
three checkpoint ranks and the fixed same-permutation orientation pair.
The concrete centering is fixed by the existing p_ref, not chosen to optimize
the new data. Source is the same original e81dd59 paths completed at9c495ab1;
conditional values reuse0d1e586d. No new simulation, conditional solve or
verification suite. The numerical joint-stratum readout is a separate result,
not assumed by this derivation.
