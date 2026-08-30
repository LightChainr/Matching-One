# Radon--Nikodym anatomy of the remaining ULC gap

Let `mu_k` be uniform measure on the fixed-line layer `F_k`, and let `nu_k^up`
and `nu_(k+1)^down` be the two marginals of a uniform internal edge. If
`h_x=x/(N-k)` and `h_beta=b/k`, then the exact density changes are

`d nu_k^up / d mu_k = u/E_k[u] = (1-h_x)/(1-xi_k)`,

`d nu_(k+1)^down / d mu_(k+1) = d/E_(k+1)[d]`

`                              = (1-h_beta)/(1-beta_(k+1))`.

The lower marginal necessarily underweights high exit hazard. Expanding both
changes of measure turns the previous opaque degree bias into

`xi_(k+1)-xi_k = Delta_edge`

`  - Var_k(h_x)/(1-xi_k)`

`  + Cov_(k+1)(h_beta,h_x)/(1-beta_(k+1))`.

This proves the following conditional lemma:

> If `Delta_edge >= Var_k(h_x)/(1-xi_k)` and the upper-layer aggregate
> `Cov(h_beta,h_x)` is nonnegative, then the uniform exit hazard is
> nondecreasing.

Together with the dual-hazard lemma of `58768f0`, these two inequalities on
both complementary carriers imply fixed-line ULC.

Both hypotheses pass every one of the existing 984 carrier-layer pairs. The
degree bias is negative in 484 pairs, but its largest exact fraction of edge
slack is only `5/14`. Eight equivalent `N=12` matching realizations attain that
ratio, all at lower layer 6 with variance penalty `1/21`, edge slack `2/15`,
and uniform increment `3/35`. None saturates the desired inequality. The 78
zero uniform increments are all trivial: both edge slack and degree bias vanish.

Three tempting ways to prove the conditional hypotheses are nevertheless too
strong:

1. **First-order stochastic transport fails.** At matching `diag(3,4)`, line
   `(0,1)`, layer `4->5`, the lower uniform exit-hazard tail at threshold `3/4`
   is `4/19`, while the upper edge-weighted tail is `1/5`. The tail difference
   is `-1/95`, although the mean inequality remains positive.
2. **Pointwise birth/exit comonotonicity fails.** At primal `diag(2,4)`, line
   `(1,0)`, layer 5, masks 87 and 91 have `(b,x)=(0,2)` and `(2,1)`. The layer
   covariance is nevertheless the positive aggregate `1/90`.
3. **Worst-case Cauchy control fails.** Replacing the actual association by its
   most negative Cauchy value passes 980/984 pairs but fails four `N=12`
   realizations. In the first, the squared available margin is `4/103041`,
   smaller than the Cauchy variance product `50545/1179716409`; the actual
   positive covariance restores a uniform increment `4/321`.

Thus a general proof, if true, must be genuinely aggregate and mean-level. The
two promising statements are now explicit: edge slack dominates the lower
variance penalty, and birth/exit fragilities have nonnegative aggregate
covariance. Neither can be replaced by a full monotone coupling or a pointwise
ordering.
