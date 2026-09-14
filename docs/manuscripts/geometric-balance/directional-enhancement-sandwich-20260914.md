# Direction-uniform enhancement sandwich for the full inverse-correlation norm

2026-09-14.  Consequence of the two-terminal pivotal conversion in `matching-enhancement-mass-gap-20260914.md`.  The axial strict inequality is already closed there.  This note records the stronger direction-uniform comparison that does **not** require a separate angle-by-angle computation.

## 1. Endpoint-independent pivotal conversion

The bounded surgery in Lemma 4.1 of the mass-gap note uses only:

1. a two-terminal connectivity event;
2. an original pivotal vertex lying on a shortest/induced open path joining the two terminals;
3. the Balister--Bollobas--Riordan square-lattice two-arm rerouting in a fixed ball;
4. insertion of the same facial-site gadget.

Nothing in the local map uses that the target is `n e1`.  Therefore, for every compact

\[
I\Subset(0,p_c(G8)),
\]

there exists `delta_I>0` and a finite endpoint cutoff such that for **every sufficiently distant original vertex** `x in Z^2` and every `p in I`,

\[
\boxed{
P_p^{G8}(0\leftrightarrow x)
\ge
P_{p+\delta_I}^{G4}(0\leftrightarrow x).}                      \tag{1.1}
\]

The same `delta_I` works for all directions.  Endpoint-near pivotal cases are still only finitely many local types and are absorbed into the same finite-energy constant.

This is stronger than comparing a finite set of axial/diagonal cylinder masses.

## 2. Full directional norm comparison

Let `tau_{G,p}` denote the subcritical inverse-correlation norm, so along any sequence `x_n/|x_n| -> e`,

\[
-\frac1{|x_n|}\log P_p^G(0\leftrightarrow x_n)
\to \tau_{G,p}(e).
\]

Taking logarithmic rates in (1.1) gives, uniformly in direction,

\[
\boxed{
\tau_{8,p}(e)
\le
\tau_{4,p+\delta_I}(e),
\qquad p\in I,\ e\in S^1.}                                    \tag{2.1}
\]

Equivalently for homogeneous vectors `x`,

\[
\tau_{8,p}(x)\le\tau_{4,p+\delta_I}(x).                        \tag{2.2}
\]

The existing graph-inclusion comparison only gives

\[
\tau_{8,p}(x)\le\tau_{4,p}(x).
\]

Equation (2.2) is a genuine strict-parameter improvement: one full unit of local matching enhancement beats a positive ordinary-site sprinkling uniformly over direction.

## 3. What is still needed for a same-p strict directional inequality

If one has strict parameter monotonicity of the NN norm in every direction,

\[
\tau_{4,p+\delta}(e)<\tau_{4,p}(e),                            \tag{3.1}
\]

then (2.1) immediately yields

\[
\boxed{\tau_{8,p}(e)<\tau_{4,p}(e)}                            \tag{3.2}
\]

uniformly on compact `p` intervals and all directions.

The present #739 branch proves a quantitative strict parameter inequality for the axial mass `kappa_4`; that is enough for the strict centre theorem already recorded.  This note does **not** silently promote (3.1) to every direction without either:

- a direct directional analogue of the branch's Friedgut--Kalai argument, or
- a precise site-percolation theorem giving strict `p` monotonicity of the full norm.

The non-strict enhancement sandwich (2.1) itself is already rigorous at the same author-proof level as the two-terminal pivotal map and is useful without (3.1).

## 4. Interface to #765

For a tilted exponential torus with shortest direction `e`, the centre equation should use `tau_{G,p}(e)`, not the axial mass.  Equation (2.1) supplies a model comparison:

\[
\tau_{8,p}(e)=d
\quad\Longrightarrow\quad
\tau_{4,p+\delta_I}(e)\ge d.                                  \tag{4.1}
\]

Together with monotonicity in `p`, this constrains the relative black/matching directional centre locations before any numerical directional mass is estimated.

Combined with the deterministic lattice-vector separation already proved in `structural-consequences-20260914.md`, the #765 architecture becomes:

1. shortest Euclidean period selects the only relevant homology direction in exponential elongation;
2. `tau_{G,p}(e)` sets the birth centre in that direction;
3. the matching enhancement sandwich compares the two graph norms uniformly over `e`.

The remaining hard step is then the actual SITE directional seam/connection estimate, not multi-direction competition.

## 5. Interface to #766

The first-exit exponential-moment domain in #766 is designed to give rigorous finite inner approximations to the Wulff/correlation body.  Equation (2.2) gives an independent inclusion check for any such certificates.

Write the polar/Wulff body schematically as

\[
K_{G,p}=\{t:t\cdot x\le\tau_{G,p}(x)\ \forall x\}.
\]

From (2.2),

\[
\boxed{K_{8,p}\subseteq K_{4,p+\delta_I}.}                    \tag{5.1}
\]

Thus any certified first-exit inner body for `K_{8,p}` that exits a rigorous outer body for `K_{4,p+delta_I}` would signal a normalization/implementation error.  Conversely, a directional certificate need not rediscover the graph-enhancement ordering numerically.

## 6. Possible stronger theorem

The most useful next theoretical closure is a direction-uniform strict parameter inequality

\[
\tau_{4,p}(e)-\tau_{4,q}(e)
\ge c_I(q-p)\tau_{4,q}(e),
\qquad p<q,\ e\in S^1,                                        \tag{6.1}
\]

on compact subcritical `I`.  This is the natural directional analogue of the axial inequality already derived from the birth/Friedgut--Kalai argument.  If (6.1) is proved, then (2.1) upgrades immediately to a compact-uniform strict graph gap

\[
\tau_{4,p}(e)-\tau_{8,p}(e)
\ge c_I\delta_I\tau_{4,p+\delta_I}(e)>0.                       \tag{6.2}
\]

Such a result would simultaneously strengthen #765 and provide a clean consistency constraint for #766's numerical Wulff-body certificates.
