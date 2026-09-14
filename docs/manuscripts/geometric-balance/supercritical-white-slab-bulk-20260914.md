# Bulk law inside the giant complementary white component

2026-09-14.  Continuation of the alternating-barrier / Poisson-tessellation analysis.  This note contains one exact infinite-volume identity and a concrete conditional law for the huge white matching component between rare black winding barriers.

## 1. An exact infinite-cluster boundary-density identity

Consider independent site percolation with open probability `q` on any locally finite graph, in a regime with a positive probability of an infinite open cluster.  Let

\[
\theta(q)=P_q(0\hbox{ belongs to an infinite open cluster})
\]

and

\[
\beta(q)=P_q(0\hbox{ is closed and has an open neighbour connected to infinity}).
\]

Then exactly

\[
\boxed{\beta(q)=\frac{1-q}{q}\,\theta(q).}
\]

**Proof.**  Freeze the states of all sites except the origin.  If the origin is open and lies in an infinite cluster, delete the origin.  The infinite cluster minus one finite-degree vertex has at least one infinite connected component: otherwise it would be a finite union of finite neighbour-components plus the origin and hence finite.  Thus in the frozen outside configuration some neighbour of the origin is connected to infinity without using the origin.  Conversely, if such a neighbour exists, opening the origin puts it in an infinite cluster.  Therefore the two events have the **same outside configurations** and differ only in whether the origin is open or closed.  Their product-measure weights have ratio `(1-q)/q`.  End of proof.

No uniqueness or one-endedness of the infinite cluster is required for this identity.

For the white matching graph at white density `q=1-p`, write `theta_8(q)` for the infinite-cluster density.  The density of black sites forming the external vertex boundary of that infinite white cluster is therefore

\[
\boxed{\beta_8(q)=\frac{p}{1-p}\theta_8(q).}
\]

This independently explains the boundary/volume ratio forced by the exact complementary component-Palm score in `structural-consequences-20260914.md`.

## 2. Geometry supplied by the black-barrier process

Fix black NN density `p<p_c(G4)` and white matching density `q=1-p>p_c(G8)`.  Let `nu_w=nu_w^4(p)` be the black winding-component row intensity.  The preceding notes give:

- black essential component thickness is `O(w)` under component Palm;
- black anchors, after vertical rescaling by `nu_w`, converge locally to a unit-rate Poisson process;
- between consecutive black essential components there is exactly one white matching essential component;
- its span `L_w` satisfies
  \[
  \nu_wL_w\Rightarrow Exp(1),\qquad \nu_wE L_w\to1.
  \]

Since `nu_w` is exponentially small in `w`, the intervening white component occupies a slab whose vertical length is exponentially larger than both the width and the subcritical black boundary layers.

## 3. Bulk-density conjecture with an exact target

Let `N_w` be the number of white sites in the intervening white essential component and `B_w` the number of **distinct black external boundary sites** of that component.

The natural supercritical slab statement is

\[
\boxed{
\frac{N_w}{wL_w}\xrightarrow{P}\theta_8(q),
\qquad
\frac{B_w}{wL_w}\xrightarrow{P}\beta_8(q)
=\frac{p}{1-p}\theta_8(q).}
\]

This is not yet proved here.  Its target constants, however, are not free amplitudes: they are the ordinary infinite-volume supercritical cluster density and the exact one-site-flip boundary density above.

A proof should use a bulk/boundary decomposition.  Remove from the white interval boundary layers whose thickness grows faster than the supercritical mixing/correlation scale but remains `o(L_w)`.  In the remaining slab, local events should be asymptotically unaffected by the two rare black barriers and by the conditioning that no additional black winding barrier occurs.  Spatial ergodicity/mixing then gives the two densities above.  The already established rare-event Poisson description should quantify why the no-extra-barrier conditioning has vanishing effect on a fixed or slowly growing bulk window.

This is a cleaner target than attempting to infer `N_w` or `B_w` from small-width component tables.

## 4. Joint macroscopic law if the bulk-density statement holds

Combine the bulk law with the exponential span law.  If `E~Exp(1)`, then jointly

\[
\boxed{
\left(
\nu_wL_w,
\frac{\nu_wN_w}{w\theta_8(q)},
\frac{\nu_wB_w}{w\beta_8(q)}
\right)
\Rightarrow(E,E,E).}
\]

Equivalently,

\[
\frac{B_w}{N_w}\xrightarrow{P}\frac{p}{1-p}.
\]

Thus almost all macroscopic randomness of `(L,N,B)` on the **white giant-component side** comes from the Poisson gap length.  Internal density fluctuations are lower order.

This predicts asymptotic correlations

\[
Corr(L_w,N_w)\to1,
\qquad
Corr(L_w,B_w)\to1,
\qquad
Corr(N_w,B_w)\to1,
\]

provided the second moments are uniformly integrable at the stated scaling.

That behaviour is deliberately different from the fixed-subcritical **black** complete-component conjecture `J` in `span-resolvent-frontier.md`, where Brownian transverse range and additive occupancy/boundary marks may become asymptotically independent after their own centering.  The complementary pair therefore supplies an internal positive/negative control:

- rare black component: width-`sqrt(w)` diffusive shape with `O(w)` additive marks;
- huge white complementary component: length-`1/nu_w` slab whose additive marks are dominated by the same exponential gap.

## 5. A sharper fluctuation conjecture

Condition on the white interval length `L_w`.  Standard supercritical mixing suggests

\[
N_w-\theta_8(q)wL_w=O_P(\sqrt{wL_w}),
\]

with a similar joint CLT for `B_w`.  Since `L_w` is of order `1/nu_w`,

\[
\frac{\sqrt{wL_w}}{wL_w}
=O_P\left(\sqrt{\frac{\nu_w}{w}}\right),
\]

which is exponentially smaller than the order-one randomness of `nu_wL_w`.

If a slab CLT is proved, then after subtracting the random gap-length contribution one should see a two-dimensional Gaussian bulk fluctuation:

\[
\frac1{\sqrt{wL_w}}
\begin{pmatrix}
N_w-\theta_8wL_w\\
B_w-\beta_8wL_w
\end{pmatrix}
\Longrightarrow N(0,\Sigma_{bulk}(q)).
\]

This is a secondary target; the law of large numbers already supplies the main geometric test.

## 6. Practical consequence for morphology work

A future #762-style joint `(L,K,B)` measurement on the complementary white side has parameter-free checks before any Brownian or OZ interpretation:

\[
N/(wL)\to\theta_8(1-p),
\qquad
B/N\to p/(1-p),
\qquad
\nu_wL\Rightarrow Exp(1).
\]

The second target is especially cheap because it needs no numerical value of `theta_8`.

Failure of `B/N -> p/(1-p)` at widths where black thickness is already negligible would localize a definition/sampling problem (wrong external-boundary convention, wrong component Palm, or white component not being the intended complementary giant) before it is interpreted as a failure of a continuum shape theory.

## 7. Claim boundary

The one-site-flip identity for `beta/theta` is exact.  The black Poisson-tessellation and reciprocal-span statements use the author-level #739 process theorem.  The white bulk law and conditional slab CLT are conjectural proof targets.  No new sampling is required to formulate them, and no unknown amplitude is introduced.
