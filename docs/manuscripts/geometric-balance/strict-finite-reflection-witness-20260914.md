# A matching-only essential cycle makes finite reflection dominance strict

2026-09-14.  Deterministic witness completing the strictness clause in `finite-reflection-dominance-20260914.md` for every nondegenerate honest square-cell torus with sufficiently separated local lifts (in particular throughout the asymptotic regime of #739).

No percolation estimate is used: one explicit occupied set has positive product probability for every `0<p<1`.

## 1. Goal

Graph inclusion always gives

\[
r_8(\omega)\ge r_4(\omega).                                  \tag{1.1}
\]

To make the finite reflection inequalities strict it suffices to construct one configuration with

\[
\boxed{r_8(\omega)>r_4(\omega).}                              \tag{1.2}
\]

We will in fact construct

\[
r_4(\omega)=0,
\qquad
r_8(\omega)=1.                                                 \tag{1.3}
\]

## 2. Non-axis shortest period

Let `u=(a,b)` be any nonzero period with both coordinates nonzero.  Apply an exact square-lattice reflection/quarter-turn so that

\[
a\ge b>0.                                                      \tag{2.1}
\]

Consider the lifted vertex sequence

\[
(0,0),(1,1),\ldots,(b,b),(b+1,b),\ldots,(a,b)=u.              \tag{2.2}
\]

Project its distinct vertices modulo the period lattice, identifying only the final endpoint `u` with the starting vertex `0`.  For a shortest period in an honest quotient, no two other vertices in the monotone rectangle can be period translates: their Euclidean difference is strictly shorter than `|u|`.

In the matching graph, consecutive diagonal steps followed by horizontal steps form a closed path whose lift displacement is exactly `u`.  Hence the occupied set has nonzero ambient homology.

In the NN graph, every diagonal step is absent.  The diagonal vertices before `(b,b)` are isolated from one another in NN connectivity, while the horizontal tail from `(b,b)` back to the identified endpoint `u=0` is only a path.  There is no NN closed walk with nonzero gain.  Monotonicity of the lifted coordinates rules out an incidental NN chord; such a chord would require two listed vertices to differ by one axial step, which occurs only along the declared horizontal tail.

Thus (1.3) holds.

The same argument works for arbitrary sign patterns by reflection.

## 3. Axis shortest period

Suppose instead

\[
u=(a,0),                                                     \tag{3.1}
\]

with a nondegenerate circumference.  Use the lifted matching path

\[
(0,0)\to(1,1)\to(2,0)\to(3,0)\to\cdots\to(a,0)=u.           \tag{3.2}
\]

The first two edges are matching diagonals with displacements `(1,1)` and `(1,-1)`; the remaining edges are horizontal NN edges.  The total lift displacement is `(a,0)=u`, so after quotienting this is an essential matching cycle.

In the NN graph, the two diagonal connections are absent.  The occupied horizontal tail from `(2,0)` through `(a,0)=0` is a path with one missing link at the detour; `(1,1)` does not close it.  Hence the NN ambient rank is zero.

A vertical axis period is the quarter-turned version.

For the very shortest quotients one must retain the repository's existing lifted-parallel-edge convention or use the separate finite oracle.  In the #739 asymptotic honest-torus scope, where the systole tends to infinity, no local degeneracy enters.

## 4. Strict finite endpoint inequality

Fix `0<p<1`.  The explicit occupied set above, together with every site outside it declared vacant, is one finite configuration of strictly positive Bernoulli probability.  On it

\[
r_4=0,
\qquad r_8>0.                                                  \tag{4.1}
\]

Therefore the graph-inclusion event containment is strict:

\[
\boxed{P_p^{8}(r=0)<P_p^{4}(r=0).}                            \tag{4.2}
\]

Digital Alexander complement gives

\[
P_{1-p}^{4}(r=2)=P_p^{8}(r=0).                                \tag{4.3}
\]

Hence

\[
\boxed{P_2^{4}(1-p)<P_0^{4}(p),\qquad0<p<1.}                  \tag{4.4}
\]

Reflecting `p` gives

\[
P_2^{4}(p)<P_0^{4}(1-p).                                      \tag{4.5}
\]

Adding the two strict inequalities yields

\[
\boxed{M(p)+M(1-p)<0,\qquad0<p<1.}                            \tag{4.6}
\]

## 5. Strict finite root and quantile reflection

At `p=1/2`, equation (4.6) gives

\[
\boxed{M(1/2)<0.}                                             \tag{5.1}
\]

Since the finite matching function is strictly increasing with a unique zero,

\[
\boxed{p_\Lambda>1/2.}                                       \tag{5.2}
\]

Likewise the fair birth-mixture CDF obeys

\[
F(p)+F(1-p)<1                                                  \tag{5.3}
\]

throughout the interior, except that quantile inequalities can become equal at levels affected by finite jumps only under the usual generalized-inverse convention.  For the continuous-label birth law on a nondegenerate finite torus the distribution is continuous, and

\[
Q(u)+Q(1-u)>1                                                  \tag{5.4}
\]

whenever the two corresponding quantiles lie in the strict interior regime.

The universal non-strict statement remains the safe formulation for any tiny quotient outside the honest-cell scope.

## 6. Interpretation

The negativity of `M(1/2)` seen in the exact `L=3,4` censuses is therefore not an accidental small-size feature.  It is forced by a simple structural fact:

> the matching graph admits essential cycles on vertex sets whose NN graph has no essential cycle.

This finite topological asymmetry is the zero-scale precursor of the strict inverse-correlation mass gap `kappa_8(p)<kappa_4(p)` proved separately on this continuation branch.

The two theorems live at different levels:

- this note: one explicit finite configuration gives strict sign/root asymmetry;
- enhancement theorem: exponentially many long connections give a strict difference in the infinite-volume directional mass and hence strict separated-window centre asymmetry.
