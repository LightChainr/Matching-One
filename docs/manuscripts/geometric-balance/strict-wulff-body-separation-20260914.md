# Strict same-parameter separation of the NN and matching exponential-moment bodies

2026-09-14.  Convex-geometric consequence of the all-direction strict matching mass gap proved on this continuation branch.

## 1. Correlation norms and dual bodies

For `G=G4` (NN) or `G8` (matching), fixed subcritical `p`, let

\[
\tau_{G,p}(x)
\]

be the homogeneous inverse-correlation norm.  Define its dual/exponential-moment body

\[
K_{G,p}
=\{t\in\mathbb R^2:t\cdot x\le\tau_{G,p}(x)\ \forall x\}.    \tag{1.1}
\]

The support function is exactly

\[
\boxed{h_{K_{G,p}}(x)=\tau_{G,p}(x).}                          \tag{1.2}
\]

The first-exit theorem on this branch identifies the interior of `K_{G,p}` with the exponential-susceptibility domain.

## 2. Pointwise strict support-function gap

`varying-direction-exponential-centres-20260914.md` combines endpoint-direction-independent matching enhancement with all-direction strict `p`-monotonicity and proves

\[
\boxed{
\tau_{8,p}(e)<\tau_{4,p}(e)
\qquad
(e\in S^1,\ 0<p<p_c(G8)).}                                   \tag{2.1}
\]

Hence

\[
K_{8,p}\subsetneq K_{4,p}.                                    \tag{2.2}
\]

The inclusion is in the correct direction: matching connectivity is easier, so its inverse-correlation norm and therefore its dual support body are smaller.

## 3. Uniform Euclidean buffer on compact parameter intervals

Fix

\[
I\Subset(0,p_c(G8)).                                           \tag{3.1}
\]

The directional norms are continuous in `(p,e)` on `I x S^1` under the same standard subcritical norm regularity used in the directional centre theorem.  Therefore

\[
\Delta(p,e)=\tau_{4,p}(e)-\tau_{8,p}(e)                       \tag{3.2}
\]

is continuous and strictly positive on the compact set `I x S^1`.  Put

\[
\boxed{
\eta_I=\min_{p\in I,e\in S^1}\Delta(p,e)>0.}                  \tag{3.3}
\]

By homogeneity, for every `x in R^2`,

\[
\tau_{8,p}(x)+\eta_I|x|\le\tau_{4,p}(x).                     \tag{3.4}
\]

Let `B_2` be the Euclidean unit disk.  Since support functions add under Minkowski sum,

\[
h_{K_{8,p}+\eta_I B_2}(x)
=h_{K_{8,p}}(x)+\eta_I|x|.                                   \tag{3.5}
\]

Equations (1.2), (3.4) and the support-function characterization of convex-body inclusion give

\[
\boxed{
K_{8,p}+\eta_I B_2\subset K_{4,p}
\qquad(p\in I).}                                              \tag{3.6}
\]

This is stronger than strict inclusion: the smaller matching body can be thickened by a fixed positive Euclidean radius, uniformly over the whole compact parameter interval, and still remain inside the NN body.

## 4. Equivalent support and gauge statements

Equation (3.6) is equivalent to the uniform support gap

\[
\boxed{
h_{K_{4,p}}(e)-h_{K_{8,p}}(e)\ge\eta_I
\quad(e\in S^1,p\in I).}                                     \tag{4.1}
\]

If `B_{G,p}={x:tau_{G,p}(x)<=1}` is the primal correlation-norm unit ball, then the easier matching model has the larger primal ball,

\[
B_{4,p}\subsetneq B_{8,p}.                                    \tag{4.2}
\]

The dual-body buffer (3.6) is usually the cleaner quantitative statement for first-exit certification, because the finite certificates approximate `K_{G,p}` directly through exponential tilts.

## 5. Finite first-exit certification consequence

Let `C^{(G)}_S(p)={t:B^{(G)}_S(t;p)<1}` be a rigorous first-exit certified region.  The large-box exhaustion theorem says that for every compact subset of `K_{G,p}^circ`, sufficiently large boxes certify it.

Fix `I` and choose any `epsilon in (0,eta_I/4)`.  Then

\[
K_{8,p}+ (\eta_I-\epsilon)B_2
\Subset K_{4,p}                                               \tag{5.1}
\]

uniformly for `p in I` after an arbitrarily small inward shrink if one wants compact containment in the interior.

Thus a numerical/certified Wulff computation has a strong same-parameter cross-graph control:

1. certify a compact inner approximation to `K_{8,p}`;
2. thicken it by any radius safely below `eta_I` once a rigorous lower bound on the support gap is available;
3. the result must still lie inside the true NN body and should eventually be certifiable by large enough NN first-exit boxes.

Conversely, a claimed pair of certified bodies that violates `K8 subset K4` has an adjacency, orientation, support-function or outward-rounding error before it indicates physics.

## 6. Enhancement sandwich gives a second nested-body relation

The local enhancement argument actually gives, on a compact interval `I`, a positive site-density sprinkling `delta_I` such that

\[
\tau_{8,p}(x)\le\tau_{4,p+\delta_I}(x)                        \tag{6.1}
\]

for all directions.  Therefore

\[
\boxed{K_{8,p}\subset K_{4,p+\delta_I}.}                      \tag{6.2}
\]

Since `p -> tau_{4,p}` is strictly decreasing, the NN bodies shrink strictly as `p` increases, and

\[
K_{4,p+\delta_I}\subsetneq K_{4,p}.                           \tag{6.3}
\]

So the same-parameter strict gap can be viewed as the convex-body sandwich

\[
\boxed{
K_{8,p}\subset K_{4,p+\delta_I}\subsetneq K_{4,p}.}           \tag{6.4}
\]

This is the geometric version of “full matching enhancement beats a positive ordinary-site sprinkling.”

## 7. Dilute asymptotic interpretation

The dilute directional results on this branch show that the support gap is smallest, to leading order, near the coordinate axes.

Axially,

\[
\tau_{4,p}(e_1)-\tau_{8,p}(e_1)
=\log3+O(p)                                                    \tag{7.1}
\]

as `p downarrow 0`.

For a genuinely tilted fixed rational direction the coefficient of `log(1/p)` is already strictly larger for NN than matching, so the gap diverges as `p downarrow0`.

This suggests, and the explicit fixed-rational formulas verify directionwise, that the compact-interval buffer remains macroscopic in the dilute regime rather than collapsing there.

A uniform asymptotic statement `eta_p -> log3` would additionally require controlling the minimising direction uniformly as `p->0`; it is not asserted here.

## 8. Claim boundary

The Minkowski-buffer conclusion is an exact convex consequence of: (i) all-direction strict mass gap; and (ii) continuity of the directional norms on compact subcritical parameter-direction sets.  No differentiability in direction or in `p`, no OZ amplitude, and no numerical Wulff construction is required.