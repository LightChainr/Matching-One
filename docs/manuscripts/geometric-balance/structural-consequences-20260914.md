# Structural consequences from the geometric-balance programme — 2026-09-14

Status: owner-authorized continuation from PR #739 head `907a9d94`. This note consolidates deterministic consequences, consequences that additionally use the existing author-level #739 probability arguments, and clearly marked conjectural interfaces. It does **not** upgrade the parent manuscript to independent referee acceptance or certify literature priority.

The point of this note is to reduce, not enlarge, the open claim surface. Several questions that were being treated as separate analyses collapse once the finite 4/8 topology, the component-intensity point process, and the directional mass are used together.

## 1. Persistent 4/8 birth reflection

Let `v_1,...,v_N` be any strict ordering of the sites of an honest square-cell torus and let `S_k={v_1,...,v_k}`. Let `K_1^4,K_2^4` be the first indices at which the NN occupied ambient rank reaches one and two. On the matching graph occupy the same sites in reverse order and define `K_1^8,K_2^8` analogously.

The configurationwise digital-Alexander identity

\[
r_4(S)+r_8(S^c)=2
\]

implies the pathwise identities

\[
\boxed{K_1^8=N+1-K_2^4,\qquad K_2^8=N+1-K_1^4.}
\]

Proof: put `R_j=V\setminus S_{N-j}`. Then `r_8(R_j)=2-r_4(S_{N-j})`. Immediately before `K_2^4` the NN rank is at most one and at `K_2^4` it is two, so the reverse matching process first reaches rank at least one at `N+1-K_2^4`. The other identity is the rank-zero/positive version. A direct rank jump `0 -> 2` makes the two identities coincide rather than invalidating them.

For continuous labels `U_v` and reflected labels `V_v=1-U_v`, this becomes

\[
\boxed{T_1^8(V)=1-T_2^4(U),\qquad T_2^8(V)=1-T_1^4(U).}
\]

At every complementary rank-one pair the black and white ambient homology lines are the same. The existing digital-Alexander proof gives the image relation `C=A^perp`; in the two-dimensional symplectic space `H_1(T^2;Q)`, a one-dimensional isotropic line equals its symplectic orthogonal.

Finite controls committed with this note exhaust all 512 configurations and all `9!` site orders on `3 x 3`, and all 65,536 configurations plus 20,000 fixed-seed random site orders on `4 x 4`. They test rank duality, the winding-component count identity, rank-one line matching, and the birth-index reflection. These controls are not the proof.

## 2. Dual-even and dual-odd birth coordinates

Define, configuration by configuration,

\[
G=T_2-T_1,\qquad C=\frac{T_1+T_2-1}{2}.
\]

Under the 4/8 reflection `(T_1,T_2) -> (1-T_2,1-T_1)`, one has

\[
\boxed{G\mapsto G,\qquad C\mapsto-C.}
\]

Thus `G` is the digital-complement-even coordinate measuring separation of the two births, whereas `C` is the odd coordinate measuring their common displacement from `1/2`.

The two coordinates have exact integral representations which require no quantile fit:

\[
E T_1=\int_0^1P_0(p)\,dp,
\qquad
E T_2=\int_0^1[1-P_2(p)]\,dp,
\]

hence

\[
\boxed{E G=\int_0^1P_1(p)\,dp,}
\qquad
\boxed{E C=-\frac12\int_0^1M(p)\,dp.}
\]

If `K_1,K_2` are the two birth occupation indices in a uniform random permutation, then conditional order-statistic means give

\[
E G=\frac{E(K_2-K_1)}{N+1},
\qquad
E C=\frac{E[K_1+K_2-(N+1)]}{2(N+1)}.
\]

Existing rank-birth archives can therefore recover both coordinates without rescanning `p`.

If an independent fair selector chooses one of the two births, write it as a sign `S=+1/-1`. Then

\[
T=\frac12+C+\frac S2G,
\]

so exactly

\[
\boxed{\operatorname{Var}(T)=\operatorname{Var}(C)+\frac14E[G^2].}
\]

There is no `Cov(C,G)` term because the selector is independent and centered. This gives a clean decomposition of mixture broadness into configuration-to-configuration centre motion and within-configuration birth separation.

## 3. Exact reduction of same-parameter black/white winding counts

Let `W_4` be the number of black NN components with nonzero ambient homology and `W_8` the corresponding white matching count on the same labelled finite torus. The wrapping-component classification already used in the branch implies exactly one of

\[
(W_4,W_8)=(0,1),\quad (1,0),\quad (K,K),\ K\ge1.
\]

Consequently the full same-parameter joint probability generating function has the exact form

\[
\boxed{E[s^{W_4}t^{W_8}]
=P_2(p)s+P_0(p)t+P_1(p)H_p(st),}
\]

where `H_p(z)=E[z^K\mid r=1]`.

Thus a common critical-window model does **not** require an arbitrary two-dimensional joint count law. It is completely specified by the two rank endpoint probabilities and one positive-integer law in the rank-one sector. Independent or generically correlated bivariate Poisson ansatzes contain redundant degrees of freedom and usually violate the exact support.

Combining this finite support statement with the existing lower-window result `W_4 => Poi(lambda)` and `P(r=2)->0` gives the constrained same-parameter limit

\[
P[(W_4,W_8)=(0,1)]\to e^{-\lambda},
\]

\[
P[(W_4,W_8)=(k,k)]\to e^{-\lambda}\lambda^k/k!,\qquad k\ge1.
\]

The upper-window statement is the colour-reflected version. This is compatible with the exact obstruction to independent black/white Poisson counts at one common parameter, and distinct from the asymptotic independence of the two separated birth windows proved in the branch.

At the first-birth median `lambda=log 2`. Conditional on winding having occurred, the limiting probability of at least two black winding components is

\[
P(K\ge2\mid K\ge1)=0.3068528194\ldots.
\]

Therefore an event-conditioned winding configuration cannot be treated as containing a unique `birth cluster` by default. Component-Palm sampling is the correct microscopic protocol when a single component shape is required.

## 4. Alternating barriers and the reciprocal mean white span

Fix `p<p_c(NN)` on the infinite width-`w` cylinder. Order complete black NN essential components by their vertical topological order and anchor each at its minimum row using a fixed horizontal tie-break. Let `nu_w` be their stationary intensity per row. Let `B_i,B_{i+1}` be consecutive black essential components, `G_i` the difference of their anchor rows, and `W_i` the intervening white matching essential component.

The existing digital-Alexander construction gives more than equality of homology images. The facewise reduced embedded white graph is a CW 1-skeleton of the complement of a regular neighbourhood of the black graph. CW two-cells do not merge distinct 1-skeleton connected components, so white matching components correspond componentwise to connected complementary regions.

For consecutive black essential neighbourhoods the region between them is one essential annular complementary component, hence contains exactly one white essential component. With any fixed lattice regular-neighbourhood convention there is a universal lattice constant `C_0` such that

\[
G_i-L(B_i)-C_0\le L(W_i)\le G_i+L(B_{i+1})+C_0,
\]

and therefore

\[
|L(W_i)-G_i|\le L(B_i)+L(B_{i+1})+C_0.
\]

Publication-grade presentation should make the regular-neighbourhood thickness and one explicit admissible `C_0` visible, but no new probabilistic mechanism is needed for this comparison.

For any stationary ordered point process on the line, mass transport/Kac gives the Palm spacing identity

\[
\boxed{E^{Palm}G_i=1/\nu_w.}
\]

The existing fixed-subcritical component-activity estimate gives `E^{Palm}L(B_i)=O(w)`, whereas the established component intensity satisfies

\[
\nu_w=\exp[-\kappa(p)w+o(w)].
\]

Taking expectations in the deterministic span comparison therefore yields

\[
\boxed{\nu_w E[L(W_i)]\to1.}
\]

No renewal or Poisson hypothesis is needed for this mean identity.

The process Chen--Stein result can then be used only for the distributional upgrade. Taking torus height `m_w=lambda_w/nu_w` with `lambda_w -> infinity` but `log lambda_w=o(w)` preserves the vanishing absolute Poisson error. After scaling vertical coordinates by `nu_w`, the local barrier process tends to unit-rate Poisson. Together with `nu_w L(B_i)->0`, this gives the natural next lemma

\[
\nu_w L(W_i)\Rightarrow Exp(1).
\]

Sampling an interval by a uniform stationary vertical location instead of component Palm length-biases the exponential gap and gives `Gamma(2,1)`. The two observed laws are therefore different sampling protocols for one gap process, not separate shape mechanisms.

## 5. Marked-Poisson transport: do morphology once at component Palm

The anchor proof in `poisson-birth-windows.md` already works with local indicators. Any finite-valued mark that is measurable inside the same localization window can be attached to an anchor without changing the dependency graph. Decomposing marks into finitely many types gives a vector of locally dependent indicators and the same Chen--Stein estimates apply componentwise and jointly.

Therefore the efficient architecture for span/core/occupancy/boundary investigations is:

1. prove or measure the law of one complete winding component under component Palm;
2. localize/quantize the desired mark if necessary;
3. lift that marked Palm law to the long torus using the existing marked-Poisson machinery.

There is no reason to resimulate the whole exponentially long torus separately for each morphology observable.

This also clarifies the #762 operational core issue. The mark must be a precisely defined local functional of the full component. A `bridge-deleted winding connected part` and the union of nonzero-winding vertex-biconnected blocks are not generally the same object; contractible loops attached at an articulation provide a deterministic counterexample. The core definition must be fixed before any conditional shape result is interpreted.

## 6. Complementary Palm score identity

For a complete component `C` at site probability `p`, with `n(C)=|C|` and `b(C)` the number of **distinct external boundary sites**, the activity is

\[
q_p(C)=p^{n(C)}(1-p)^{b(C)}.
\]

Differentiating the log activity gives the exact component-Palm score

\[
S_p(C)=\frac{n(C)}p-\frac{b(C)}{1-p}.
\]

The exact complementary intensity identity `nu_w^4(p)=nu_w^8(1-p)` therefore gives, wherever finite differentiation is justified (and exactly at every fixed finite transfer representation),

\[
\boxed{
E_{4,p}\!\left[\frac n p-\frac{b_4}{1-p}\right]
+E_{8,1-p}\!\left[\frac n{1-p}-\frac{b_8}{p}\right]=0.}
\]

If the alternating-barrier result above is used on the white giant essential component at fixed black-subcritical `p`, then `E L_{white}\asymp1/nu_w` and connectedness with vertical step at most one implies `E n_{white}\ge E L_{white}\gg w`. The black score is only `O(w)` under its ordinary component Palm law. Consequently the white score must cancel internally to leading order, giving the testable prediction

\[
\boxed{\frac{E b_{white}}{E n_{white}}\to\frac{p}{1-p}.}
\]

This boundary/volume ratio is independent of the Brownian-range conjecture and is a useful checksum for future #762-type joint `(L,K,B)` work.

## 7. Convex structure of the loop--branch rate candidate

The current #758 candidate Palm-span rate is

\[
I_p(A)=\min_{0\le r\le A}
\{\tau_p(1,2r)-\kappa+\kappa(A-r)\},
\qquad \kappa=\tau_p(1,0)=\tau_p(0,1).
\]

This section analyzes the candidate only; it does not prove that the actual SITE Palm LDP equals it.

Put

\[
\phi(r)=\tau_p(1,2r),\qquad
 g(r)=\phi(r)-\kappa-\kappa r.
\]

Then

\[
I_p(A)=\kappa A+\min_{0\le r\le A}g(r).
\]

If the directional norm is strictly convex, `g` is strictly convex, `g'(0)=-kappa<0` by reflection symmetry, and norm bounds imply `g(r)->infinity`. Hence there is a unique finite minimizer `r_*>0`, independent of `A`, and

\[
\boxed{
I_p(A)=
\begin{cases}
\tau_p(1,2A)-\kappa,&0\le A\le r_*,\\
\kappa A+c_*,&A\ge r_*,
\end{cases}}
\]

with `c_*=g(r_*)`. When differentiable,

\[
2\,\partial_y\tau_p(1,2r_*)=\kappa.
\]

Square symmetry gives `0<r_*<1/2`: at `(1,1)`, exchange symmetry and Euler homogeneity give `partial_y tau(1,1)=tau(1,1)/2`, and strict convexity gives `tau(1,1)>kappa`.

Thus the candidate mechanism has a mandatory geometric transition: for small target span all extra height is supplied by bulging the closed essential loop; after `r_*` the optimal bulge saturates and every further unit of span is paid by a linear branch of slope `kappa`. If an eventual measured/rigorous large-`A` rate has slope strictly below `kappa`, the candidate topology is falsified by a cheaper network mechanism.

## 8. Diffusion coefficient from Wulff curvature: a three-way consistency relation

Assume the directional norm is twice differentiable at the horizontal direction and write

\[
\tau_p(1,y)=\kappa+\frac12a_p y^2+O(y^4),
\qquad a_p=\partial_{yy}\tau_p(1,0)>0.
\]

On the first branch of the candidate rate,

\[
I_p(A)=\tau_p(1,2A)-\kappa
=2a_pA^2+O(A^4).
\]

If the complete-component Brownian-bridge range picture is also valid,
`L/sqrt(D_p w) => R_BB`, and the Brownian bridge range tail satisfies

\[
-\log P(R_{BB}>r)=2r^2+O(\log r).
\]

In the moderate-deviation overlap `A->0`, `A sqrt(w)->infinity`, the Brownian prediction is

\[
-\frac1w\log P(L\gtrsim Aw)\sim\frac{2A^2}{D_p}.
\]

Matching the two descriptions forces

\[
\boxed{D_p^{-1}=\partial_{yy}\tau_p(1,0).}
\]

If the direction mass is written in angular form `tau(r cos theta,r sin theta)=r kappa(theta)`, reflection gives `kappa'(0)=0` and a direct expansion gives

\[
\boxed{D_p^{-1}=\kappa_p(0)+\kappa_p''(0).}
\]

This is not yet a theorem for the actual complete SITE component. It is a sharp three-way consistency condition connecting: (i) the #758 small-`A` loop-deformation LDP, (ii) the closed Brownian/HK transverse limit, and (iii) the directional mass/Wulff curvature targeted by #766. Proving any two in an overlapping regime determines and tests the third.

## 9. Directional competition in exponentially elongated tori

Let `Lambda` have shortest Euclidean period `u`, `|u|=ell`, area `N`, and transverse height `h=N/ell`. For any nonparallel `v in Lambda`, the determinant condition gives `|det(u,v)|>=N`, so `|v|>=h`.

For a square-symmetric norm `tau_p`, with `kappa=tau_p(1,0)`, convexity and symmetry give

\[
\kappa\|z\|_\infty\le\tau_p(z)\le\kappa\|z\|_1.
\]

Therefore

\[
\frac{\tau_p(v)}{\tau_p(u)}\ge\frac{h}{2\ell}.
\]

Whenever `h/ell -> infinity` (in particular under fixed positive `log h/ell`), all nonparallel homology classes are automatically separated in correlation cost. The only parallel primitive minimizers are `+/-u`. This removes a generic multi-direction competition from the default #765 proof architecture; remaining work is the actual SITE directional connection/seam estimate and parameter inversion.

## 10. Claim boundary

The following are finite/deterministic consequences once the already-committed digital-Alexander and wrapping-component classifications are accepted: persistent birth reflection, the `(rank,K)` support reduction, the dual-even/odd integral identities, and the convex algebra of the #758 candidate functional.

The reciprocal mean white span additionally uses the existing #739 author-level fixed-subcritical component-intensity exponent and component-volume control. The exponential white-span law additionally uses the process Poisson/Palm passage. The `D^{-1}` curvature relation is a consistency implication until the SITE LDP/Brownian overlap is proved. None of these statements identifies a continuum field or changes the original-U source contract.
