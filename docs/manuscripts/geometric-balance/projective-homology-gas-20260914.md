# Projective-homology gas: the exact finite state behind directional winding counts

2026-09-14.  This note refines the `(rank,K)` reduction by retaining the one piece of directional information that survives exactly at rank one: a projective homology line.  It is deterministic topology first, with a conjectural dilute-gas interpretation second.

## 1. One-colour exact classification

Let an occupied graph on an honest torus have connected components `C_j`.  For each component define

\[
A_j=\operatorname{im}[H_1(C_j;\mathbb Q)\to H_1(T^2;\mathbb Q)].
\]

The global ambient image is

\[
A=\operatorname{span}_j A_j.
\]

Because `H_1(T^2;Q)` is two-dimensional, exactly the following possibilities occur.

### Rank zero

`A=0`.  No component is essential.

### Rank one

`A` is one line `ell` in `P(H_1(T^2;Q))`.  Every essential component has nonzero image contained in `A`, hence

\[
\boxed{A_j=A=\ell\quad\hbox{for every essential component}.}
\]

Thus a rank-one configuration is described topologically by one projective rational slope `ell` and an integer count `K>=1` of parallel essential components.

### Rank two

At least one component must have two-dimensional image.  Indeed, if every essential component had rank one, then two components with distinct lines would be disjoint essential subsets.  Representatives of nonparallel homology classes on a torus have nonzero algebraic intersection and cannot be supported in disjoint components.  Hence all rank-one component images would have to be the same line, contradicting global rank two.

A rank-two component excludes every other essential component: any essential cycle disjoint from it has zero intersection with every class carried by the rank-two component, hence its ambient class is zero by nondegeneracy of the torus intersection form.  Therefore

\[
\boxed{r=2\quad\Longrightarrow\quad\text{there is exactly one essential component, of rank two}.}
\]

This recovers the single/cross component part of the wrapping classification without referring to a probabilistic model.

## 2. Complementary 4/8 classification with slope

For the black NN graph and the complementary white matching graph, digital Alexander gives

\[
A_{white}=A_{black}^{\perp}.
\]

Hence:

- black rank `0`: white rank `2`, with one white rank-two component;
- black rank `2`: white rank `0`, with one black rank-two component;
- black rank `1`: `A_black=ell` and, because every one-dimensional subspace of the two-dimensional symplectic torus homology is self-orthogonal,
  \[
  A_{white}=ell.
  \]
  The alternating complementary-region argument gives equal black/white essential-component counts `K`.

So the exact same-parameter two-colour state is

\[
\boxed{
\begin{cases}
D=-1:&(r_4,r_8)=(0,2),\ K=0,\\
D=0:&(r_4,r_8)=(1,1),\ (\ell,K),\ K\ge1,\\
D=+1:&(r_4,r_8)=(2,0),\ K=0,
\end{cases}}
\]

where `D=r_4-1=W_4-W_8`.

The previous `(rank,K)` reduction is the marginal obtained after forgetting `ell`.

## 3. Exact generating object

Let `chi` be any test function/character on projective rational homology lines.  Define the rank-one marked generating function

\[
H_p(z;\chi)
=E[z^K\chi(\ell)\mid r=1].
\]

Then every same-parameter topological statistic that depends only on essential-component count and slope is encoded by

\[
\boxed{
Z_p(h,z;\chi)
=P_2(p)e^h+P_0(p)e^{-h}+P_1(p)H_p(z;\chi).}
\]

At `z=1, chi=1`, the derivative in `h` at zero is exactly the matching observable

\[
\partial_hZ_p(0,1;1)=P_2-P_0=M(p).
\]

Thus the finite matching root is a **zero mean topological charge** condition.  The potentially extensive neutral rank-one gas `K` is orthogonal to that charge at the level of this exact decomposition; only the unpaired rank-zero/rank-two endpoint sectors set the sign of `M`.

## 4. Consequence for same-parameter count fluctuations

The exact identities may be written

\[
W_4=K+1_{D=+1},
\qquad
W_8=K+1_{D=-1},
\]

with `|D|<=1`.  Hence

\[
|W_4-W_8|\le1
\]

configuration by configuration.  At any fixed width and long height, any extensive CLT or LDP for one colour transfers to the other with the **same** fluctuation rate.  In particular, if `Var(K)` is of order height, then

\[
Corr(W_4,W_8)\to1.
\]

The two colours are therefore maximally correlated at the extensive same-parameter count scale, in deliberate contrast with the asymptotic independence of counts in the two **separated parameter windows**.

This is also why the joint pressure found earlier depends only on the sum of the two count sources at leading order: the topological charge is bounded and disappears after division by height.

## 5. Directional dilute-gas conjecture

For a primitive lattice period/slope `ell` represented by `u`, let the type fugacity be schematically

\[
\lambda_u(p)
\asymp \frac{N}{|u|}\,A(p,\hat u)|u|^{-\beta(p,\hat u)}e^{-\tau_p(u)}.
\]

The `homological-free-energy` variable in `research-frontier-20260914.md` is the logarithmic leading part of this expression.

When one type has much larger fugacity than all others, the rank-one slope should concentrate on that type; the deterministic direction-separation lemma proves this in the fixed-positive exponential elongation regime at the level of correlation cost.

When finitely many primitive slopes have comparable fugacity, the correct rank-one object is **not** a collection of independent Poisson species.  Rank one permits only one projective line.  A better picture is a projective hard-core gas:

1. the first winding event is a race among slope types;
2. after a slope `ell` wins, additional rank-one components may accumulate only parallel to `ell`;
3. the appearance of a nonparallel homology class forces the system into rank two and ends the rank-one sector.

This supplies a topological state space for the multi-direction version of #765 without inventing a bivariate/multivariate independent-Poisson approximation that the torus cannot support.

## 6. First-birth slope selection: a falsifiable leading rule

In a regime where all type fugacities are small and simultaneous nonparallel events are of higher order, the first rank birth should select slope `u` with probability approximately

\[
\boxed{
P(\ell=u\mid T_1\hbox{ occurs in the window})
\approx \frac{\lambda_u}{\sum_v\lambda_v}.}
\]

The first-birth void probability should depend at leading order on the sum of candidate fugacities, while the rank-two birth is sensitive to their incompatibility/intersection structure.  This is a conjectural asymptotic rule, not an exact finite formula.

A useful test geometry is a family with two intentionally near-degenerate primitive period classes.  The prediction is not merely a 50/50 split: the ratio should be set by the directional masses, transverse opportunity factors, and eventually the sewing amplitudes.  The rank-one plateau should be shorter than in a one-minimizer geometry because a nonparallel second type is an efficient route to rank two.

## 7. Relation to continuum/modular work

No continuum field identification is needed for this state.  On fixed-aspect tori, a scaling limit of the marked slope distribution would be a natural lattice object to compare with known homology/wrapping formulas.  On exponentially elongated tori, the deterministic cost separation collapses it to the shortest direction.  These are two limits of the same exact finite projective-homology state, rather than two unrelated observables.
