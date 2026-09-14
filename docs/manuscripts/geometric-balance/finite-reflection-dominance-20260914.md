# Finite reflection dominance from graph inclusion plus persistent Alexander duality

2026-09-14.  Exact finite theorem for every honest torus.  No sharpness, RSW, mass, OZ or asymptotic limit is used.

The result says that the NN two-birth law is intrinsically shifted to the right of its reflection about `1/2`.  The strict large-`d` centre theorem on this branch is a quantitative asymptotic strengthening, not the source of the sign.

## 1. Rank monotonicity under matching enhancement

For the **same occupied vertex set** `omega`, the matching graph contains all NN edges and adds diagonals.  Therefore the ambient homology image can only grow:

\[
\boxed{r_8(\omega)\ge r_4(\omega).}                           \tag{1.1}
\]

Under the same continuous labels `U_v`, let

\[
T_1^G\le T_2^G
\]

be the first parameters at which the ambient rank in graph `G` reaches one and two.  Equation (1.1) gives pathwise

\[
\boxed{T_j^8(U)\le T_j^4(U),\qquad j=1,2.}                    \tag{1.2}
\]

## 2. Persistent duality turns graph inclusion into reflection dominance

The persistent digital-Alexander identity proved in `structural-consequences-20260914.md` gives, for reflected labels `V=1-U`,

\[
T_1^8(V)=1-T_2^4(U),
\qquad
T_2^8(V)=1-T_1^4(U).                                          \tag{2.1}
\]

Since `V` has the same iid uniform law as `U`,

\[
T_1^8\ \overset d=\ 1-T_2^4,
\qquad
T_2^8\ \overset d=\ 1-T_1^4.                                \tag{2.2}
\]

Combine (1.2) and (2.2).  In stochastic order,

\[
\boxed{1-T_2^4\preceq_{st}T_1^4,}                            \tag{2.3}
\]

\[
\boxed{1-T_1^4\preceq_{st}T_2^4.}                            \tag{2.4}
\]

These are exact finite inequalities for the NN birth pair.

## 3. Endpoint-sector form

Write

\[
f_1(p)=P(T_1\le p)=1-P_0(p),
\qquad
f_2(p)=P(T_2\le p)=P_2(p).                                   \tag{3.1}
\]

Equation (2.3) is equivalent to

\[
f_1(p)+f_2(1-p)\le1.                                         \tag{3.2}
\]

Substituting (3.1) gives the particularly simple topological inequality

\[
\boxed{P_2(1-p)\le P_0(p),\qquad0\le p\le1.}                  \tag{3.3}
\]

Reflecting `p` gives also

\[
P_2(p)\le P_0(1-p).                                           \tag{3.4}
\]

There is an even shorter proof of (3.3): digital Alexander gives

\[
P_2^4(1-p)=P_0^8(p),                                          \tag{3.5}
\]

while graph inclusion gives `P_0^8(p)<=P_0^4(p)`.

## 4. Matching-function reflection inequality

Recall

\[
M(p)=P_2(p)-P_0(p).                                           \tag{4.1}
\]

Adding (3.3) and (3.4) yields

\[
\boxed{M(p)+M(1-p)\le0.}                                     \tag{4.2}
\]

In particular

\[
\boxed{M(1/2)\le0.}                                          \tag{4.3}
\]

Since `M` is strictly increasing on every honest nontrivial torus, its unique zero obeys

\[
\boxed{p_\Lambda\ge1/2.}                                     \tag{4.4}
\]

Whenever the matching enhancement is strict at the rank level with positive probability—for example on the ordinary `L x L` square torus with `L>2`, where occupying the diagonal orbit gives a matching essential cycle but no NN edge—the inequalities are strict in the interior and

\[
p_\Lambda>1/2.                                                \tag{4.5}
\]

For arbitrary unusual period quotients, (4.4) is the unconditional statement; strictness only needs one configuration of positive product probability with `r_8>r_4`.

## 5. Birth-mixture reflection dominance

The fair birth mixture has CDF

\[
F(p)=\frac12[f_1(p)+f_2(p)]=\frac12[1+M(p)].                  \tag{5.1}
\]

Equation (4.2) is exactly

\[
\boxed{F(p)+F(1-p)\le1.}                                     \tag{5.2}
\]

For continuous birth labels, the reflected random variable `1-T` has CDF

\[
F_{1-T}(p)=1-F(1-p).                                          \tag{5.3}
\]

Thus

\[
\boxed{T\succeq_{st}1-T.}                                   \tag{5.4}
\]

The full finite mixture law, not merely its median, is reflection-shifted toward the high-`p` side.

## 6. Quantile and moment consequences

Let `Q` be the inverse CDF of the fair mixture.  From (5.2), for every `u in (0,1)`,

\[
\boxed{Q(u)+Q(1-u)\ge1.}                                     \tag{6.1}
\]

In particular

\[
Q(1/2)\ge1/2,                                                  \tag{6.2}
\]

\[
Q(1/4)+Q(3/4)\ge1.                                            \tag{6.3}
\]

Stochastic dominance also gives

\[
\boxed{E T\ge1/2,}                                            \tag{6.4}
\]

and, in the dual-odd coordinate

\[
C=\frac{T_1+T_2-1}{2},
\]

\[
\boxed{E C\ge0.}                                              \tag{6.5}
\]

Using the exact area identity

\[
E C=-\frac12\int_0^1M(p)\,dp,                                \tag{6.6}
\]

we obtain the finite integral sign

\[
\boxed{\int_0^1M(p)\,dp\le0.}                                \tag{6.7}
\]

Equivalently, in permutation birth indices,

\[
\boxed{E[K_1+K_2]\ge N+1.}                                   \tag{6.8}
\]

This last inequality can also be read directly from
`E K_1^8=N+1-E K_2^4` and `K_1^8<=K_1^4` in stochastic order.

## 7. Relationship to the strict mass-gap theorem

The finite reflection theorem gives the **sign** of the complement-odd centre for every size but not a uniform positive limiting gap.

In the fixed-positive exponential-aspect regime, `matching-enhancement-mass-gap-20260914.md` proves the stronger deterministic-centre statement

\[
a(d)+b(d)>1.                                                  \tag{7.1}
\]

The dilute directional calculation further gives an explicit large-`d` asymptotic for that positive displacement.  These are quantitative refinements of the exact finite reflection dominance, not independent sign coincidences.

## 8. A useful finite-data control

Any exact/Monte-Carlo rank-birth archive on an honest torus should satisfy, within its declared statistical errors,

\[
Q(u)+Q(1-u)\ge1                                               \tag{8.1}
\]

for every symmetric quantile pair, and

\[
E(K_1+K_2)\ge N+1.                                            \tag{8.2}
\]

Violations indicate a rank dictionary, matching-complement, weighting, or quantile-reconstruction error before they indicate new physics.
