# Persistent rank-one slope as a birth mark

2026-09-14.  Exact finite persistence consequences of the rank filtration on an honest torus.  No scaling limit, OZ input, continuum identification or numerical fit is used.

## 1. Monotone homology images freeze the rank-one line

Let `U_v` be iid continuous labels and let `X_G(p)` be the occupied subgraph at parameter `p` for graph `G` (NN or matching).  Write

\[
A_G(p)=\operatorname{im}\bigl(H_1(X_G(p))\to H_1(T^2)\bigr),
\qquad r_G(p)=\dim A_G(p).
\]

As `p` increases the occupied graph only grows, hence

\[
A_G(p)\subseteq A_G(q),\qquad p\le q.                         \tag{1.1}
\]

Let `T_1^G<=T_2^G` be the first parameters at which the rank reaches one and two.  On every realization with `T_1^G<T_2^G`, for every

\[
p,q\in(T_1^G,T_2^G)
\]

both `A_G(p)` and `A_G(q)` are one-dimensional and nested.  Therefore they are equal.  There is a single projective rational line

\[
\boxed{L_G\in\mathbb P(H_1(T^2;\mathbb Q))}                  \tag{1.2}
\]

such that

\[
A_G(p)=L_G\qquad\text{for all }p\in(T_1^G,T_2^G).             \tag{1.3}
\]

Thus the rank-one slope is born at the first homology birth and remains frozen until the second birth.  It cannot rotate, jump between primitive directions, or be re-sampled at intermediate `p` while the rank stays one.

If `T_1=T_2`, the rank-one plateau has zero Lebesgue width and no slope mark is needed.

## 2. Plateau width is the occupation time of the slope mark

Put

\[
G_G=T_2^G-T_1^G.
\]

For every bounded function `f` on projective primitive lines, the pathwise identity is

\[
\boxed{
\int_0^1 1_{\{r_G(p)=1\}}f(L_G)\,dp
=G_G f(L_G).}                                                  \tag{2.1}
\]

Taking expectations gives

\[
\boxed{
\int_0^1 E\bigl[1_{\{r_G(p)=1\}}f(L_G)\bigr]dp
=E[G_Gf(L_G)].}                                                \tag{2.2}
\]

For a particular primitive line `ell`,

\[
\boxed{
\int_0^1 P(r_G(p)=1,L_G=\ell)\,dp
=E[G_G1_{\{L_G=\ell\}}].}                                     \tag{2.3}
\]

For the projective spin-four readout of `projective-slope-harmonic-control-20260914.md`,

\[
\boxed{
\int_0^1 A_{4,G}(p;\tau)\,dp
=E\bigl[G_G Z_4^{(\tau)}(L_G)\bigr].}                          \tag{2.4}
\]

This gives an exact archive interface: integrating rank-one slope histograms over `p` measures the first-birth slope distribution biased by the plateau width, without any root fitting.

## 3. Uniform-p sampling produces a width-biased birth-mark law

Let `P` be uniform on `[0,1]`, independent of the labels, and condition on `r_G(P)=1`.  Then for any bounded `f`,

\[
E[f(L_G)\mid r_G(P)=1]
=\frac{E[G_G f(L_G)]}{E[G_G]}.                                \tag{3.1}
\]

Thus a slope histogram obtained by sampling a uniform occupation parameter inside the rank-one sector is **not** the unweighted law of the first-birth direction.  It is exactly the plateau-width-biased law.

This distinction matters whenever direction and plateau width are correlated.  The two distributions coincide only if `G_G` is independent of the birth slope (or is asymptotically deterministic independently of it).

## 4. Persistent 4/8 reflection with the slope mark

Let `V=1-U` be the reflected labels.  The persistent digital-Alexander theorem gives

\[
T_1^8(V)=1-T_2^4(U),
\qquad
T_2^8(V)=1-T_1^4(U),                                           \tag{4.1}
\]

and in the rank-one interval the complementary homology line is the same projective line.  Therefore, pathwise,

\[
\boxed{
G_8(V)=G_4(U),
\qquad
L_8(V)=L_4(U).}                                                \tag{4.2}
\]

With

\[
C_G=\frac{T_1^G+T_2^G-1}{2},
\]

we also have

\[
\boxed{
C_8(V)=-C_4(U).}                                               \tag{4.3}
\]

Since `V` has the same iid-uniform law as `U`, this becomes the exact distributional identity

\[
\boxed{
(G_8,L_8,C_8)\overset d=(G_4,L_4,-C_4).}                       \tag{4.4}
\]

In particular,

\[
\boxed{(G_8,L_8)\overset d=(G_4,L_4).}                        \tag{4.5}
\]

So **NN and matching have exactly the same finite-lattice joint law of rank-one plateau width and birth slope**.  Their difference is entirely in the complement-odd centre coordinate.

Every moment or transform of `(G,L)` is therefore shared between the two graphs.  For example

\[
E_4[G^m f(L)]=E_8[G^m f(L)]                                   \tag{4.6}
\]

for every nonnegative integer `m` and bounded `f`.

## 5. Direction competition should be formulated at the first birth

Because `L` is persistent, a multi-direction crossover does not consist of a rank-one system continuously rotating its winding direction as `p` changes.  The correct finite state is:

1. rank zero before `T_1`;
2. at `T_1`, one projective slope `L` is selected;
3. that same line persists throughout the whole rank-one plateau;
4. the second birth destroys the one-dimensional description by producing rank two.

Thus #765 direction competition and #767 common-window modelling should attach direction to the **first-birth mark**, not to each occupation parameter independently.

This also gives a clean falsification check for simulations that reconstruct the full rank process: any reported rank-one slope change before the second rank birth is a homology-tracking error.

## 6. A useful asymptotic conjecture for separated direction centres

Suppose an exponential-aspect sequence has finitely many candidate primitive directions `ell` with first-birth activities on comparable scales.  Let `lambda_ell(x)` denote their candidate intensities in a centred lower-birth window.  The natural conjecture is that the first-birth mark is chosen by a competing rare-event mechanism with weights proportional to these intensities, while the subsequent plateau width is then controlled by the gap to the second-rank mechanism.

The exact result above imposes two non-negotiable constraints on any such crossover model:

- the selected slope is constant through the entire rank-one interval;
- after complement/reflection, the joint `(plateau width, slope)` law must be unchanged.

A proposed direction-window law that violates either condition cannot be repaired by a different prefactor fit.

## 7. Claim boundary

Sections 1--5 are exact finite statements.  Section 6 is a research conjecture about how the exact persistent mark enters an exponential-aspect rare-event limit.  No independence between `G` and `L` is asserted.