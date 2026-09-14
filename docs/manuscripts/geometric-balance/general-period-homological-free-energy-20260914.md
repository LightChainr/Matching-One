# General-period homological free energy at fixed subcritical parameter

2026-09-14.  This note closes the two-sided fixed-`p` version of the homological-free-energy picture for **arbitrary honest period lattices**.

The result does not assume a fixed direction, a rectangular basis, or an OZ prefactor.  It combines:

- the full-period first-exit theta upper bound;
- correlation-norm lattice packing;
- a finite angular net of fixed local directional connection seeds;
- Harris positive association and finite-group translation packing.

## 1. Statement

Fix a translation-invariant finite-range independent SITE graph `G` on `Z^2` and a parameter

\[
0<p<p_c(G).                                                    \tag{1.1}
\]

Let `Lambda_n` be honest rank-two integer period lattices with

\[
N_n=[\mathbb Z^2:\Lambda_n].                                 \tag{1.2}
\]

Let `tau_p` be the planar inverse-correlation norm and define

\[
\boxed{
\rho_n=ho_1(\Lambda_n;p)
=\min_{\lambda\in\Lambda_n\setminus0}\tau_p(\lambda).}       \tag{1.3}
\]

Assume

\[
\rho_n\to\infty                                               \tag{1.4}
\]

and

\[
\frac{\log N_n}{\rho_n}\to\alpha\in[0,\infty].              \tag{1.5}
\]

Let

\[
f_n(p)=P_p^{T_{\Lambda_n}}(r>0).                             \tag{1.6}
\]

**Theorem (fixed-p general-period free energy).**

\[
\boxed{
\lim_{n\to\infty}\frac1{\rho_n}\log f_n(p)
=-\max\{1-\alpha,0\}.}                                       \tag{1.7}
\]

Equivalently,

\[
\boxed{
\log f_n(p)
=-(\rho_n-\log N_n)_+ + o(\rho_n).}                           \tag{1.8}
\]

In particular:

- if `alpha<1`, positive homology is exponentially rare with exact exponent `rho_n-log N_n` at the `rho_n` scale;
- if `alpha>1`, positive homology occurs with probability tending to one;
- at `alpha=1`, `log f_n=o(rho_n)`; the boundary probability itself may depend on subexponential geometry.

## 2. Preliminary norm equivalence and a cheapest period

At fixed subcritical `p`, all norms on `R^2` are equivalent.  There are constants

\[
0<c_p\le C_p<\infty                                           \tag{2.1}
\]

such that

\[
c_p|x|\le\tau_p(x)\le C_p|x|.                              \tag{2.2}
\]

Choose

\[
\lambda_n\in\Lambda_n\setminus0,
\qquad
\tau_p(\lambda_n)=\rho_n.                                    \tag{2.3}
\]

Let

\[
L_n=|\lambda_n|.                                               \tag{2.4}
\]

Then

\[
\rho_n/C_p\le L_n\le\rho_n/c_p.                             \tag{2.5}
\]

So `L_n->infinity` and

\[
\log L_n=o(\rho_n).                                           \tag{2.6}
\]

Moreover the Euclidean shortest period also tends to infinity, so every fixed local seed used below eventually injects into the torus.

## 3. Upper bound

For any fixed `epsilon in (0,1)`, `first-exit-torus-winding-upper-20260914.md` gives

\[
f_n(p)
\le
N_n C_{p,\epsilon}
\sum_{\lambda\in\Lambda_n\setminus0}
 e^{-(1-\epsilon)\tau_p(\lambda)}.                            \tag{3.1}
\]

`correlation-norm-successive-minima-20260914.md` proves the norm-ball packing estimate

\[
\sum_{\lambda\ne0}
 e^{-(1-\epsilon)\tau_p(\lambda)}
\le
\widetilde C_{p,\epsilon}(\rho_n)
 e^{-(1-\epsilon)\rho_n},                                    \tag{3.2}
\]

with

\[
\log\widetilde C_{p,\epsilon}(\rho_n)=o(\rho_n).             \tag{3.3}
\]

Hence

\[
\limsup\frac1{\rho_n}\log f_n(p)
\le
\min\{0,\alpha-(1-\epsilon)\}.                               \tag{3.4}
\]

Let `epsilon downarrow0`:

\[
\boxed{
\limsup\frac1{\rho_n}\log f_n(p)
\le-\max\{1-\alpha,0\}.}                                    \tag{3.5}
\]

## 4. A finite family of fixed local directional seeds

The lower bound must produce a winding ring of class `lambda_n` at cost `rho_n+o(rho_n)` even though the minimizing direction may vary.

Fix `eta>0`.  Directional continuity of the norm on the compact unit circle permits a finite set of primitive integer directions

\[
w_1,\ldots,w_J                                              \tag{4.1}
\]

such that every `e in S^1` lies within a sufficiently small angular neighbourhood of some

\[
e_j=w_j/|w_j|,                                                \tag{4.2}
\]

and within that neighbourhood

\[
|\tau_p(e)-\tau_p(e_j)|<\eta.                                \tag{4.3}
\]

For each `j`, choose a fixed integer multiple

\[
z_j=M_jw_j                                                    \tag{4.4}
\]

so large that a finite-box connection seed from `0` to `z_j`, conditional on the initial site open, has probability `q_j` satisfying

\[
-\frac1{|z_j|}\log q_j
\le\tau_p(e_j)+\eta.                                         \tag{4.5}
\]

Then choose one fixed finite box supporting that seed.  The number, lengths and support sizes of all seed types are finite constants depending on `p,eta`, never on `n`.

Because the torus shortest Euclidean period tends to infinity, all these fixed seed boxes inject for large `n`.

## 5. Approximating the minimizing period by repeated seeds

Let

\[
e_n=\lambda_n/L_n.                                           \tag{5.1}
\]

Choose a seed direction `e_j` from the finite net close to `e_n`.  Let

\[
k_n=\operatorname{round}
\left(\frac{\lambda_n\cdot z_j}{|z_j|^2}\right)              \tag{5.2}
\]

and set

\[
r_n=\lambda_n-k_nz_j.                                        \tag{5.3}
\]

By making the angular net sufficiently fine as a function of `eta`,

\[
|r_n|\le c\eta L_n+O_{p,\eta}(1).                            \tag{5.4}
\]

Repeat the fixed connection seed `k_n` times along `z_j`, then append a deterministic NN path from `k_nz_j` to `lambda_n`.  The latter uses at most

\[
|r_n|_1\le\sqrt2|r_n|                                        \tag{5.5}
\]

sites up to an endpoint constant.  For a general finite-range graph `G`, use any fixed generating set path; its length is at most `C_G|r_n|` because the graph is connected and periodic.  For the actual NN/matching pair, the NN path is available in both.

The projected concatenation begins and ends at the same torus vertex and has lift displacement `lambda_n`, so it forces positive homology.

All seed and connector events are increasing.  Their supports may overlap after projection, but Harris positive association gives the product lower bound; overlaps can only reduce the number of distinct required open sites.

Using (4.3)--(5.5) and norm equivalence,

\[
\boxed{
P_p(\mathcal R_n)
\ge
\exp[-(1+C_p\eta)\rho_n-o(\rho_n)]}                          \tag{5.6}
\]

for one prescribed winding-ring attempt `R_n`.

## 6. The support costs only a polynomial number of translation centres

The union of the repeated fixed seed boxes and the deterministic connector lies in a Euclidean region of diameter `O_{p,eta}(L_n)`.

Therefore its projected support `S_n` obeys the crude but sufficient bound

\[
\boxed{|S_n-S_n|\le C_{p,\eta}L_n^2.}                         \tag{6.1}
\]

The finite-group translation-packing lemma supplies at least

\[
M_n
\ge
\frac{N_n}{C_{p,\eta}L_n^2}                                  \tag{6.2}
\]

pairwise site-disjoint translated attempts.  These attempts are independent.

By (2.6),

\[
\frac{\log M_n}{\rho_n}
=\frac{\log N_n}{\rho_n}+o(1)
\to\alpha.                                                    \tag{6.3}
\]

Thus the polynomial support footprint does not change the free-energy scale.

## 7. Lower bound for alpha<1

Let

\[
r_n=P_p(\mathcal R_n).                                       \tag{7.1}
\]

The probability that at least one disjoint attempt succeeds is

\[
1-(1-r_n)^{M_n}.                                               \tag{7.2}
\]

If `alpha<1`, choose `eta` so small that

\[
\alpha<1-C_p\eta.                                             \tag{7.3}
\]

Then `M_nr_n->0` exponentially.  Using `1-(1-r)^M >= Mr/2` when `Mr` is small,

\[
\liminf\frac1{\rho_n}\log f_n(p)
\ge\alpha-(1+C_p\eta).                                       \tag{7.4}
\]

Let `eta downarrow0`:

\[
\boxed{
\liminf\frac1{\rho_n}\log f_n(p)
\ge-(1-\alpha).}                                              \tag{7.5}
\]

Together with (3.5), this gives the exact rare-event exponent.

## 8. Lower bound for alpha>1

If `alpha>1`, choose `eta` so small that

\[
1+C_p\eta<\alpha.                                             \tag{8.1}
\]

Then

\[
M_nr_n\to\infty                                               \tag{8.2}
\]

exponentially, so

\[
\boxed{f_n(p)\to1.}                                          \tag{8.3}
\]

This proves the supercritical-in-opportunity side of the fixed-`p` homological transition.

## 9. Boundary alpha=1

For every fixed `eta>0`, the same construction gives

\[
\liminf\frac1{\rho_n}\log f_n(p)\ge-C_p\eta.                 \tag{9.1}
\]

Let `eta downarrow0`, while trivially `log f_n<=0`:

\[
\boxed{\frac1{\rho_n}\log f_n(p)\to0.}                       \tag{9.2}
\]

No universal limit of `f_n(p)` itself follows.  Polynomial/support factors and finer activity amplitudes live exactly at this boundary.

## 10. Interpretation

The free energy is

\[
\boxed{
\mathcal F_n(p)=\rho_1(\Lambda_n;p)-\log N_n.}               \tag{10.1}
\]

At fixed subcritical `p`, its sign determines the first homology event at exponential scale:

\[
\mathcal F_n\gg0\Rightarrow r=0\text{ whp},                  \tag{10.2}
\]

\[
\mathcal F_n\ll0\Rightarrow r>0\text{ whp}.                  \tag{10.3}
\]

More sharply,

\[
\log P(r>0)=-[\mathcal F_n]_+ +o(\rho_n).                    \tag{10.4}
\]

This is the rigorous fixed-`p` form of the earlier heuristic “connection energy minus log opportunities.”  The raw number of translations is `N`, while the need to reserve a support of diameter `O(L_n)` costs only `O(log L_n)` in the logarithm and therefore disappears at the `rho_n~L_n` scale.

## 11. Relation to previous geometry theorems

### Fixed-direction exponential torus

If `Lambda=<nu,mv>` and `log m/n->d`, then

\[
\rho_n=n\xi_u(p),
\qquad
\log N=\log(nm)=dn+o(n).                                     \tag{11.1}
\]

Hence

\[
\alpha=d/\xi_u(p),                                            \tag{11.2}
\]

and (1.7) becomes

\[
\frac1n\log P(r>0)
=-\max\{\xi_u(p)-d,0\},                                      \tag{11.3}
\]

recovering the fixed-direction theorem.

### Varying shortest direction

When the shortest period direction converges and the transverse height is exponential, the correlation-norm minimizer is asymptotically that shortest projective line and the same reduction gives the varying-direction rate.

### Geometric full-law criterion

At fixed `p`, norm equivalence makes `rho_1` comparable with the Euclidean shortest period `ell`.  Thus `log N/ell->0` implies `alpha->0`, but the present theorem is quantitatively sharper because it retains the actual anisotropic correlation norm and the exact constant one in the energy/entropy balance.

## 12. What this theorem does not yet give

`rho_1(Lambda;p)` depends on `p` and may change its minimizing projective direction as `p` varies.  Turning (1.7) into a universal birth-centre formula for an arbitrary changing lattice sequence requires control of the `p`-dependence of these minima and, near ties, the projective hard-core state space.

So (1.7) closes the fixed-parameter free energy for arbitrary shapes; it does not erase genuine direction-crossover questions in a moving-parameter window.

## 13. Claim boundary

The upper side uses the author-level first-exit/Wulff domain theorem.  The lower side uses only fixed local finite seeds, Harris, deterministic connectors and finite-group packing.  No OZ prefactor or analyticity is used.