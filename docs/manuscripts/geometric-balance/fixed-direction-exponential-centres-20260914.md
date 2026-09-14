# Fixed integer directions: exact exponential-aspect birth centres

2026-09-14.  Author-level continuation of #739/#765.  This closes the fixed-primitive-integer-direction case without rotating the microscopic interaction and without importing a bond OZ theorem.

The key new upper bound is a **cycle splitting + site-BK** argument.  It retains the full directional correlation-norm cost of the period vector and avoids the factor-two loss of a half-radius arm estimate.

## 1. Directional mass along a fixed integer vector

Fix a nonzero primitive integer vector

\[
u=(a,b)\in\mathbb Z^2,
\]

and choose `v in Z^2` with

\[
\det(u,v)=1.                                                   \tag{1.1}
\]

For `G=G4` (NN) or `G8` (matching), let

\[
\tau^G_p(x)=P_p^G(0\leftrightarrow x),
\qquad
s^G_p(x)=\tau^G_p(x)/p.                                       \tag{1.2}
\]

Conditioning a common endpoint open and using Harris gives

\[
s_p(x+y)\ge s_p(x)s_p(y).                                    \tag{1.3}
\]

Therefore

\[
\boxed{
\xi_{G,u}(p)
=\lim_{k\to\infty}-\frac1k\log s_p^G(ku)
=\inf_{k\ge1}-\frac1k\log s_p^G(ku)}                          \tag{1.4}
\]

exists.  It is the correlation-norm cost per integer multiple of `u`.  If

\[
e=u/|u|,
\]

write

\[
\tau_{G,p}(e)=\xi_{G,u}(p)/|u|.                               \tag{1.5}
\]

Subcritical sharpness makes `xi>0` for `p<pc(G)`.  The stable extension

\[
\xi_{G,p}(x)
=\lim_{k\to\infty}-k^{-1}\log s_p^G(kx)                       \tag{1.6}
\]

is a norm on `R^2` in the subcritical regime, and

\[
s_p(x)\le e^{-\xi_{G,p}(x)}                                  \tag{1.7}
\]

for integer `x`, because (1.4) is the infimum of the finite-step costs.

## 2. The exponential torus adapted to u

For integers `n,m` define

\[
\Lambda_{n,m}=\langle nu,mv\rangle.                           \tag{2.1}
\]

Its index is

\[
N_{n,m}=nm.                                                     \tag{2.2}
\]

Assume

\[
n\to\infty,
\qquad
m=m_n\to\infty,
\qquad
\frac{\log m}{n}\to d\in(0,\infty).                          \tag{2.3}
\]

Because the component of `mv` perpendicular to `u` has magnitude `m/|u|`, every period with nonzero `v` coefficient has Euclidean length at least `m/|u|`.  Hence for exponentially large `m`, the only periods on the `O(n)` scale are the multiples of `nu`.

The physical shortest length is

\[
\ell_n=n|u|,                                                    \tag{2.4}
\]

and the transverse height is

\[
h_n=N/\ell_n=m/|u|.                                           \tag{2.5}
\]

Thus if #765 uses

\[
\frac{\log h_n}{\ell_n}\to\delta,                            \tag{2.6}
\]

then

\[
d=\delta|u|.                                                   \tag{2.7}
\]

## 3. A full-period upper bound by splitting a simple homology cycle

Let

\[
f_G(n,m;p)=P_p^G(r_G>0)                                      \tag{3.1}
\]

on `Z^2/Lambda_{n,m}`.

If `r_G>0`, choose a nonzero-homology occupied closed walk with the fewest edges.  It cannot repeat a torus vertex except at its start/end: otherwise a repeated vertex splits it into two shorter closed walks, at least one of which has nonzero homology.

Hence the chosen walk is a graph-simple cycle.  Lift it from a starting vertex `z` to a path ending at

\[
z+\lambda,
\qquad
\lambda\in\Lambda_{n,m}\setminus\{0\}.                       \tag{3.2}
\]

Split the cycle at an internal vertex.  If the first lifted arc has displacement `x`, the second has displacement `lambda-x`.  Their internal torus vertices are disjoint.

Condition the two common endpoints open.  On all other sites, the two arc-connection events occur disjointly.  The site-BK inequality therefore gives, up to the exact endpoint normalization,

\[
P(\text{the two arcs occur with displacements }x,\lambda-x)
\le s_p(x)s_p(\lambda-x).                                    \tag{3.3}
\]

Using (1.7) and the triangle inequality for the correlation norm,

\[
\boxed{
s_p(x)s_p(\lambda-x)
\le e^{-\xi_p(x)-\xi_p(\lambda-x)}
\le e^{-\xi_p(\lambda)}.}                                    \tag{3.4}
\]

This is the point of the split: the **full homology vector lambda** remains in the exponent.

### 3.1 Summing over the split does not change the exponential rate

Define

\[
C_p(\lambda)=\sum_{x\in\mathbb Z^2}
 e^{-\xi_p(x)-\xi_p(\lambda-x)}.                              \tag{3.5}
\]

Let `L=xi_p(lambda)`.  For `k>=0`, the set

\[
\{x:L+k\le\xi_p(x)+\xi_p(\lambda-x)<L+k+1\}
\]

lies inside a correlation-norm ball of radius `L+k+1`.  A two-dimensional norm ball contains `O((L+k+2)^2)` lattice points.  Therefore

\[
\boxed{
C_p(\lambda)\le C_p'(1+\xi_p(\lambda))^2e^{-\xi_p(\lambda)}.} \tag{3.6}
\]

The polynomial factor is harmless on the `n` exponential scale.

### 3.2 The period-lattice theta sum is dominated by +-nu

For `lambda=alpha nu+beta mv`, if `beta=0`, homogeneity gives

\[
\xi_p(\lambda)=|\alpha|n\xi_{G,u}(p).                         \tag{3.7}
\]

If `beta!=0`, norm equivalence and the perpendicular component give

\[
\xi_p(\lambda)\ge c_p\frac{|\beta|m}{|u|}.                   \tag{3.8}
\]

Summing over the parallel coordinate `alpha` only contributes a bounded one-dimensional lattice factor.  Consequently

\[
\sum_{\lambda\in\Lambda_{n,m}\setminus\{0\}}
(1+\xi_p(\lambda))^2e^{-\xi_p(\lambda)}
\le C_p n^2e^{-n\xi_{G,u}(p)}+e^{-c_pm/C_p}.                  \tag{3.9}
\]

There are `N=nm` possible starting torus vertices.  Combining (3.3)--(3.9),

\[
\boxed{
f_G(n,m;p)
\le C_p\,nm\,n^2e^{-n\xi_{G,u}(p)}
   +nm\,e^{-c_pm/C_p}.}                                       \tag{3.10}
\]

Thus the opportunity entropy is exactly `log m`; no second factor of `m` is lost to endpoint enumeration.

## 4. Lower bound: repeat one finite directional seed and close the seam

Fix `p<pc(G)` and `epsilon>0`.  From (1.4), choose a fixed integer `L` such that

\[
-\frac1L\log s_p(Lu)
<\xi_{G,u}(p)+\epsilon/4.                                    \tag{4.1}
\]

Exhaust the full-plane connection by finite boxes around the segment from `0` to `Lu`.  Choose one fixed finite box `B` for which

\[
q=P_p(0\leftrightarrow Lu\text{ in }B\mid0\text{ open})
>e^{-(\xi_{G,u}(p)+\epsilon/2)L}.                             \tag{4.2}
\]

Write

\[
n=kL+r,\qquad0\le r<L.                                       \tag{4.3}
\]

Translate this seed box `k` times by `Lu`.  Requiring the successive seed connections and a fixed finite connector for the remaining displacement `ru` produces a connected lifted path from `0` to `nu`.  Its projection is a nonzero winding ring on the torus.

The translated seed events overlap only through bounded regions and common endpoints.  Conditional Harris, exactly as in the axial proof, gives for all large `n`

\[
P_p(\text{one prescribed u-ring})
\ge e^{-(\xi_{G,u}(p)+\epsilon)n}.                            \tag{4.4}
\]

The complete ring support lies in a fixed-width tube around `[0,nu]`.  Therefore its projected difference set has cardinality

\[
|S-S|\le C_{p,\epsilon,u}n.                                   \tag{4.5}
\]

The finite-group packing lemma supplies at least

\[
\frac{N}{|S-S|}\ge c_{p,\epsilon,u}m                         \tag{4.6}
\]

pairwise site-disjoint translates.  Their ring events are independent.  Hence

\[
\boxed{
1-f_G(n,m;p)
\le
\exp\left[-c_{p,\epsilon,u}m
 e^{-(\xi_{G,u}(p)+\epsilon)n}\right].}                       \tag{4.7}
\]

## 5. Exact directional winding-event rate

Combining (3.10) and (4.7), for every fixed subcritical `p`,

\[
\boxed{
\lim_{n\to\infty}\frac1n\log f_G(n,m_n;p)
=-\max\{\xi_{G,u}(p)-d,0\}.}                                 \tag{5.1}
\]

At equality `xi=d`, only the zero exponential rate is asserted; no limit probability is claimed.

In the Euclidean normalization of #765, using (2.7), this becomes

\[
\boxed{
\lim\frac1{\ell_n}\log f_G
=-\max\{\tau_{G,p}(e)-\delta,0\}.}                            \tag{5.2}
\]

So the fixed-direction centre is governed by the actual directional correlation norm, not by Euclidean angle times the axial mass.

## 6. Strict p-monotonicity in every fixed integer direction

Continuity of `p -> xi_{G,u}(p)` on `(0,pc(G))` follows by the same site arguments used in `exponential-birth-centres.md`:

- monotonicity plus the finite-seed infimum (1.4) gives left continuity;
- the uniform subcritical cluster-volume tail gives right continuity after truncating the connecting cluster at size `O(n)`.

For strict decrease, suppose

\[
0<p<q<pc(G),
\qquad
\xi_{G,u}(p)=\xi_{G,u}(q)=A>0.                               \tag{6.1}
\]

Choose `d<A` sufficiently close to `A`, and then `zeta>0`, so the Friedgut--Kalai transitive sharp-threshold increment

\[
\rho\frac{A-d+\zeta}{d}<q-p.                                 \tag{6.2}
\]

Take `m=ceil(e^{dn})`.  By (5.1), at parameter `p`,

\[
f_G(n,m;p)>\epsilon_n=e^{-(A-d+\zeta)n}                       \tag{6.3}
\]

for all large `n`.  The event `r_G>0` is increasing and invariant under the transitive torus translation group on `N=nm` sites.  Friedgut--Kalai therefore forces

\[
f_G(n,m;q)>1-\epsilon_n\to1.                                 \tag{6.4}
\]

But (3.10) with `xi(q)=A>d` gives `f_G(n,m;q)->0`, a contradiction.

Thus

\[
\boxed{p\mapsto\xi_{G,u}(p)\text{ is strictly decreasing on }(0,pc(G)).} \tag{6.5}
\]

The small-`p` path count makes `xi->infinity` as `p->0`.  As `p->pc(G)`, square symmetry and the norm triangle inequality give

\[
\xi_{G,u}(p)
\le (|a|+|b|)\kappa_G(p)\to0,                                \tag{6.6}
\]

where the axial mass tends to zero by the already proved site argument.  Hence `xi_{G,u}` continuously and strictly maps `(0,pc(G))` onto `(infinity,0)`.

## 7. Birth centres for the fixed direction

Let `T_1,T_2` be the two NN rank births on `Lambda_{n,m}`.  Since

\[
P(T_1\le p)=f_{G4}(n,m;p),                                    \tag{7.1}
\]

(5.1) and strict invertibility imply

\[
\boxed{T_1\xrightarrow{P}a_u(d),
\qquad
\xi_{4,u}(a_u(d))=d.}                                        \tag{7.2}
\]

Digital Alexander duality gives

\[
P_p^{G4}(r=2)=P_{1-p}^{G8}(r=0),                              \tag{7.3}
\]

so the second birth satisfies

\[
\boxed{T_2\xrightarrow{P}b_u(d),
\qquad
b_u(d)=1-c_u(d),
\qquad
\xi_{8,u}(c_u(d))=d.}                                        \tag{7.4}
\]

Equivalently, with Euclidean direction `e=u/|u|` and `delta=d/|u|`,

\[
\boxed{
\tau_{4,a(e,\delta)}(e)=\delta,
\qquad
\tau_{8,1-b(e,\delta)}(e)=\delta.}                           \tag{7.5}
\]

This is the fixed-integer-direction centre theorem requested in #765.

## 8. Strict matching mass gap in every fixed direction

`directional-enhancement-sandwich-20260914.md` gives, on compact subintervals of `(0,pc(G8))`, a positive `delta_p` independent of the endpoint direction such that

\[
\tau_{8,p}(x)\ge\tau_{4,p+\delta_p}(x)                        \tag{8.1?}
\]

at the probability level in the direction that yields, after taking logarithmic rates,

\[
\xi_{8,u}(p)\le\xi_{4,u}(p+\delta_p).                         \tag{8.1}
\]

By the strict directional monotonicity just proved,

\[
\xi_{4,u}(p+\delta_p)<\xi_{4,u}(p).                           \tag{8.2}
\]

Therefore

\[
\boxed{
\xi_{8,u}(p)<\xi_{4,u}(p),
\qquad0<p<pc(G8),}                                            \tag{8.3}
\]

for every fixed primitive integer direction `u`.

Consequently the two fixed-direction exponential-aspect centres satisfy the strict complement-odd displacement analogue of the axial theorem.

## 9. What remains of #765

This note closes the case where the short period direction is a fixed primitive integer line and the transverse period is an exponentially large multiple of one fixed Bezout complement.

The genuinely remaining problem is **varying directions** `u_n/|u_n| -> e` with arithmetic changing with `n`.  There one must control the correlation norm and finite-seed constants uniformly along the direction sequence.  The issue's warning about several comparable period classes remains relevant when the period lattice itself is not of the separated form (2.1).

No angle scan is needed for the fixed-direction theorem.