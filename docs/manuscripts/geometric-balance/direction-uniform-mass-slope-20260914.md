# A direction-uniform quantitative occupation slope for the correlation norm

2026-09-14.  Quantitative strengthening of the all-direction strict `p`-monotonicity proved in the directional centre theorem.

The same Friedgut--Kalai constant that yielded the axial mass-slope inequality gives a multiplicative comparison **uniform in direction**.

## 1. Statement

Let `G` be NN or matching SITE percolation on `Z^2`.  Let

\[
\tau_p(e),\qquad e\in S^1,\quad0<p<p_c(G)                    \tag{1.1}
\]

be the subcritical inverse-correlation norm on Euclidean unit directions.

Let `rho_FK>0` be a universal constant for the transitive Friedgut--Kalai sharp-threshold inequality in the convention used in `exponential-birth-centres.md`.

**Theorem.**  For every

\[
0<p<q<p_c(G)                                                   \tag{1.2}
\]

and every `e in S^1`,

\[
\boxed{
\tau_p(e)-\tau_q(e)
\ge
\frac{q-p}{\rho_{FK}}\,\tau_q(e).}                           \tag{1.3}
\]

Equivalently,

\[
\boxed{
\tau_p(e)
\ge
\left(1+\frac{q-p}{\rho_{FK}}\right)\tau_q(e).}              \tag{1.4}
\]

By homogeneity, the same inequality holds for every `x in R^2`.

## 2. Proof from the varying-direction torus rate theorem

Fix `e` and abbreviate

\[
A=\tau_p(e),
\qquad
B=\tau_q(e).                                                  \tag{2.1}
\]

Monotonicity gives `A>=B>0`.  Suppose for contradiction that

\[
A-B<\frac{q-p}{\rho_{FK}}B.                                  \tag{2.2}
\]

Then

\[
\rho_{FK}\frac{A-B}{B}<q-p.                                  \tag{2.3}
\]

Choose a number

\[
d<B                                                        \tag{2.4}
\]

sufficiently close to `B`, and then `zeta>0` sufficiently small, so that

\[
\boxed{
\rho_{FK}\frac{A-d+\zeta}{d}<q-p.}                           \tag{2.5}
\]

Take the sequence of honest integer-period tori used in `varying-direction-exponential-centres-20260914.md`, with shortest directions tending to `e`, short length `ell_n`, and

\[
\frac{\log h_n}{\ell_n}\to d.                                \tag{2.6}
\]

At parameter `p`, the directional rate theorem gives

\[
P_p(r>0)
=\exp[-(A-d+o(1))\ell_n].                                    \tag{2.7}
\]

Therefore for large `n`,

\[
P_p(r>0)>
\epsilon_n:=e^{-(A-d+\zeta)\ell_n}.                          \tag{2.8}
\]

The event `r>0` is increasing and invariant under the transitive translation action on the `N_n` torus sites.  Also

\[
\frac{\log N_n}{\ell_n}\to d,                                \tag{2.9}
\]

because the extra `log ell_n/ell_n` vanishes.  The Friedgut--Kalai threshold increment for probability `epsilon_n` therefore tends to

\[
\rho_{FK}\frac{A-d+\zeta}{d},                                \tag{2.10}
\]

which is strictly less than `q-p` by (2.5).

Hence at parameter `q`,

\[
P_q(r>0)>1-\epsilon_n\to1.                                   \tag{2.11}
\]

But `d<B=tau_q(e)`, so the same directional rate theorem gives

\[
P_q(r>0)\to0.                                                 \tag{2.12}
\]

Contradiction.  Thus (1.3) holds.

No differentiability of `tau` is used.

## 3. Direction-uniform logarithmic slope

Rearrange (1.4):

\[
\frac{\tau_p(e)}{\tau_q(e)}
\ge1+\frac{q-p}{\rho_{FK}}.                                  \tag{3.1}
\]

At a differentiability point in `p`, let `q downarrow p`.  Then

\[
\boxed{
-\partial_p\tau_p(e)
\ge\frac{\tau_p(e)}{\rho_{FK}}>0.}                            \tag{3.2}
\]

Equivalently,

\[
\boxed{
-\partial_p\log\tau_p(e)\ge\frac1{\rho_{FK}}.}               \tag{3.3}
\]

This is uniform in direction.

The inequality supplies a positive lower bound on every regular directional Gumbel scale slope once the centre mass is fixed.

## 4. Minimum period spectrum inherits the same inequality

For any fixed period lattice `Lambda`, define

\[
\rho_1(p)
=\min_{\lambda\in\Lambda\setminus0}\tau_p(\lambda).          \tag{4.1}
\]

Since (1.4) holds for every period vector,

\[
\tau_p(\lambda)
\ge
\left(1+\frac{q-p}{\rho_{FK}}\right)\tau_q(\lambda).          \tag{4.2}
\]

Taking minima gives

\[
\boxed{
\rho_1(p)
\ge
\left(1+\frac{q-p}{\rho_{FK}}\right)\rho_1(q).}              \tag{4.3}
\]

Thus the scalar period-spectrum cost has a quantitative strict decrease even when its minimizing projective direction switches.

## 5. Consequence for normalized period-spectrum limits

Let `N_n` be any period-lattice sequence and suppose on an interval `J`

\[
\frac{\rho_{1,n}(p)}{\log N_n}\to g(p)                       \tag{5.1}
\]

pointwise (local uniform convergence is useful for centre trapping but not needed for this inequality).

Passing to the limit in (4.3),

\[
\boxed{
g(p)
\ge
\left(1+\frac{q-p}{\rho_{FK}}\right)g(q).}                   \tag{5.2}
\]

Whenever `g(q)>0`,

\[
\boxed{g(p)>g(q)\quad(p<q).}                                 \tag{5.3}
\]

So a positive limiting minimum-period spectrum cannot contain a flat segment.

In particular, if `g` takes values above and below one, the crossing

\[
g(p)=1                                                       \tag{5.4}
\]

is automatically unique.

This removes the “unique crossing” hypothesis from `period-spectrum-birth-centres-20260914.md` once positivity and existence of the limiting spectrum are known.

## 6. Quantitative movement of inverse centres

Suppose `a(d)` is defined by

\[
\tau_{a(d)}(e)=d.                                             \tag{6.1}
\]

For `d_1>d_2>0`, put

\[
p_i=a(d_i),                                                   \tag{6.2}
\]

so `p_1<p_2`.  Applying (1.4),

\[
d_1
\ge
\left(1+\frac{p_2-p_1}{\rho_{FK}}\right)d_2.                 \tag{6.3}
\]

Therefore

\[
\boxed{
p_2-p_1
\le\rho_{FK}\left(\frac{d_1}{d_2}-1\right).}                 \tag{6.4}
\]

This is a rough universal modulus for the inverse directional centre map, independent of angle and microscopic adjacency within the two site models.

## 7. Interaction with the matching enhancement gap

On a compact parameter interval inside `(0,pc(G8))`, the matching enhancement gives a density sprinkling `delta_I>0` with

\[
\tau_{8,p}(e)\le\tau_{4,p+\delta_I}(e).                       \tag{7.1}
\]

Apply (1.3) to `G4`:

\[
\tau_{4,p}(e)-\tau_{4,p+\delta_I}(e)
\ge
\frac{\delta_I}{\rho_{FK}}\tau_{4,p+\delta_I}(e).             \tag{7.2}
\]

Hence the same-parameter mass gap has the quantitative lower bound

\[
\boxed{
\tau_{4,p}(e)-\tau_{8,p}(e)
\ge
\frac{\delta_I}{\rho_{FK}}\tau_{4,p+\delta_I}(e)>0.}          \tag{7.3}
\]

uniformly in direction on the compact parameter interval.

This supplies an explicit conceptual source for the Wulff-body Minkowski buffer in `strict-wulff-body-separation-20260914.md`.

## 8. Claim boundary

The only quantitative constant is the universal Friedgut--Kalai constant already imported in the parent manuscript; it is not numerically optimized here.  The theorem uses the author-level varying-direction winding-rate theorem on PR #771, but no OZ or p-analyticity input.