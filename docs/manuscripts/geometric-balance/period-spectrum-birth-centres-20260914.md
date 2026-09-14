# Birth centres from the minimum correlation-norm period spectrum

2026-09-14.  Abstract arbitrary-shape centre theorem obtained from the fixed-`p` homological-free-energy law.  No convergence of the minimizing period direction is required.

## 1. Minimum period cost as the geometry summary

For an honest period lattice `Lambda_n` with `N_n` sites, define for NN at `p<pc(G4)`

\[
\rho_{4,n}(p)
=\min_{\lambda\in\Lambda_n\setminus0}\tau_{4,p}(\lambda).    \tag{1.1}
\]

For matching at `q<pc(G8)`, define

\[
\rho_{8,n}(q)
=\min_{\lambda\in\Lambda_n\setminus0}\tau_{8,q}(\lambda).    \tag{1.2}
\]

At every finite `n`, `rho_{G,n}` is strictly decreasing in the occupation parameter: if `p<q` and `lambda_p` minimizes at `p`, then

\[
\rho_n(q)
\le\tau_q(\lambda_p)
<\tau_p(\lambda_p)
=\rho_n(p).                                                    \tag{1.3}
\]

The minimizer itself may change with `p` and need not be unique.

## 2. Lower birth: an abstract centre criterion

Assume

\[
\log N_n\to\infty.                                            \tag{2.1}
\]

Let `J4` be a compact subinterval of `(0,pc(G4))`.  Suppose

\[
\boxed{
\frac{\rho_{4,n}(p)}{\log N_n}
\longrightarrow g_4(p)}                                      \tag{2.2}
\]

locally uniformly for `p in J4`, where `g_4` is continuous and has a unique crossing

\[
\boxed{g_4(a)=1}                                              \tag{2.3}
\]

at an interior point `a in J4`, with

\[
g_4(p)>1\quad(p<a),
\qquad
g_4(p)<1\quad(p>a).                                        \tag{2.4}
\]

Assume also `rho_{4,n}(p)->infinity` on compact subsets of `J4`, which follows in the usual growing-shortest-period regimes.

Let `T_{1,n}` be the first NN ambient-rank birth.  Then

\[
\boxed{T_{1,n}\xrightarrow P a.}                             \tag{2.5}
\]

### Proof

For fixed `p<a`, (2.2)--(2.4) give

\[
\frac{\log N_n}{\rho_{4,n}(p)}<1-\delta                       \tag{2.6}
\]

for large `n`.  `general-period-homological-free-energy-20260914.md` gives

\[
P_p(r_4>0)\to0.                                               \tag{2.7}
\]

Thus

\[
P(T_{1,n}\le p)\to0.                                         \tag{2.8}
\]

For fixed `p>a`, the ratio is `>1+delta`, so

\[
P_p(r_4>0)\to1,                                               \tag{2.9}
\]

hence

\[
P(T_{1,n}\le p)\to1.                                         \tag{2.10}
\]

Two fixed points `a-epsilon,a+epsilon` trap the birth.

No direction label was used.

## 3. Upper birth from the complementary matching spectrum

Let

\[
q=1-p.                                                        \tag{3.1}
\]

Choose a compact interval `J8 subset (0,pc(G8))` and suppose

\[
\frac{\rho_{8,n}(q)}{\log N_n}
\to g_8(q)                                                    \tag{3.2}
\]

locally uniformly, with a unique crossing

\[
\boxed{g_8(c)=1}                                              \tag{3.3}
\]

at `c in J8`, decreasing from `>1` to `<1` as `q` increases.

For the second NN rank birth, digital Alexander gives

\[
\{T_{2,n}\le p\}
=\{r_4(p)=2\}
=\{r_8(1-p)=0\}.                                              \tag{3.4}
\]

The same free-energy argument for white matching therefore gives

\[
\boxed{T_{2,n}\xrightarrow P b=1-c.}                         \tag{3.5}
\]

So the two centre equations are simply

\[
\boxed{g_4(a)=1,\qquad g_8(1-b)=1.}                           \tag{3.6}
\]

## 4. Direction switching does not obstruct centre convergence

At a finite size,

\[
\rho_{G,n}(p)
=\min_{\lambda\ne0}\tau_{G,p}(\lambda)                       \tag{4.1}
\]

is a lower envelope of directional mass curves.  Its minimizing projective line can switch as `p` changes.  It can also vary with `n`.

The centre theorem needs none of the following:

- convergence of the minimizer direction;
- differentiability of the lower envelope;
- a unique cheapest period class;
- an OZ amplitude.

Only the scalar normalized minimum spectrum (2.2)/(3.2) and a unique crossing are needed.

This separates the problem into two layers:

1. **centre layer:** minimum period cost versus `log N`;
2. **window/mark layer:** which projective classes realize or nearly realize that minimum.

## 5. Unique versus competing slope at the centre

Let `a` be the lower centre.  Define at parameters `p_n->a` a minimizing line `ell_n` and its cheapest nonparallel cost

\[
\rho_{\perp,n}(p_n).                                          \tag{5.1}
\]

If

\[
\rho_{\perp,n}(p_n)-\log N_n\to+\infty,                      \tag{5.2}
\]

then `correlation-norm-successive-minima-20260914.md` gives

\[
P(r=1,L\ne\ell_n)\to0,                                      \tag{5.3}
\]

so the first-birth slope is asymptotically deterministic.

If instead several nonparallel projective classes satisfy

\[
\tau_{p_n}(\lambda)-\log N_n=O(1)                            \tag{5.4}
\]

simultaneously, centre convergence can still hold, but the birth mark requires a multi-direction crossover description.  The projective Poisson hard-core closure on this branch is designed for exactly that layer.

Thus direction degeneracy changes the **mark law**, not necessarily the centre.

## 6. Exponential shortest-period geometry as a special case

If a shortest period direction converges to `e` and

\[
\frac{\log N_n}{\ell_n}\to d,                                \tag{6.1}
\]

with all nonparallel periods much more expensive, then

\[
\rho_{4,n}(p)
=\ell_n\tau_{4,p}(e)+o(\ell_n).                              \tag{6.2}
\]

Hence

\[
g_4(p)=\tau_{4,p}(e)/d,                                      \tag{6.3}
\]

and `g_4(a)=1` is exactly

\[
\tau_{4,a}(e)=d.                                              \tag{6.4}
\]

The matching side is identical after complement.

So the directional centre theorem is one explicit realization of the abstract minimum-spectrum criterion.

## 7. A practical finite-size diagnostic

For a family of arbitrary period lattices, one can avoid guessing a direction by computing or certifying the scalar quantity

\[
R_{G,n}(p)=\frac{\rho_{G,n}(p)}{\log N_n}.                    \tag{7.1}
\]

A first-exit/Wulff inner certificate gives lower bounds on every `tau_p(lambda)` and hence on `rho_n`.  A finite connection construction gives upper bounds.

If the resulting intervals show a stable unique crossing of one, they certify the birth centre without resolving the full Wulff shape or choosing the minimizing direction in advance.

The slope/class should then be analysed separately through the certified near-minimizer set.

## 8. Boundary and regularity

At a point where the limiting minimum spectrum touches one without crossing, or where local uniform convergence fails, the theorem does not force a unique birth centre.

Likewise a centre theorem does not imply a `1/log N` or `1/ell` Gumbel window.  An affine window requires local intensity/derivative information about the near-minimizing classes.

This note is therefore a centre-level theorem, not a replacement for #763/#767.