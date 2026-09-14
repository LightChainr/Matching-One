# First-exit certificates as a convergent block-renormalization hierarchy

2026-09-14.  This note identifies the one-site first-exit certificate with the elementary all-walk exponential-tilt bound, and places the larger finite-box calculations of #761/#766 in one rigorous hierarchy converging to the true Wulff/exponential-moment domain.

## 1. The smallest block reproduces the random-walk Green-function bound

Take

\[
S=\{0\}.
\]

For every external neighbour `v`, the first-exit event in the definition of `b_S(v;p)` is simply “the origin is occupied”, because the origin is the only internal neighbour available and the external site `v` itself is excluded from the event.

Thus

\[
\boxed{b_{\{0\}}(v;p)=p.}                                    \tag{1.1}
\]

The one-site first-exit polynomial is therefore

\[
B_{\{0\}}(t;p)
=p\sum_{v\in S_G}e^{t\cdot v},                                \tag{1.2}
\]

where `S_G` is the graph step set.

The theorem in `vector-first-exit-domain-20260914.md` says

\[
B_{\{0\}}(t;p)<1
\quad\Longrightarrow\quad
\sum_xP(0\leftrightarrow x)e^{t\cdot x}<\infty.               \tag{1.3}
\]

But (1.2) is exactly the denominator condition in the all-walk Green-function/exponential-tilt bound used in the dilute directional analysis.  The two arguments are the same certificate written in path-skeleton and random-walk language.

## 2. Explicit NN and matching one-site domains

### NN

For steps `(+/-e_1,+/-e_2)`,

\[
\boxed{
B^{(4)}_{\{0\}}(t_x,t_y)
=2p(\cosh t_x+\cosh t_y).}                                   \tag{2.1}
\]

Hence the certified domain is

\[
\mathcal C_0^{(4)}(p)
=\{t:2p(\cosh t_x+\cosh t_y)<1\}.                             \tag{2.2}
\]

On the horizontal ray `(t,0)`, the boundary solves

\[
2p(\cosh t+1)=1,
\]

so

\[
\boxed{
t=\operatorname{arcosh}\left(\frac{1/p-2}{2}\right),}        \tag{2.3}
\]

which is exactly the elementary lower bound on `kappa_4(p)` in `dilute-directional-mass-centres-20260914.md`.

### Matching / king graph

The eight-neighbour step generating function factorizes:

\[
\sum_{v\in S_8}e^{t\cdot v}
=(1+2\cosh t_x)(1+2\cosh t_y)-1.                              \tag{2.4}
\]

Therefore

\[
\boxed{
B^{(8)}_{\{0\}}(t_x,t_y)
=p[(1+2\cosh t_x)(1+2\cosh t_y)-1].}                          \tag{2.5}
\]

The horizontal intercept solves

\[
p(6\cosh t+2)=1,
\]

or

\[
\boxed{
t=\operatorname{arcosh}\left(\frac{1/p-2}{6}\right),}        \tag{2.6}
\]

again exactly the dilute mass lower bound already derived.

Thus #766's vector certificate contains those earlier scalar estimates as its `S={0}` base case.

## 3. The dilute geodesic entropy is the support function of the one-site body

Let

\[
h_{\mathcal C_0}(u)=\sup_{t\in\mathcal C_0}t\cdot u.          \tag{3.1}
\]

The first-exit theorem gives

\[
\tau_p(u)\ge h_{\mathcal C_0(p)}(u).                           \tag{3.2}
\]

As `p->0`, the support function of the NN body (2.2) satisfies

\[
h_{\mathcal C_0^{(4)}}(a,b)
=(|a|+|b|)\log(1/p)
-(|a|+|b|)H_2\left(\frac{|a|}{|a|+|b|}\right)+o(1).           \tag{3.3}
\]

The matching body gives

\[
h_{\mathcal C_0^{(8)}}(a,b)
=M\log(1/p)-Mh_3(m/M)+o(1),                                   \tag{3.4}
\]

with `M=max(|a|,|b|)` and `m=min(|a|,|b|)`.

The adaptive-path lower-probability constructions in `dilute-directional-geodesic-entropy-20260914.md` prove matching upper bounds on `tau`.  Therefore the true Wulff support function is asymptotic to the **smallest-block first-exit certificate** in the dilute limit.

This explains why graph distance and shortest-path entropy are the first two dilute terms: before larger blocks matter, the exponential-moment domain is controlled by the bare step generating function.

## 4. Larger S resums local connectivity before the next exit

For a nontrivial finite `S`, the coefficient

\[
b_S(v;p)
\]

already sums every way the origin can connect inside `S` to an internal neighbour of `v`.  Thus moving from `S={0}` to a larger block replaces a single bare step by an **exact local connected passage**.

The first-exit skeleton then concatenates those passages using BK.  In renormalization language:

```text
one-site certificate    = bare walk kernel,
finite S certificate    = exact local connected block kernel,
large S                 = increasingly complete local resummation,
S -> infinity           = true exponential-moment/Wulff domain on compact interiors.
```

No equality between successive block kernels is asserted; different shapes can outperform each other in different directions.

## 5. A monotone accumulated certificate even though individual boxes need not nest

If `S_1 subset S_2`, it is **not** necessary or generally safe to assume

\[
\mathcal C_{S_1}\subseteq\mathcal C_{S_2}.                    \tag{5.1}
\]

The finite first-exit polynomials use different decompositions and may cross.

The true domain is convex, however, and every individual `C_S` is certified.  Therefore define the accumulated body

\[
\boxed{
K_R^{cert}
=\operatorname{conv}\left(
\bigcup_{S\in\mathfrak S_R}\mathcal C_S
\right),}                                                     \tag{5.2}
\]

where `mathfrak S_R` is any increasing family of finite block shapes/budgets.  Then

\[
K_R^{cert}\subseteq K_{R'}^{cert}\subseteq\mathcal D_p
\quad(R<R'),                                                   \tag{5.3}
\]

and if the family contains sufficiently large balls (or comparable exhausting shapes),

\[
\bigcup_RK_R^{cert}=\mathcal D_p                              \tag{5.4}
\]

at the level of compact interior exhaustion.

This gives a genuinely monotone numerical/theoretical hierarchy without demanding monotonicity of one chosen box sequence.

## 6. Shape selection should target the tilted deficit

For a target point `t` strictly inside the true domain, the boundary contribution of a large first-exit block at external displacement `v` has logarithmic form

\[
b_S(v)e^{t\cdot v}
\approx \exp[-\tau_p(v)+t\cdot v].                            \tag{6.1}
\]

Define the tilted deficit

\[
\Phi_t(v)=\tau_p(v)-t\cdot v.                                 \tag{6.2}
\]

The slowest boundary pieces are those with the smallest `Phi_t`.  Therefore a natural asymptotic block shape for certifying a fixed `t` is a discrete approximation to a level set

\[
\boxed{S_R(t)\approx\{x:\Phi_t(x)\le R\}.}                    \tag{6.3}
\]

Such a tilted-Wulff block equalizes the exponential difficulty around its boundary.  Rectangles, crosses or directionally elongated boxes are useful approximations when the true norm is not yet known, but they need not be optimal at a fixed state budget.

This is a design principle rather than a theorem that the dynamic-programming state count is minimized by (6.3); computational complexity depends on the chosen exact algorithm.  It does identify the correct **probability geometry** to target.

## 7. Adaptive certificate loop

A practical rigorous programme can iterate:

1. start with `S={0}` and obtain the bare convex body analytically;
2. add small exact boxes and convex-hull their certified domains;
3. use the current support lower bound as a proxy for `tau`;
4. construct a new integer block approximating a tilted-deficit level set in the direction where the certificate is weakest;
5. compute rigorous `b_S(v)` bounds and enlarge the convex hull;
6. stop when the desired directional support interval is narrower than the downstream centre/amplitude question requires.

Every stage remains a lower certificate even if the proxy shape is poor.  No fitted density is fed back as a rigorous mass input.

## 8. Interfaces

- #761 is the axial/small-number-of-parameters instance of this hierarchy.
- #766 is the full vector/convex implementation.
- `directional-enhancement-sandwich-20260914.md` supplies a graph-to-graph inclusion that every certified body should satisfy.
- `dilute-directional-geodesic-entropy-20260914.md` proves that the base `S={0}` body is already asymptotically sharp as `p->0`.

Together these turn finite first-exit boxes from an ad hoc numerical trick into a convergent block-renormalization scheme with a solved base case and a clear shape-selection principle.
