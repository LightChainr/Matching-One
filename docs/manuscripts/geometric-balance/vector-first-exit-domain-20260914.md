# Vector first-exit certificates for the full exponential-moment / Wulff domain

2026-09-14.  Author-level proof of the finite-vector certificate requested in #766.  It uses only independent SITE product measure, site BK, translation invariance and the standard existence/uniform directional convergence of the subcritical inverse-correlation norm.

The main point is that the finite first-exit polynomial is naturally **multivariate** and certifies points of the full exponential-moment domain, not merely one axial mass.

## 1. Definitions

Work on a translation-invariant finite-range site graph `G` on `Z^2` at fixed subcritical occupation probability `p`.  Let

\[
\pi(x)=P_p(0\leftrightarrow x).                               \tag{1.1}
\]

For `t in R^2` define the exponential susceptibility

\[
\chi(t)=\sum_{x\in\mathbb Z^2}\pi(x)e^{t\cdot x},             \tag{1.2}
\]

and its convergence domain

\[
\mathcal D_p=\{t:\chi(t)<\infty\}.                            \tag{1.3}
\]

Let `S` be a finite set containing the origin.  For every **distinct external neighbour site** `v notin S` adjacent to at least one site of `S`, define

\[
b_S(v;p)
=P_p\left(0\text{ is connected inside }S
\text{ to an internal neighbour of }v\right).                 \tag{1.4}
\]

The site `v` itself is **not** required open in (1.4).  The origin occupation is included in the event probability.  Distinct external sites, not incident edges, are summed.

Set

\[
\boxed{
B_S(t;p)=\sum_{v\in\partial_{ext}S}b_S(v;p)e^{t\cdot v}.}      \tag{1.5}
\]

## 2. First-exit skeleton decomposition

**Theorem 2.1.**  If

\[
B_S(t;p)<1,                                                     \tag{2.1}
\]

then

\[
\boxed{
\chi(t)\le
\frac{\displaystyle\sum_{s\in S}e^{t\cdot s}}
     {1-B_S(t;p)}<\infty.}                                    \tag{2.2}
\]

In particular

\[
\boxed{\{t:B_S(t;p)<1\}\subset\mathcal D_p.}                 \tag{2.3}
\]

### Proof

Fix an open self-avoiding path `gamma` from `0` to `x`.  Starting at `y_0=0`, inspect the translated set `y_0+S`.

- If the target `x` belongs to `y_0+S`, stop.
- Otherwise let `y_1` be the first path vertex outside `y_0+S`.  The prefix immediately before `y_1` lies in `y_0+S` and witnesses the translated event (1.4) with exit displacement `v_1=y_1-y_0`.
- Restart from `y_1` with the translate `y_1+S`, and continue.

Because `gamma` is self-avoiding, this procedure terminates after a finite sequence of exit displacements

\[
v_1,\ldots,v_n\in\partial_{ext}S                              \tag{2.4}
\]

and a final residual displacement `s in S`, with

\[
x=v_1+\cdots+v_n+s.                                           \tag{2.5}
\]

Crucially, the witness sets of the successive first-exit events are disjoint.  The `i`th witness uses the path vertices from `y_{i-1}` up to the internal neighbour immediately before `y_i`; the exit site `y_i` itself is excluded from that event and becomes the starting occupied site of the next witness.  Hence consecutive witnesses do not share a SITE variable.

For any fixed skeleton `(v_1,...,v_n)`, the existence of such a path implies the disjoint occurrence of the translated increasing first-exit events.  Repeated site BK therefore bounds its probability by

\[
\prod_{i=1}^n b_S(v_i;p).                                     \tag{2.6}
\]

We impose no event on the final residual `s`; replacing it by probability one only enlarges the upper bound.

Union over all skeletons and final `s` gives the pointwise renewal bound

\[
\pi(x)
\le
\sum_{n\ge0}
\sum_{v_1,\ldots,v_n}
\sum_{s\in S}
1_{\{x=v_1+\cdots+v_n+s\}}
\prod_i b_S(v_i;p).                                           \tag{2.7}
\]

Multiply by `e^{t·x}` and sum over `x`.  All terms are nonnegative, so Tonelli applies and the right side factorizes:

\[
\chi(t)
\le
\left(\sum_{s\in S}e^{t\cdot s}\right)
\sum_{n\ge0}\left(\sum_vb_S(v;p)e^{t\cdot v}\right)^n.       \tag{2.8}
\]

Under (2.1) this is exactly (2.2).  `square`

No continuation path was assumed independent of the prefix.  The self-avoiding path plus BK is the reason the decomposition is legitimate even when later translated copies of `S` overlap geometrically.

## 3. Every finite certificate is convex

The function `B_S(t)` is a finite positive sum of exponentials and therefore convex.  Hence

\[
\mathcal C_S:=\{t:B_S(t)<1\}                                  \tag{3.1}
\]

is an open convex certified subset of `D_p`.

More importantly, the true domain itself is convex.  If `t_0,t_1 in D_p` and `0<lambda<1`, Holder gives

\[
\chi(\lambda t_0+(1-\lambda)t_1)
\le
\chi(t_0)^\lambda\chi(t_1)^{1-\lambda}.                       \tag{3.2}
\]

Therefore

\[
\boxed{\mathcal D_p\text{ is convex}.}                        \tag{3.3}
\]

Consequently certificates from different finite sets may be combined safely by

\[
\boxed{
\operatorname{conv}\left(\bigcup_i\mathcal C_{S_i}\right)
\subset\mathcal D_p.}                                        \tag{3.4}
\]

This is stronger than taking their union and answers the finite-domain combination question in #766.  One should take the convex hull **after** certifying each contributing point/set; convexifying the raw first-exit coefficients is a different operation and is not needed.

## 4. Relation to the directional inverse-correlation norm

Let `tau_p(x)` be the homogeneous subcritical connection norm,

\[
\tau_p(x)=\lim_{n\to\infty}-\frac1n\log\pi(\lfloor nx\rfloor).\tag{4.1}
\]

If `t in D_p`, each term in the convergent sum (1.2) is bounded, so along every direction

\[
t\cdot x\le\tau_p(x).                                        \tag{4.2}
\]

Thus

\[
\overline{\mathcal D_p}\subseteq
K_p:=\{t:t\cdot x\le\tau_p(x)\ \forall x\}.                  \tag{4.3}
\]

Conversely, under the standard uniform directional exponential estimate associated with the subcritical norm, every compact subset of the strict polar interior

\[
K_p^\circ=\{t:t\cdot e<\tau_p(e)\ \forall e\in S^1\}         \tag{4.4}
\]

has finite exponential susceptibility.  Hence

\[
\boxed{\mathcal D_p=K_p^\circ,\qquad
\overline{\mathcal D_p}=K_p.}                                \tag{4.5}
\]

In convex-analysis language, `K_p` is the dual unit ball associated with the norm `tau_p`, and

\[
\boxed{h_{K_p}(x)=\sup_{t\in K_p}t\cdot x=\tau_p(x).}          \tag{4.6}
\]

This is the correct support-function dictionary.

### Radial intercept is not the support function

For a Euclidean unit vector `e`, the ray intercept

\[
\rho_{K_p}(e)=\sup\{r\ge0:re\in K_p\}                         \tag{4.7}
\]

satisfies

\[
\rho_{K_p}(e)
=\inf_{x:e\cdot x>0}\frac{\tau_p(x)}{e\cdot x}.              \tag{4.8}
\]

In general

\[
\boxed{\rho_{K_p}(e)\ne\tau_p(e).}                            \tag{4.9}
\]

Equality requires the boundary normal at the radial point to line up with `e`, which holds on symmetry axes but not for a generic anisotropic direction.  A one-dimensional first-exit ray search must therefore not be labelled a direct support-function measurement unless this geometry is accounted for.

## 5. Large finite sets exhaust compact subsets of the true domain

Take any compact

\[
T\Subset\mathcal D_p=K_p^\circ.                               \tag{5.1}
\]

By compactness there is a margin `eta>0` such that

\[
\tau_p(e)-t\cdot e\ge3\eta
\qquad(t\in T,e\in S^1).                                     \tag{5.2}
\]

Let `S_R` be a large Euclidean or norm ball of radius `R`, enlarged by the finite interaction range.  For `v in partial_ext S_R`, the first-exit event is contained in a connection from the origin to a point at distance `R+O(1)`.  Uniform directional convergence to the subcritical norm gives, for all large `R`,

\[
b_{S_R}(v;p)
\le C e^{-\tau_p(v)+\eta|v|}.                                \tag{5.3}
\]

The number of boundary sites grows only polynomially in `R`.  Combining (5.2)--(5.3), uniformly for `t in T`,

\[
B_{S_R}(t;p)
\le \operatorname{poly}(R)e^{-\eta R}\longrightarrow0.       \tag{5.4}
\]

Therefore, for all sufficiently large `R`,

\[
\boxed{T\subset\mathcal C_{S_R}.}                             \tag{5.5}
\]

So the finite first-exit method is not only sound but **complete on compact interior subsets**: appropriately growing finite boxes eventually certify every point strictly inside the true exponential-moment/Wulff domain.

This theorem does not give a cheap state complexity for doing so; it answers the mathematical convergence question.

## 6. Numerical/certification workflow

A rigorous implementation can therefore proceed as follows.

1. For a finite `S`, compute outward-rounded upper bounds on each distinct-site coefficient `b_S(v;p)`.
2. The certified set is
   \[
   \sum_v \overline b_S(v;p)e^{t\cdot v}<1.
   \]
3. Combine several shapes/orientations by the convex hull of their certified regions.
4. For a direction `e`, obtain a **lower bound on `tau(e)`** from the support function of the certified convex body:
   \[
   \tau_p(e)\ge h_{C_{cert}}(e).                              \tag{6.1}
   \]
5. Keep this distinct from the radial intercept of the body.
6. Pair it with an independently proved finite-cylinder/connection upper bound on `tau(e)` when an interval is needed.

Because the support function of a polytope/convex certified body is a deterministic convex program, no angular finite difference is needed to obtain directional lower bounds.

## 7. Interfaces to the other 2026-09 notes

- `directional-enhancement-sandwich-20260914.md` gives the independent graph inclusion
  \[
  K_{8,p}\subseteq K_{4,p+\delta_I},
  \]
  which every numerical first-exit body should respect.
- `matrix-sewing-unit-residue-20260914.md` gives the conditional curvature/diffusion relation
  \[
  D^{-1}=\partial_{yy}\tau(1,0),
  \]
  once a common Markov-additive sewing kernel is established.
- `dilute-directional-geodesic-entropy-20260914.md` supplies an elementary small-`p` asymptotic for the support function in every fixed rational direction, giving another strong control on finite-box certificates in the dilute regime.

## 8. Claim boundary

The first-exit/BK theorem and convexity are self-contained product-measure arguments.  Equality of the domain closure with the polar body and the exhaustion statement use the standard uniform directional exponential-rate estimate for the subcritical norm; when citing them in the final manuscript, the precise SITE reference/hypotheses should be stated.  No numerical Wulff body is manufactured in this note.
