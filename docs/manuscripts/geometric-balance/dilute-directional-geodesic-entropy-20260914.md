# Dilute directional norm = graph distance cost minus geodesic entropy

2026-09-14.  Elementary fixed-direction small-`p` theorem for the square NN and matching graphs.  This generalizes the axial `1` versus `3` dilute centre calculation in `dilute-directional-mass-centres-20260914.md`.

No OZ theorem or continuum limit is used.  The lower probability bound follows one adaptive forward path through independent layers; the upper probability bound sums all walks and uses a multivariate exponential tilt.

## 1. Entropy functions

For `alpha in [0,1]`, let

\[
H_2(\alpha)
=-\alpha\log\alpha-(1-\alpha)\log(1-\alpha)                   \tag{1.1}
\]

with the usual `0 log 0=0` convention.

For `beta in [-1,1]`, let

\[
h_3(\beta)
=\inf_{t\in\mathbb R}
\{\log(1+2\cosh t)-\beta t\}.                                \tag{1.2}
\]

This is the maximum entropy of a probability law on `{-1,0,+1}` with prescribed mean `beta`.  Equivalently, if `J_3` is the Cramer rate function of the uniform law on those three steps,

\[
h_3(\beta)=\log3-J_3(\beta).                                  \tag{1.3}
\]

It is even, strictly concave on `(-1,1)`, and

\[
h_3(0)=\log3,\qquad h_3(\pm1)=0.                              \tag{1.4}
\]

## 2. Theorem

Fix a nonzero integer vector `u=(a,b)`.  Put `A=|a|`, `B=|b|`.

For NN square-site percolation,

\[
\boxed{
\tau_{4,p}(u)
=(A+B)\log(1/p)
-(A+B)H_2\!\left(\frac{A}{A+B}\right)
+o(1),\qquad p\downarrow0.}                                   \tag{2.1}
\]

For matching/king adjacency, let

\[
M=\max(A,B),\qquad m=\min(A,B),\qquad \beta=m/M.
\]

Then

\[
\boxed{
\tau_{8,p}(u)
=M\log(1/p)-M h_3(\beta)+o(1),\qquad p\downarrow0.}            \tag{2.2}
\]

The theorem is homogeneous: these are costs for displacement `n u` divided by `n`.  Reflections and coordinate exchange reduce the proofs to `a>=b>=0` as needed.

The leading coefficient is the graph distance per copy of `u`; the constant correction is exactly the exponential entropy of shortest directed paths with that slope.

## 3. A multivariate all-walk upper bound on connection probability

For a finite step set `S`, any open connection from `0` to `n u` contains a self-avoiding path.  Replacing self-avoiding paths by all walks and using a vector exponential tilt `t in R^2` gives

\[
P_p(0\leftrightarrow n u)
\le p e^{-n t\cdot u}
\sum_{L\ge0}[p M_S(t)]^L,                                     \tag{3.1}
\]

where

\[
M_S(t)=\sum_{s\in S}e^{t\cdot s}.                             \tag{3.2}
\]

Thus, whenever `p M_S(t)<1`,

\[
\tau_p(u)\ge t\cdot u.                                       \tag{3.3}
\]

Taking the supremum gives

\[
\boxed{
\tau_p(u)\ge
\sup\{t\cdot u:pM_S(t)<1\}.}                                 \tag{3.4}
\]

We now evaluate this support function asymptotically as `p->0`.

## 4. NN lower mass bound

For NN,

\[
M_4(t_x,t_y)=2\cosh t_x+2\cosh t_y.                           \tag{4.1}
\]

Assume first `a,b>0`.  At the maximizing scale both tilts tend to `+infinity`, so writing `X=e^{t_x}`, `Y=e^{t_y}` reduces the leading constraint to

\[
X+Y=p^{-1}(1+o(1)).                                           \tag{4.2}
\]

Maximizing `a log X+b log Y` under `X+Y=p^{-1}` gives

\[
X=\frac{a}{a+b}p^{-1},
\qquad
Y=\frac{b}{a+b}p^{-1}.                                       \tag{4.3}
\]

The negative-exponential terms in the coshes contribute `o(1)` to the optimized objective.  Hence

\[
\sup_{pM_4<1}t\cdot u
=(a+b)\log(1/p)
+a\log\frac{a}{a+b}
+b\log\frac{b}{a+b}
+o(1),                                                        \tag{4.4}
\]

which is the right side of (2.1).  If one coordinate is zero, the same formula follows by taking the boundary optimizer; the entropy term is zero.

Therefore (3.4) gives the lower bound on `tau_4` required for (2.1).

## 5. NN adaptive forward path gives the matching upper mass bound

Assume `a,b>=0` and put `L=a+b`.  Use the layer coordinate `x+y`.  Starting from an occupied origin, at every step inspect the two forward neighbours

\[
(x+1,y),\qquad(x,y+1).
\]

If neither is open, fail; otherwise choose one open candidate uniformly using independent auxiliary randomness.

Each step inspects a fresh layer, so the two-site occupation pairs are independent across steps.  The one-step survival probability is

\[
q_2(p)=1-(1-p)^2=2p-p^2.                                     \tag{5.1}
\]

Conditional on survival, symmetry makes the chosen move exactly uniform between the two forward directions.  After `nL` successful steps, the endpoint is `n(a,b)` exactly when a Binomial `(nL,1/2)` count equals `na`.

Hence

\[
P_p(0\leftrightarrow n(a,b))
\ge p\,q_2(p)^{nL}
P\{\operatorname{Bin}(nL,1/2)=na\}.                           \tag{5.2}
\]

The binomial local large-deviation formula gives

\[
-\frac1n\log P\{\operatorname{Bin}(nL,1/2)=na\}
=L\,[\log2-H_2(a/L)]+o(1).                                   \tag{5.3}
\]

Since `-log q_2=-log(2p)+o(1)`, equations (5.2)--(5.3) yield

\[
\tau_{4,p}(u)
\le L\log(1/p)-L H_2(a/L)+o(1),                               \tag{5.4}
\]

matching (4.4) and proving (2.1).

The auxiliary randomization causes no issue: for every fixed occupation configuration, algorithmic success implies the existence of an open path, so averaging the success probability is a valid lower bound on the percolation event.

## 6. Matching adaptive path

Assume `a>=b>=0`; the other sectors follow by square symmetry.  Use columns as independent layers.  At each successful step increase `x` by one and inspect the three sites at vertical offsets `-1,0,+1` in the next column.

The survival probability is

\[
q_3(p)=1-(1-p)^3=3p-3p^2+p^3.                                \tag{6.1}
\]

Conditional on survival, the selected vertical increment is exactly uniform on `{-1,0,+1}`.  After `na` steps the endpoint is `(na,nb)` precisely when the increment sum equals `nb`.

If `J_3(beta)` is the Cramer rate for the uniform three-step law,

\[
P\{S_{na}=nb\}
=\exp[-na J_3(b/a)+o(n)].                                     \tag{6.2}
\]

Therefore

\[
\tau_{8,p}(a,b)
\le a[-\log q_3(p)+J_3(b/a)]                                 \tag{6.3}
\]

and, using `J_3=log3-h_3` and `-log q_3=-log(3p)+o(1)`,

\[
\tau_{8,p}(a,b)
\le a\log(1/p)-a h_3(b/a)+o(1).                               \tag{6.4}
\]

## 7. Matching all-walk support function

For king steps,

\[
M_8(t_x,t_y)
=(1+2\cosh t_x)(1+2\cosh t_y)-1.                             \tag{7.1}
\]

When `a>b>=0`, the optimal `t_x->+infinity` while `t_y` stays at the finite Cramer dual parameter associated with `beta=b/a`.  Uniformly for bounded `t_y`,

\[
M_8(t_x,t_y)
=e^{t_x}(1+2\cosh t_y)(1+o(1)).                               \tag{7.2}
\]

The constraint `pM_8=1` therefore gives

\[
t_x=\log(1/p)-\log(1+2\cosh t_y)+o(1).                       \tag{7.3}
\]

Substitution into `a t_x+b t_y` and optimization over `t_y` gives

\[
\sup_{pM_8<1}t\cdot(a,b)
=a\log(1/p)
-a\inf_t\{\log(1+2\cosh t)-(b/a)t\}
+o(1),                                                        \tag{7.4}
\]

which is

\[
a\log(1/p)-a h_3(b/a)+o(1).                                  \tag{7.5}
\]

At the boundary `b=a`, the dual parameter tends to `+infinity`; direct optimization of (7.1) gives the same leading value `a log(1/p)+o(1)`, consistent with `h_3(1)=0`.  Thus (3.4) matches the adaptive upper bound in every direction and proves (2.2).

## 8. Large-d centres for a fixed rational direction

Let

\[
r=\sqrt{a^2+b^2},
\qquad e=(a,b)/r.
\]

By homogeneity,

\[
\tau_p(e)=\tau_p(a,b)/r.
\]

### NN

Put `L=A+B` and `alpha=A/L`.  Solving

\[
\tau_{4,p}(e)=d
\]

with (2.1) gives

\[
\boxed{
p_{4,e}(d)
=\exp[-H_2(\alpha)]
\exp[-(r/L)d]\,[1+o(1)].}                                    \tag{8.1}
\]

### Matching

Put `M=max(A,B)`, `beta=min(A,B)/M`.  Solving the matching equation gives

\[
\boxed{
p_{8,e}(d)
=\exp[-h_3(\beta)]
\exp[-(r/M)d]\,[1+o(1)].}                                    \tag{8.2}
\]

These formulas are for `d->infinity` along a fixed rational direction.  They do not require a directional OZ amplitude.

## 9. Examples

### Axis

For `e=(1,0)`,

\[
H_2(1)=0,
\qquad h_3(0)=\log3.
\]

Hence

\[
p_{4,e}(d)\sim e^{-d},
\qquad
p_{8,e}(d)\sim\frac13e^{-d},                                 \tag{9.1}
\]

recovering `dilute-directional-mass-centres-20260914.md`.

### Diagonal

For `e=(1,1)/\sqrt2`,

\[
H_2(1/2)=\log2,
\qquad h_3(1)=0.
\]

Thus

\[
\boxed{
p_{4,e}(d)\sim\frac12 e^{-d/\sqrt2},
\qquad
p_{8,e}(d)\sim e^{-\sqrt2 d}.}                                \tag{9.2}
\]

The matching centre is exponentially smaller than the NN centre in `d`; the graph-enhancement separation is far stronger here than the axial factor-three difference.

### Generic genuinely tilted direction

If `A>B>0`, then

\[
\frac{r}{M}>\frac{r}{A+B}.
\]

Therefore

\[
\boxed{
\frac{p_{8,e}(d)}{p_{4,e}(d)}
=\exp[-c(e)d+O(1)]\to0}                                      \tag{9.3}
\]

for an explicit `c(e)>0`.  So at very large exponential aspect, matching enhancement creates exponentially separated occupation centres for every non-axial fixed rational direction.

## 10. Interface to the homological free energy

The leading dilute directional cost can be written

\[
\tau_{G,p}(u)
=d_G(u)\log(1/p)-s_G(u)+o(1),                                 \tag{10.1}
\]

where `d_G(u)` is graph distance per period and `s_G(u)` is geodesic entropy.  In the conjectural homological free energy

\[
\Psi_\Lambda(p)=\min_u\left\{\tau_{G,p}(u)-\log\frac{N}{|u|}\right\},
\]

the `p->0` race between primitive slopes is therefore controlled first by graph-distance geometry and then by geodesic entropy, before any closed-component sewing amplitude enters.

This gives a concrete asymptotic ordering of candidate homology classes in the large-`d` regime and a useful analytic control for #765-style tilted birth centres.
