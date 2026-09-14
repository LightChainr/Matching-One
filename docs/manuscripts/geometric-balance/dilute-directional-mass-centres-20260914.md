# Elementary dilute directional masses and the large-d birth centres

2026-09-14.  A rigorous small-`p` comparison using only one adaptive forward path and an all-walk exponential-tilt bound.  No Ornstein--Zernike theorem, component sewing, or prefactor assumption is used.

## 1. Result

Let `kappa_4(p)` and `kappa_8(p)` be the horizontal inverse correlation lengths of square-site NN and matching (king) adjacency.  As `p downarrow 0`,

\[
\boxed{\kappa_4(p)=-\log p+O(p),}                              \tag{1.1}
\]

\[
\boxed{\kappa_8(p)=-\log(3p)+O(p).}                           \tag{1.2}
\]

Consequently, for the exponential-aspect centres

\[
a(d)=\kappa_4^{-1}(d),
\qquad
c(d)=\kappa_8^{-1}(d)=1-b(d),
\]

one has as `d->infinity`,

\[
\boxed{a(d)=e^{-d}+O(e^{-2d}),}                               \tag{1.3}
\]

\[
\boxed{c(d)=\frac13e^{-d}+O(e^{-2d}),}                        \tag{1.4}
\]

and therefore

\[
\boxed{a(d)+b(d)-1
=a(d)-c(d)
=\frac23e^{-d}+O(e^{-2d}).}                                   \tag{1.5}
\]

The separated birth gap also has

\[
\boxed{b(d)-a(d)=1-\frac43e^{-d}+O(e^{-2d}).}                 \tag{1.6}
\]

Thus the strict centre asymmetry from the enhancement theorem has a concrete elementary extreme-elongation scale.

## 2. A general all-walk exponential-tilt bound

Let a finite-range graph have step set `S subset Z^2` and independent site density `p`.  If `0` is connected to `ne_1`, there is a self-avoiding open path from `0` to `ne_1`.  A path of `L` edges uses `L+1` open vertices, so a union bound followed by replacing self-avoiding paths by all walks gives, for every `t>0`,

\[
P_p(0\leftrightarrow ne_1)
\le p e^{-tn}\sum_{L\ge0}
\left[p\sum_{s\in S}e^{t s_x}\right]^L.                       \tag{2.1}
\]

Whenever

\[
pM_S(t)<1,
\qquad
M_S(t)=\sum_{s\in S}e^{t s_x},                                \tag{2.2}
\]

this yields

\[
P_p(0\leftrightarrow ne_1)
\le\frac{p}{1-pM_S(t)}e^{-tn}.                                \tag{2.3}
\]

Hence

\[
\boxed{\kappa(p)\ge t_*},                                    \tag{2.4}
\]

where `t_*` is the positive solution of `pM_S(t_*)=1`.  Taking `t<t_*` and then `t upward t_*` avoids using the divergent geometric series at equality.

This is only a path-count upper bound on connectivity probability; it is not an OZ statement.

## 3. NN square graph

For NN steps

\[
S_4=\{(\pm1,0),(0,\pm1)\},
\]

so

\[
M_4(t)=2\cosh t+2.                                             \tag{3.1}
\]

Equation `pM_4(t_4)=1` gives

\[
\boxed{
t_4(p)=\operatorname{arcosh}\left(\frac{1/p-2}{2}\right).}    \tag{3.2}
\]

Therefore

\[
\kappa_4(p)\ge t_4(p).                                        \tag{3.3}
\]

For the converse, the straight horizontal path from `(0,0)` to `(n,0)` uses `n+1` sites, so

\[
P_p^{G4}(0\leftrightarrow ne_1)\ge p^{n+1},                   \tag{3.4}
\]

and

\[
\kappa_4(p)\le-\log p.                                        \tag{3.5}
\]

As `p downarrow0`,

\[
t_4(p)=-\log p+O(p),                                          \tag{3.6}
\]

which proves (1.1).

A more detailed Taylor expansion can be extracted from the explicit arcosh bound if useful, but the `O(p)` bracket is enough for the centre inversion.

## 4. Matching graph: a constructive forward-path lower probability

For the matching graph, a monotone-in-`x` step may change the vertical coordinate by `-1,0,+1`.  Starting at an occupied origin, run the following auxiliary randomized algorithm.

At column `x+1`, if the current height is `y`, inspect the three sites

\[
(x+1,y-1),\quad(x+1,y),\quad(x+1,y+1).
\]

If none is open, fail.  Otherwise choose one of the open sites uniformly and move there.

Different steps inspect different columns, so their occupation triples are independent.  Put

\[
q(p)=1-(1-p)^3=3p-3p^2+p^3.                                  \tag{4.1}
\]

At each step the survival probability is `q(p)`.  Conditional on survival, symmetry implies that the chosen vertical increment is exactly uniform on `{-1,0,+1}`: the three unconditional probabilities of choosing each increment are equal and sum to `q`.

Therefore, with `S_n` a sum of iid uniform `{-1,0,+1}` increments,

\[
P_{p,aux}(\text{algorithm reaches }(n,0))
=p\,q(p)^n P(S_n=0).                                          \tag{4.2}
\]

For every fixed occupied configuration, algorithmic success implies the existence of an open matching path.  Hence averaging over the auxiliary randomness gives the valid percolation lower bound

\[
P_p^{G8}(0\leftrightarrow ne_1)
\ge p\,q(p)^n P(S_n=0).                                      \tag{4.3}
\]

The lazy symmetric walk has

\[
P(S_n=0)=e^{o(n)}
\]

(in fact order `n^{-1/2}`).  Consequently

\[
\boxed{\kappa_8(p)\le-\log q(p).}                            \tag{4.4}
\]

This lower probability bound uses no independence between multiple candidate paths; the algorithm follows only one path and uses fresh columns.

## 5. Matching graph: all-walk upper probability

For king/matching steps

\[
S_8=\{(\pm1,0),(0,\pm1),(\pm1,\pm1)\}.
\]

The horizontal exponential moment is

\[
M_8(t)=6\cosh t+2.                                             \tag{5.1}
\]

Indeed the two horizontal axial steps contribute `2 cosh t`, the four diagonals contribute `4 cosh t`, and the two vertical steps contribute `2`.

Thus the all-walk bound gives

\[
\boxed{
\kappa_8(p)\ge
t_8(p):=\operatorname{arcosh}\left(\frac{1/p-2}{6}\right) }    \tag{5.2}
\]

for sufficiently small `p` so the displayed positive solution exists.

Now

\[
t_8(p)=-\log(3p)+O(p),                                       \tag{5.3}
\]

while

\[
-\log q(p)
=-\log(3p)+O(p).                                               \tag{5.4}
\]

Equations (4.4) and (5.2) prove (1.2).

The factor `3` has a transparent geometric meaning: to leading dilute order, each forward column offers three matching-compatible vertical increments, while NN has only the straight forward geodesic.  The proof above shows that this entropy survives the `n->infinity` limit at the logarithmic level; it is not inferred from a finite path count alone.

## 6. Inverting the masses

If

\[
d=-\log p+O(p),                                               \tag{6.1}
\]

then

\[
\log(pe^d)=O(p),
\]

so `p=e^{-d}(1+O(p))`; substituting once gives

\[
p=e^{-d}+O(e^{-2d}).                                          \tag{6.2}
\]

Applying this to (1.1) proves (1.3).

Likewise if

\[
d=-\log(3p)+O(p),                                             \tag{6.3}
\]

then

\[
p=\frac13e^{-d}+O(e^{-2d}),                                  \tag{6.4}
\]

proving (1.4).  Equations (1.5)--(1.6) follow algebraically.

## 7. Consequences for the dual-even/odd coordinates

The limiting separated-window coordinates are

\[
G_\infty(d)=b(d)-a(d),
\qquad
C_\infty(d)=\frac{a(d)+b(d)-1}{2}.
\]

Hence

\[
\boxed{G_\infty(d)=1-\frac43e^{-d}+O(e^{-2d}),}                \tag{7.1}
\]

\[
\boxed{C_\infty(d)=\frac13e^{-d}+O(e^{-2d}).}                  \tag{7.2}
\]

Using the exact integrated identities from `neutral-gas-topological-charge-20260914.md`,

\[
\int_0^1P_1(p)\,dp
\to 1-\frac43e^{-d}+O(e^{-2d}),                               \tag{7.3}
\]

and

\[
\boxed{
\int_0^1M(p)\,dp
\to-\frac23e^{-d}+O(e^{-2d}).}                                \tag{7.4}
\]

Thus even in the very elongated regime, where the unscaled birth mixture is nearly the endpoint split, the complement-odd area under the matching curve retains a calculable exponentially small asymmetry.

## 8. Relation to the fixed-p dilute component crossover

`dilute-winding-crossover.md` gives the more refined complete-component formula

\[
\nu_w(p)\sim p^w I_0(2wp)
\]

in the joint regime `wp^2->0`.  The present note has a different limit order: first define the infinite-plane mass by `n->infinity` at fixed `p`, then take `p->0`.  Its elementary `O(p)` mass bracket is consistent with the Bessel exponent `-log p-2p+...` on the NN side but does not assert the `-2p` coefficient at fixed `p`.

For matching, the factor `3` already appears at the shortest-path entropy level and is likewise compatible with the leading minimal matching-cycle coefficient in the existing dilute note.

No interchange of those two limits is made here.
