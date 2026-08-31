# The canonical joint response is one linear original-U source

**Result.** For the fixed homogeneous vacant tensor
`1+(epsilon/N)Kreg`, the complete requested response is

\[
 \boxed{J_2(N)=\left.\partial_{\log Q}\partial_\epsilon^2U
       \right|_{Q=1,\epsilon=0}=\mathcal L[s_2],\qquad
 s_2(A)=\frac{2}{N^2}\sum_{x<y}g_{xy}(A).}
 \tag{1}
\]

All original normalization, induced root motion and slope terms are
retained in L. No same-site term or covariance of separately closed
first-Q marks is to be added. Translation reduces the required moments
to an anchored source with normalization `1/N`, not `1/N²`.

A new configuration-level result is also available: the adjacent and
nonadjacent parts of s2 are linearly independent modulo
`span{1,K,a}`, where a is the old site-average one-insertion activation.
Thus neither their contact structure nor their joint source can be
universally removed by a density clock or a reparameterization of the
old one-insertion source. This does not predetermine either numerical
U response.

Base: `a237968f1d7a82d26b46e83c58179dbba7f1a908`.
The canonical completion is fixed in
[two-insertion algebra](local-pair-two-insertion-algebra.md), and g is
the unchanged [spatial kernel](regular-pair-spatial-kernel.md).
The four-term functional is the one already implemented in
`scripts/p337_local_pair_insertion_score.py::response`.
No old score, kernel table or population was recomputed.

## 1. The first-Q logarithm fixes the second-epsilon source

For fixed original occupation A, divide its perturbed colour weight
by its unperturbed weight and write

\[
 F_A(Q,\epsilon)=1+\frac\epsilon N\sum_x\beta_x(A,Q)
 +\frac{\epsilon^2}{N^2}\sum_{x<y}\beta_{xy}(A,Q)
 +\cdots .
 \tag{2}
\]

Each beta includes its vacant-mark indicators. All ranks, q and E
are computed on A, not on virtual projector joins. The finite-network
theorem gives `F_A(1,epsilon)=1` for every epsilon and every A:
every nonempty canonical insertion coefficient vanishes pointwise.
Consequently

\[
 \left.\partial_{\log Q}\log F_A\right|_1
 =\epsilon a(A)+\frac{\epsilon^2}{2}s_2(A)+O(\epsilon^3),
 \quad
 a=\frac1N\sum_x\left.\partial_Q\beta_x\right|_1,
 \quad g_{xy}=\left.\partial_Q\beta_{xy}\right|_1.
 \tag{3}
\]

Products of one-insertion betas have zero first-Q derivative at Q1;
their first nonzero possible order is `(Q-1)^2`. The factor two in
(1) is the second derivative of epsilon², not an additional ordered-
pair multiplicity. Equivalently
`s2=N^(-2) sum_(x!=y) g_xy` because the two-mark contraction is symmetric.
There is no x=y term: a site factor is linear in epsilon and cannot
supply two insertions at that site. An exponential local tensor would
be a different frozen model.

Write the complete unperturbed first-Q occupation score as S0(A).
It may include the declared background rank weight and activity-path
Jacobian. The first-Q score of the perturbed law is
`S0+epsilon a+epsilon² s2/2+...`. At Q1 its underlying normalized law
and thermal root are independent of epsilon. Differentiating twice
in epsilon therefore removes S0 and a; for any fixed observer O,

\[
 \left.\partial_{\log Q}\partial_\epsilon^2\langle O\rangle
       \right|_{1,0}=\operatorname{Cov}_{iid}(O,s_2).
 \tag{4}
\]

Background-Q, thermal and root terms involving a product of two
nonempty insertion responses vanish at this derivative order. This
does not remove the normalization or root response induced by s2
itself.

## 2. The four original-root terms remain

Normalize the two geometries separately. At the original pooled root
v0, define

\[
 M=\tfrac12(\langle q\rangle_f+\langle q\rangle_s),\quad
 Y=(\langle E\rangle_f-\langle E\rangle_s)/\Delta,
 \quad D=M_v,\ B=Y_v,\ T=M_{vv},\ H=Y_{vv}.
\]

For any source s use
`jM=(Cov_f(q,s_f)+Cov_s(q,s_s))/2` and
`jY=(Cov_f(E,s_f)-Cov_s(E,s_s))/Delta`. Then

\[
 \boxed{\frac{\mathcal L[s]}{A_N}
 =\frac{jY_v}{D}-\frac{H\,jM}{D^2}
       -\frac{B\,jM_v}{D^2}+\frac{B T\,jM}{D^3}.}
 \tag{5}
\]

The four pieces are direct centered response, root motion, slope-source
response and slope-root response. The induced mixed root displacement
is `v_(logQ,epsilon,epsilon)=-jM[s2]/D`. This proves (1), either by
differentiating the root equation or by applying the existing first-
variation functional to the first-Q score (3). The functional remains
linear in s; no new root is fitted for each pair or contact class.

## 3. Translation and a fixed contact split

At each geometry, g is translation-covariant and the iid law is
translation-invariant. For every translation-invariant f,

\[
 \boxed{\mathbb E[f s_2]
 =\mathbb E\!\left[f\,\frac1N\sum_{y\ne0}g_{0y}\right].}
 \tag{6}
\]

Indeed, all N anchored sums in the ordered-pair expression have the
same expectation after translating A and f. This is a moment identity,
not a pointwise replacement of s2 by one anchor. It applies to f=1,q,E,
K and their products, and holds for every v. It therefore supplies
all source moments, covariances and thermal derivatives in (5).

Define contact to mean a distinct NN pair and noncontact to mean every
other distinct pair. Use that same partition of pairs in each geometry:

\[
 s_2=s_c+s_n,\qquad
 \mathbb E[f s_c]=\mathbb E\!\left[f\,\frac1N
                 \sum_{y\sim0}g_{0y}\right],\qquad
 J_2=\mathcal L[s_c]+\mathcal L[s_n].
 \tag{7}
\]

On the N25 quotients this means four contact and twenty noncontact
displacements. Adjacent vacant marks share their actual isolated
edge-node; it must receive one common component label. Assigning two
different singleton labels would change the contact source. Both parts
come from the same occupation archive and retain their common
covariance/dependency group; they are not independent evidence blocks.
Here noncontact includes diagonally separated sites and does not mean
a macroscopic separation.

## 4. Contact and noncontact are not hidden copies of the old source

The following paper calculation works on both N25 quotients: their
periods introduce no new cycles of length four or less. All sites not
listed are vacant. In each row K=2, all occupied components are
contractible, and the site-average first activation is
`a=(N-2)/N`.

| occupied set | s_c | s_n | a |
|---|---:|---:|---:|
| nearest-neighbour domino `{(0,0),(1,0)}` | `1/N²` | 0 | `(N-2)/N` |
| plaquette diagonal `{(0,0),(1,1)}` | 0 | `1/(2N²)` | `(N-2)/N` |
| straight distance-two pair `{(0,0),(2,0)}` | 0 | 0 | `(N-2)/N` |

For the domino, precisely two adjacent vacant pairs, immediately above
and below it, share both the occupied component and their common
isolated edge-node. At each mark those two shared colours occupy
adjacent ports with two private colours. The known two-shared-component
contrast is -1/2 at each end, hence each pair has g=1/4. All other
pairs share at most one component. Two unordered pairs and the factor
`2/N²` give the first row.

For the diagonal occupation, the other two plaquette corners are the
unique vacant pair sharing both occupied singleton components. Their
shared ports are again adjacent at both ends, so g=1/4. The marks are
nonadjacent, giving the second row. For a straight distance-two pair,
there is only one common vacant neighbour; no pair of distinct marks
can share both components, giving the third row.

Finally, every vacant site in these three configurations has four
distinct exterior port components. For the domino no third site is
NN-adjacent to both occupied endpoints; in the other two cases two
occupied neighbours, when present, are different components. The
canonical single-insertion all-free activation is 1, while occupied
marks contribute zero. This proves the common a value in the table.

If `alpha s_c+beta s_n` were in `span{1,K,a}`, its value would be the
same on all three rows. The zero third row and the two nonzero rows
force alpha=beta=0. Thus the two joint sources are linearly independent
modulo that span. In fact any proposed function of only (K,a) fails
to distinguish these three configurations. This excludes an exact
occupation-level reduction to a density clock plus a reparameterized
old first-insertion source; it does not assert that the particular
linear functional L is nonzero on each remaining direction.

## 5. What the fixed global discriminator excludes

A model whose first-Q effective log weight is additive and linear in
the **fixed** epsilon, `Xi_A(epsilon)=Xi_A(0)+epsilon a(A)`, predicts
s2=0 and hence J2=0. Adding an epsilon² source proportional to 1 and
K still predicts zero after the common thermal quotient. Nonzero J2
therefore excludes this stated global closure. Conversely J2=0 would
not erase the already demonstrated conditional two-mark interaction:
occupation averaging and the original-U functional can cancel it.

Additivity without linearity in epsilon is a weaker, different claim.
Even one-dimensional parameter renaming `epsilon_old=epsilon_new+
c epsilon_new²+...` changes the displayed second response by `2c J1`.
Therefore a single nonzero J2 cannot exclude arbitrary nonlinear
renaming. The experiment fixes epsilon by the literal local tensor
`1+(epsilon/N)Kreg`; the three configurations above separately rule
out absorbing its joint occupation source into a common old-source
renaming plus a thermal score.

Finally, total J2 may be carried wholly by contact terms. A nonzero
noncontact part would still include short diagonal separations. Neither
number alone establishes a long-range field, a continuum operator
identity, or an additional independent source population.
