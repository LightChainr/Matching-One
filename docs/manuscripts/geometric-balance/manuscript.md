# Balance without concentration: two geometric scales for square-site percolation on integer-period tori

**Working mathematical manuscript — 2026-09-13.** This text consolidates the
probability argument in #735 and the axial full-law result in #736. Its new
step is an orientation-uniform staircase corridor and translation-packing
argument, which proves the full-law criterion for arbitrary integer periods.
It is an author-supplied proof under the explicit standard inputs below, not
an independently accepted publication or a claim of literature priority.

## Abstract

Let independent nearest-neighbour site percolation be defined on
\(\mathbb Z^2/\Lambda\), where \(\Lambda\) is a rank-two integer lattice, with
\(N=[\mathbb Z^2:\Lambda]\) vertices and shortest nonzero Euclidean period
\(\ell\). Write \(r\in\{0,1,2\}\) for the rank of the ambient homology image,
and \(P_j(p)=\Pr_p(r=j)\). The topological balance root solves \(P_2=P_0\).
Using site sharpness and matching duality, its convergence to the infinite
square-site critical probability is uniform over all period lattices as
\(\ell\to\infty\), with no area, aspect, or shear constraint. In contrast,
for any sequence of honest tori with \(N\to\infty\), the entire two-birth
mixture converges to a point mass at that critical probability **if and only
if** \(\log N/\ell\to0\). Necessity uses critical square-site box crossing:
an integer staircase of axis-aligned crossing rectangles closes around the
actual shortest period, at arbitrarily small exponential cost below
criticality. A finite-group packing lemma produces sufficiently many disjoint
translates without rotating the interaction or requiring an ambient-primitive
period. We also give an elementary, all-period endpoint-splitting corollary
when \(\log N/\ell\to\infty\). No critical exponent, conformal-invariance
assumption, numerical critical probability, or growing transfer matrix enters.

## 1. Model, inputs, and statements

### 1.1 The finite quantities

The physical interaction is the nearest-neighbour (NN) square lattice with
unit edges. Quotient by a full-rank subgroup \(\Lambda\le\mathbb Z^2\), while
retaining lifted edge displacements. Assume each periodic unit square is an
embedded cell with four distinct corners; this is automatic when the shortest
period is sufficiently large. These are the **honest tori** throughout. A
change of period basis is not a rotation of the NN interaction.

For an occupied vertex set \(\omega\), let \(G_\omega\) be its induced NN
graph in \(T_\Lambda=\mathbb R^2/\Lambda\), and put
\[
 r(\omega)=\dim_{\mathbb Q}\operatorname{im}
  [H_1(G_\omega;\mathbb Q)\longrightarrow H_1(T_\Lambda;\mathbb Q)].
\]
The source variable is \(X=r-1\), and
\[
 M_\Lambda(p)=\mathbb E_pX=P_2^\Lambda(p)-P_0^\Lambda(p),\qquad
 F_\Lambda(p)=\tfrac12\mathbb E_pr=\tfrac12(1+M_\Lambda(p)).             \tag{1}
\]
No directional wrapping marginal is substituted for rank. In particular a
rank-one spiral has nonzero projections on two coordinates but still rank one.

Give each vertex an independent uniform label in \([0,1]\), and occupy labels
at most \(p\). Write \(T_1,T_2\) for the first times at which rank is at least
one and two; a simultaneous jump is allowed. Then
\[
 F_\Lambda(p)=\tfrac12\Pr(T_1\le p)+\tfrac12\Pr(T_2\le p).              \tag{2}
\]
Thus \(F\) is a CDF, that of an independent fair choice between the two birth
times. It is not the rank law at a fixed parameter. Denote its inverse by
\(Q_\Lambda\), and its median, equivalently the balance root, by
\(p_\Lambda=Q_\Lambda(1/2)\). The conditional rank odds
\[
 H_\Lambda(p)=\frac{P_2^\Lambda(p)}{P_0^\Lambda(p)+P_2^\Lambda(p)}       \tag{3}
\]
are a different function.

### 1.2 Imported inputs and deterministic duality

Only the following infinite-volume inputs are used.

**S: site sharpness.** For NN and for its matching graph NN+NNN, each fixed
subcritical parameter has a bound
\[
 \Pr_p(0\leftrightarrow \partial B_R)\le C(p)e^{-c(p)R},\quad c(p)>0.
                                                                         \tag{4}
\]
Distances can be Euclidean, with a change in constants because both step sets
are finite. Duminil-Copin--Tassion [S], Theorem 1.1(3), is printed for bonds;
its section 1.2 explicitly discusses the site adaptation. We use the adapted
site conclusion, not the value or threshold theorem for square bond percolation.

**D: the matching critical relation.** With \(p_c=p_c^{\rm site}(\mathrm{NN})\),
\[
 p_c^{\rm site}(\mathrm{NN+NNN})=1-p_c.                                \tag{5}
\]
A direct source is Grimmett--Li [D], Theorem 5.5 in the amenable section.
Their Theorem 1.1(a) gives the more general relation with \(p_u\). We use the
amenable square-lattice instance, not a claim about every planar matching pair.

**R: critical square-site box crossing.** There is \(c_0>0\) such that for
every integer \(s\ge1\), the probability at \(p_c\) of an occupied horizontal
crossing of a \(3s\)-by-\(2s\) axis-aligned rectangle is at least \(c_0\).
Zeng [R], Theorem 1.1, states the requisite square-site RSW result. It is used
as a preprint theorem; no publication status beyond the cited text is assumed.
The full-law necessity uses R. Balance-root consistency and full-law
sufficiency use S and D but not R.

The finite identity
\[
 r_{\rm NN}(\omega)+r_{\rm NN+NNN}(\omega^c)=2                         \tag{6}
\]
is the deterministic digital-Alexander lemma used in the earlier repository
proofs. For completeness, its topological content is as follows. Take a closed
regular neighbourhood \(U\) of the black NN graph and its complementary
subsurface \(V\). The images \(A,C\) of their first homology in that of the
torus satisfy \(C=A^\perp\) for the intersection form. Indeed, the annihilator
of \(A\) identifies with the kernel of restriction \(H^1(T)\to H^1(U)\);
the relative exact sequence, excision, and Poincare--Lefschetz duality identify
that kernel with the image from \(H_1(V)\). Their dimensions sum to two.
The white matching graph has the same ambient homology image as \(V\): in a
face retain the diagonal only for an opposite-white-pair pattern; every other
active white diagonal has a white boundary replacement in that contractible
face. The retained local graph is an embedded spine up to filling local faces,
which does not change the ambient image. This is the precise 4/8 adjacency
convention. The existing `notes/digital-alexander-duality-proof.md` supplies
its facewise implementation. It is not inferred from the numerical controls
in this paper. Short-period cell degeneracies are outside this lemma.

### 1.3 Main conclusions

**Theorem A (balance consistency, consolidated from #735).** The balance
root is unique on every honest torus, and
\[
 \lim_{L\to\infty}\ \sup_{\Lambda:\ell(\Lambda)\ge L}
       |p_\Lambda-p_c|=0.                                             \tag{7}
\]
More precisely, for each fixed \(p<p_c\) there are \(L(p),\kappa(p)>0\),
independent of area and period shape, such that
\[
 \frac{P_2^\Lambda(p)}{P_0^\Lambda(p)}
 \le \exp[-\kappa(p)N/\ell],\qquad \ell\ge L(p).                        \tag{8}
\]
Above \(p_c\), the inverse ratio has the same type of estimate. Every fixed
interior quantile of \(H\) in (3) also converges uniformly under \(\ell\to\infty\).
This is a sufficient geometry for balance; necessity for balance is not asserted.

**Theorem B (sharp full-law geometry).** For any sequence of honest integer-
period tori with \(N_n\to\infty\), the following are equivalent:

1. \(\log N_n/\ell_n\to0\).
2. For every fixed \(p<p_c\), \(F_n(p)\to0\), and for every fixed
   \(p>p_c\), \(F_n(p)\to1\).
3. Each of \(T_{1,n},T_{2,n}\) converges in probability to \(p_c\).
4. The birth-time mixture in (2) converges weakly to \(\delta_{p_c}\).
5. \(Q_n(u)\to p_c\) for every fixed \(u\in(0,1)\).
6. The convergence in 5 is uniform on every compact subinterval of \((0,1)\).

The extension of necessity from axial periods (#736) to **all** integer
periods is proved in sections 4--6. It uses only axis-aligned NN rectangles,
not a rotated-lattice RSW assertion.

**Corollary C (extreme elongation, all period shapes).** If instead
\(\log N_n/\ell_n\to\infty\), then for each \(p\in(0,1)\),
\[
 P_0^n(p),P_2^n(p)\to0,\qquad F_n(p)\to\tfrac12,
\]
and the mixture converges to \(\tfrac12\delta_0+\tfrac12\delta_1\).
When also \(\ell_n\to\infty\), its finite median still converges to \(p_c\)
by Theorem A. This includes genuinely oblique examples, not only rectangles.

## 2. Finite strict monotonicity

Adding sites cannot decrease the ambient homology image. Both events
\(r\ge1\) and \(r\ge2\) are increasing and nonconstant. Along an empty-to-
full occupation chain each changes value somewhere. Such a pivotal assignment
of the other sites has strictly positive product probability at every
\(0<p<1\). The finite product derivative formula therefore gives
\[
 M_\Lambda'(p)=\frac{d}{dp}\Pr(r\ge1)+\frac{d}{dp}\Pr(r\ge2)>0.
\]
The endpoint values are \(-1,+1\), proving the unique root and invertibility
of \(F\). Also \(P_2'>0\) and \(P_0'<0\), so \(H'>0\) on the interior.
A simultaneous rank jump does not invalidate either argument.

## 3. Why balance needs no area restriction

This section puts the arbitrary-period argument of #735 into the same
manuscript; it is not a second independent proof certificate for that PR.

### 3.1 Shortest period and transverse height

Choose a shortest period \(u\), so \(|u|=\ell\). It is primitive in
\(\Lambda\), since a proper lattice multiple would not be shortest. Complete
it to a basis \((u,v)\) with determinant \(N>0\), and subtract an integer
multiple of \(u\) from \(v\) so that \(|u\cdot v|\le\ell^2/2\). Since
\(|v|\ge\ell\), the transverse height is
\[
 h=N/\ell\ge\sqrt3\ell/2.                                             \tag{9}
\]
Neither vector must be primitive in ambient \(\mathbb Z^2\). The circle
coordinate \(\theta(x)=\det(u,x)/\ell\pmod h\) is well defined. On vertices,
\(\det(u,x)\pmod N\) implements it exactly. A rank-two image contains a
closed walk with nonzero transverse winding. Every physical step changes
\(\theta\) by at most \(\sqrt2\), on either adjacency.

### 3.2 Overlapping local balls: a lower bound for no winding

Assume \(\ell\ge64\), put \(R=\ell/64\), and let \(a_R(p)\) be the one-arm
probability including the occupied origin. Stopping a path on its first exit
makes this an event on vertices at Euclidean distance at most \(R+\sqrt2\).
That ball injects into the torus because its diameter is less than \(\ell\).
The corresponding local event at each of the \(N\) vertices has exactly the
infinite-lattice probability \(a_R\).

A nonzero winding walk has a lift that escapes this radius. Consequently,
absence of every local arm implies rank zero. The arm-absence events overlap
but are decreasing. Harris positive association, not independence, gives
\[
 P_0\ge(1-a_R)^N.                                                       \tag{10}
\]
Harris association for a product measure follows by induction on the sites:
condition on one Bernoulli variable, apply induction to the conditional
covariances, and use that the two conditional means are monotone in the same
direction for the remaining covariance.

### 3.3 Disjoint transverse bands: an upper bound for rank two

Use \(k=\lfloor8N/\ell^2\rfloor\) transverse bands of physical width
\(\ell/8\), leaving the residual strip unused. Shift boundaries off vertices.
In each band require an occupied path using only its vertices, joining the
lower and upper boundary layers of thickness \(\sqrt2\).

A lift with nonzero transverse winding traverses every band. For a chosen
lifted band, take its first exit through the upper side after a visit below
the lower side, and the last entrance through the lower side preceding that
exit. The intervening path lies inside the band. This last-entry construction
allows arbitrary backtracking. The endpoint separation is at least
\(\ell/8-2\sqrt2>R\), and hence witnesses a local arm from an entry vertex.

The entry layer contains at most \(B=4\lceil\ell\rceil\) vertices. To see
this uniformly in tilt, center unit squares at the lattice vertices; they tile
the quotient with area one each. Squares centered in a layer of width
\(\sqrt2\) lie in a layer of width \(2\sqrt2\). The latter has area
\(2\sqrt2\ell\), since transverse fibres have length \(\ell\).
Thus a band crossing has probability at most \(Ba_R\). The bands have disjoint
site supports, so their events **are** independent. Therefore
\[
 P_2\le(Ba_R)^k,\qquad k\ge4N/\ell^2.                                 \tag{11}
\]
This remains valid for matching diagonals: edge lengths are at most
\(\sqrt2\), and edges between unused bands are not part of the events.

### 3.4 Rate comparison

At a fixed subcritical parameter, S gives \(a_R\le C e^{-c\ell/64}\).
For sufficiently large \(\ell\), independently of \(N\),
\[
 \log(Ba_R)\le-c\ell/128,\quad a_R\le1/2,\quad
 2a_R\le c/(64\ell).
\]
Equations (10)--(11) yield
\[
 \log(P_2/P_0)
 \le -cN/(32\ell)+2Na_R
 \le -cN/(64\ell).                                                     \tag{12}
\]
For \(p>p_c\), apply this to the matching graph at \(1-p\) and use (5)--(6).
Two fixed parameters \(p_c\pm\varepsilon\) trap the unique finite root
uniformly over all lattices with large \(\ell\). This proves A, including
the conditional-odds assertion. No value of the probabilities at \(p_c\)
and no quantitative control of \(c(p)\) near \(p_c\) is needed.

## 4. An occupied staircase ring around an arbitrary integer period

We now construct the new necessity mechanism. Constants are deliberately
conservative; their optimization is irrelevant to the criterion.

**Lemma 4.1 (deterministic staircase).** Given \(s\in\mathbb N\) and a
nonzero integer vector \(u\), one can specify at most \(8\ell/s\) horizontal
or vertical rectangle-crossing events, for \(\ell=|u|\ge64s\), such that:

- every rectangle has side lengths between \(2s\) and \(3s\);
- their projected intersection forces a closed occupied walk with lift
  displacement exactly \(u\), including the periodic seam;
- the lifted vertex support \(\widetilde S\) lies within distance \(4s\)
  of the segment \([0,u]\).

When \(u\) is a shortest period of \(\Lambda\), every individual rectangle
injects into \(T_\Lambda\), since its diameter is at most \(\sqrt{13}s<\ell\).

**Construction and proof.** A lattice reflection or quarter-turn is an exact
symmetry, so for notation arrange \(u=(a,b)\), \(a\ge b\ge0\). No arbitrary
angle rotation is performed. Set \(k=\lceil a/s\rceil\), and
\[
 z_j=(\lfloor ja/k\rfloor,\lfloor jb/k\rfloor),\quad 0\le j\le k.
                                                                         \tag{13}
\]
Between \(z_j\) and \(z_{j+1}\), move horizontally to
\((z_{j+1,x},z_{j,y})\), then vertically to \(z_{j+1}\). Omit zero moves.
The resulting centers \(c_0,\ldots,c_J\) satisfy
\(c_0=0,c_J=u\), \(J\le2k\), and each axis step has length at most \(s\).

At every \(c_i\), \(0\le i<J\), require both a horizontal and a vertical
crossing of the square \(Q_i=c_i+[-s,s]^2\). For a horizontal step from
\((x,y)\) to \((x+d,y)\), require a horizontal crossing of
\([x-s,x+d+s]\times[y-s,y+s]\). For a vertical step use the corresponding
vertical rectangle. There are \(3J\le6\lceil a/s\rceil\le8\ell/s\) events.
The last connector reaches \(Q_J=Q_0+u\); its crossing conditions are those
of \(Q_0\), lifted periodically.

The two crossings in each hub meet, and form a connected set \(C_i\). A
horizontal connector meets the vertical crossing of each incident hub: its
restriction to that hub contains a full left-right subcrossing. The vertical
case is analogous. These are intersections of actual planar NN paths, hence
common vertices. Choose hub paths once, with \(C_J=C_0+u\). The connected
chain joins any chosen \(z\in C_0\) to \(z+u\in C_J\), completing an
occupied nonzero-winding closed walk on the torus. Merely reaching opposite
cut sides would not establish this; the periodic last hub does.

Each \(z_j\) is within \(\sqrt2\) of the straight segment. Each point of
a staircase step is within an additional \(s\) of an endpoint, and each
connector rectangle is the step thickened by \([-s,s]^2\). The distance to
\([0,u]\) is therefore at most
\((1+\sqrt2)s+\sqrt2<4s\). This proves the support assertion. The construction
never uses \(\gcd(a,b)=1\). In particular, a shortest period nonprimitive in
ambient \(\mathbb Z^2\) is treated exactly as a primitive one. \(\square\)

**Lemma 4.2 (arbitrarily small subcritical cost).** For every \(\eta>0\)
there are a fixed integer \(s=s_\eta\) and a fixed parameter
\(p_\eta\in(0,p_c)\), such that for every shortest integer period
\(u\) with \(\ell\ge64s\), the event \(\mathcal G_{u,s}\) of Lemma 4.1
satisfies
\[
 \Pr_{p_\eta}(\mathcal G_{u,s})\ge e^{-\eta\ell}.                        \tag{14}
\]

**Proof.** Choose \(0<c_0<1\) as in R and put \(a_0=c_0/2\).
A crossing of a \(3s\)-by-\(2s\) rectangle implies the required same-
direction crossing of any shorter connector or square by path restriction.
Reflections and quarter-turns handle the other direction. For each **fixed**
\(s\), finite-product continuity provides one \(p_s<p_c\) at which this
finite-box probability is at least \(a_0\). All translated rectangles then
have the same lower bound, independent of the torus, orientation and period
arithmetic. Their individual injection is important here.

Their supports overlap, so apply Harris rather than independence. Lemma 4.1
gives a probability at least \(a_0^{8\ell/s}\). Choose
\(s\ge8\log(1/a_0)/\eta\), then choose \(p_\eta=p_s\). This proves (14).
The order is \(\eta\), then \(s\), then \(p_\eta\), then the torus limit.
Neither \(s\) nor \(p_\eta\) is allowed to vary with that limit. \(\square\)

## 5. Many independent attempts without coordinate-dependent slicing

**Lemma 5.1 (finite-group translation packing).** If \(S\) is a nonempty
subset of a finite abelian group \(\Gamma\), there exist at least
\(|\Gamma|/|S-S|\) pairwise disjoint translates of \(S\).

**Proof.** Choose a maximal family \(c_i+S\). Any other possible center must
lie in some \(c_i+(S-S)\), or it could be added. These forbidden-center sets
cover \(\Gamma\), so \(|\Gamma|\le k|S-S|\). \(\square\)

This elementary packing statement does not assert that arbitrary maximal
packings are optimal. Only the lower bound is used.

**Lemma 5.2 (small difference support).** The projected support \(S\) of
the staircase ring satisfies
\[
 |S-S|\le64s\ell,\qquad \ell\ge64s.                                    \tag{15}
\]

**Proof.** From Lemma 4.1,
\(\widetilde S-\widetilde S\subset[-u,u]+B_{8s}\).
Center disjoint unit squares at the integer points in this capsule. Their union
is contained in \([-u,u]+B_{8s+1}\). Its area is
\[
 4\ell(8s+1)+\pi(8s+1)^2
 \le36s\ell+324s^2<64s\ell.
\]
Here \(s\ge1,\ell\ge64s\) and \(\pi<4\) suffice. Passing to the quotient
cannot increase cardinality. This proves (15) without a Bezout coordinate,
a primitive ambient vector, or a preferred direction of the transverse seam.
\(\square\)

Combining Lemmas 5.1--5.2 gives at least \(N/(64s\ell)\) disjoint translated
site supports on the actual finite quotient group \(\mathbb Z^2/\Lambda\).
Events on these supports are independent under the product site measure;
unused edges between supports do not change their measurability. Since any
successful ring forces positive rank, Lemma 4.2 yields the finite bound
\[
 \boxed{\displaystyle
 P_0^\Lambda(p_\eta)
 \le\exp\left[-\frac{N}{64s_\eta\ell}e^{-\eta\ell}\right],
 \qquad\ell\ge64s_\eta.}                                             \tag{16}
\]
This is the orientation-uniform winding-corridor estimate missing from #736.
It bounds \(P_0\) from above; one successful ring does not force rank two.

## 6. Proof of the sharp full-law criterion

### 6.1 Necessity

Suppose \(\log N_n/\ell_n\not\to0\). On a subsequence,
\(\log N_n\ge d\ell_n\) for a fixed \(d>0\).

First suppose a further subsequence has \(\ell_n\to\infty\). Choose
\(\eta=d/2\) in (16). Its exponent in absolute value is at least
\[
 \frac{e^{(d/2)\ell_n}}{64s_\eta\ell_n}\longrightarrow\infty.
\]
Thus \(P_0^n(p_\eta)\to0\) at a **fixed** \(p_\eta<p_c\), whence
\[
 \liminf_n F_n(p_\eta)\ge\tfrac12.
\]
For every fixed \(u<1/2\), eventually \(Q_n(u)\le p_\eta\) on this
subsequence. All-quantile convergence is impossible. This argument deliberately
makes no assertion that the median fails.

If instead \(\ell_n\) remains bounded on a further subsequence, a fixed
bounded-length occupied NN closed walk suffices. Choose a shortest period
\(u_n\) and a monotone axis-step path from \(0\) to \(u_n\). Its projected
support \(S_n\) has at most \(|u_n|_1\le\sqrt2\ell_n\) vertices. Requiring
all these sites occupied produces nonzero winding with probability at least
\(p^{\sqrt2\ell_n}\). Lemma 5.1, with \(|S_n-S_n|\le2\ell_n^2\), gives
\[
 P_0^n(p)\le
 \exp\left[-\frac{N_n}{2\ell_n^2}p^{\sqrt2\ell_n}\right].              \tag{17}
\]
At any fixed \(0<p<p_c\), this tends to zero when \(\ell_n\) is bounded
and \(N_n\to\infty\). The same lower-quantile obstruction follows.
Every witnessing subsequence has one of these further subsequences. This
proves that assertion 5 of B implies assertion 1.

### 6.2 Sufficiency and equivalences

If \(\log N_n/\ell_n\to0\), then \(\ell_n\to\infty\). At a fixed
\(p<p_c\), an occupied nonzero-winding cycle must produce a local arm of
radius \(\ell_n/64\) somewhere. The union bound and S give
\[
 1-P_0^n(p)\le C(p)N_n e^{-c(p)\ell_n/64}\longrightarrow0.             \tag{18}
\]
For \(p>p_c\), apply (18) to the matching complement using D and (6).
Thus \(P_2^n(p)\to1\), proving assertion 2. This sufficient argument is
already present in #613; it is not new necessity disguised as a union bound.

Monotonicity of \(F_n\) traps all \(Q_n(u)\), \(u\in[\delta,1-\delta]\),
between \(p_c-\varepsilon\) and \(p_c+\varepsilon\) for large \(n\).
This proves 2 => 6 => 5. Assertion 4 is equivalent to 2 by the CDF criterion
for weak convergence to a point mass. Finally, because \(T_1\le T_2\),
\[
 \Pr(T_1\le p)\le2F_n(p),\qquad
 \Pr(T_2>p)\le2[1-F_n(p)],
\]
so 2 implies 3. Conversely 3 and (2) imply 4. All assertions in B follow.

### 6.3 Extreme elongation

Inequality (17) does not require RSW. It is valid for every honest period
lattice. Apply it also to the matching graph, which contains the same forced
NN path, at \(1-p\). Duality gives
\[
 P_2^\Lambda(p)\le
 \exp\left[-\frac{N}{2\ell^2}(1-p)^{\sqrt2\ell}\right].                \tag{19}
\]
If \(\log N/\ell\to\infty\), both exponents in (17) and (19) diverge for
every fixed \(0<p<1\). This proves Corollary C. It extends the full/empty-row
control to arbitrary period shapes using the same elementary packing lemma.

For example, set \(u_n=(3n,4n)\), \(v_n=k_n(-4n,3n)\). Then
\(\ell_n=5n\), \(N_n=25n^2k_n\), and the shortest vector is nonprimitive
in ambient \(\mathbb Z^2\) when \(n>1\). With \(k_n=\lceil e^{n^2}\rceil\),
the birth mixture splits to the endpoints while the balance root converges
to \(p_c\). With \(k_n=\lceil e^{\sqrt n}\rceil\), all fixed quantiles
converge to \(p_c\). With \(k_n=\lceil e^{dn}\rceil\), \(d>0\), Theorem B
rules out full-law concentration but does not specify its limiting nonmedian
quantiles. No thin-continuum scaling law is used for these lattice sequences.

## 7. Prior work, claim boundaries, and what remains

The ambient homology observable, matching-function root, positive association,
RSW gluing, and translation packing are not presented as inventions of this
manuscript. The following comparison separates source scope from our deductions.

| Source | Scope inspected in the primary text | Relation to this manuscript |
|---|---|---|
| Mertens--Ziff [MZ], equations (20)--(21) and the subsequent root paragraph | Finite matching lattices, wrapping contrasts, and convergence of roots on square sequences | Establishes the existing observable/root context; not an arbitrary-period full-law necessity statement |
| Duncan--Kahle--Schweinhart [DKS], section 1.1, Theorems 1--4 | Cubical plaquettes on uniform \(N\mathbb Z^d\) quotients and scaled permutohedral systems; in 2D the basic examples are bond square/site triangular | Same ambient-image question, but different finite sequence/model hypotheses; no automatic extension to arbitrary square-site period lattices |
| Zeng [R], Theorem 1.1 | Uniform critical square-site rectangle crossings | Imported box lower bound; the staircase and its endpoint/seam argument are supplied here |
| Duminil-Copin--Tassion [S], Theorem 1.1(3), section 1.2 | Subcritical sharpness, with explicit site adaptation | Imported one-arm decay, not a computed near-critical modulus |
| Grimmett--Li [D], Theorem 5.5; Theorem 1.1 and Remark 1.4 | Site matching critical relation in the amenable case and its wider uniqueness-threshold form | Imported infinite-graph relation; finite rank duality is a separate topological statement |
| Repository #735 and #736 | Arbitrary-period root argument; axial full-law necessity | Author-level predecessors consolidated here; sections 4--6 remove the axial restriction |

This bounded source comparison does not certify novelty. It did not include
systematic citation-graph traversal, books, theses, or all strip-percolation
literature. In particular, earlier results on elongated strips may package
related rare-opportunity arguments in different language. No claim of absence
from print follows from the searches made for this delivery.

The manuscript does not solve the original-U candidate-map problem, identify
an irrelevant field, prove an \(L^{-4}\) shift, furnish a new numerical bound
for \(p_c\), establish fixed-width-to-continuum interchange, or predict the
nonmedian limiting law at finite positive \(\log N/\ell\). Theorem A does
not assert its systole condition is necessary for the median. Theorem B's
necessity is for the entire law, not a failure of every selected quantile.

The real next mathematical comparison is the consolidated theorem versus its
closest existing strip and homological statements. An independent reader can
challenge the seam construction, the finite-group support packing, or the
imported site inputs directly in this one text. No further width table,
source-jet hierarchy, or unrelated Jordan example is a dependency.

## 8. Executable control and integration boundary

`scripts/oblique_winding_corridor.py` implements only the new finite geometry:
integer quotient coordinates, staircase rectangles, planar rectangle-crossing
BFS, and a separate spanning-forest winding detector on the physical NN graph.
The tiny exhaustive controls use injecting rectangles, not the asymptotic
\(\ell\ge64s\) constants. They test the deterministic gluing rather than RSW.
The larger full-support cases check actual reduced oblique bases, including
ambient-nonprimitive shortest vectors; rational distances verify the tube
bound. The packing check uses independent finite cyclic-group examples.

Execution details are in
`results/research-control-20260913/oblique-corridor-controls.json`.
Five local mathematical tests passed. The three tiny tori comprise 135,168
configurations with no ring-without-winding failure. These checks do not prove
the all-size probability inputs or supply independent publication review.
No Monte Carlo, paid compute, GPU, or full Matching-One repository CI was run.
The code depends only on the standard library, not on the unmerged width-four
certificate stack. Old proofs and frozen results are not overwritten.

## References

[S] H. Duminil-Copin and V. Tassion, *A new proof of the sharpness of the phase
transition for Bernoulli percolation and the Ising model*. arXiv:1502.03050v3,
Theorem 1.1 and section 1.2. Primary HTML read:
https://arxiv.org/html/1502.03050v3

[D] G. Grimmett and Z. Li, *Hyperbolic site percolation*. arXiv:2203.00981,
Theorem 5.5 (amenable matching pairs), also Theorem 1.1 and Remark 1.4.
Primary theorem text read: https://arxiv.org/html/2203.00981
Their *Percolation critical probabilities of matching lattice-pairs*,
arXiv:2205.02734v3, introduction (1.3), gives the companion context:
https://arxiv.org/html/2205.02734v3

[R] X. Zeng, *A Russo Seymour Welsh Theorem for critical site percolation on
\(\mathbb Z^2\)*, arXiv:1309.2273v1, Theorem 1.1. Primary theorem and setup
read: https://arxiv.org/html/1309.2273
Kohler-Schindler--Tassion, *Crossing probabilities for planar percolation*,
arXiv:2011.04618, is general RSW background, not a replacement for checking
the critical square-site input: https://arxiv.org/html/2011.04618

[MZ] S. Mertens and R. M. Ziff, *Percolation in Finite Matching Lattices*,
arXiv:1603.07289v2. Primary HTML definitions and root discussion read:
https://arxiv.org/html/1603.07289v2

[DKS] P. Duncan, M. Kahle and B. Schweinhart, *Homological percolation on a
torus: plaquettes and permutohedra*, arXiv:2011.11903v4. Primary section 1.1
and theorem statements read: https://arxiv.org/html/2011.11903v4

Repository predecessors (not independent sources of validation):
#735 at `9d29d014df28af7c635e6859d98a95ffe2b34d06`;
#736 at `64d809b4404f80ff3f9adf9713337cc76008e92d`;
the finite digital-Alexander note on main. Navigation reset #738 was read
on main together with its subsequent focus clarification.
